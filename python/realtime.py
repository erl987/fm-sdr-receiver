#  FM-radio - software defined radio using RTL-SDR
#  Copyright (C) 2019-2020 Ralf Rettig
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

import asyncio
import os
import queue

from rtlsdr import *
import sounddevice as sd
import numpy as np
from scipy import signal
import psutil

from compiled_utils import _prepare_downsample, _perform_dropout, _baseband_mixer_impl, _polar_discriminator_impl, \
    _limit_impl

STATION_FREQ = 103.3e6  # in Hz (Radio P4 Stockholm)
DESIRED_AUDIO_SAMPLE_RATE = 48.0e3  # in Hz

REL_DISTANCE_TO_CENTER_FREQ = 0.2
FM_BANDWIDTH = 200.0e3  # in Hz
DEEMPHASIS_DECAY_TIME = 50.0e-6  # in s
SDR_SAMPLE_RATE = 2.304e6  # in Hz
DESIRED_SAMPLE_PERIOD = 0.0133  # in s

DC_OFFSET_FILTER_ORDER = 3
DC_OFFSET_CUTOFF_FREQ = 300  # in Hz
DOWNSAMPLE_TO_CHANNEL_FILTER_ORDER = 3
DOWNSAMPLE_TO_AUDIO_FILTER_ORDER = 3

LIMITER_THRESHOLD = 0.9
LIMITER_DELAY = 40
LIMITER_RELEASE_COEFF = 0.9999
LIMITER_ATTACK_COEFF = 0.9


class ContinuousDownsampler(object):
    def __init__(self, dec_rate, sample_rate, output_size, filter_order, dtype):
        if dec_rate > 1:
            b, a = ContinuousDownsampler._design_filter(dec_rate, sample_rate, filter_order)
            self._do_downsample = True
            self._b = b
            self._a = a
            self._z = signal.lfilter_zi(b, a)
            self._dec_rate = dec_rate
            self._output_size = output_size
            self._num_processed_samples = self._dec_rate * self._output_size
            self._unprocessed_samples_buffer = np.empty(0, dtype=dtype)
        else:
            self._do_downsample = False

    @staticmethod
    def _design_filter(dec_rate, sample_rate, filter_order):
        b, a = signal.iirfilter(filter_order, 1 / (2 * dec_rate) * sample_rate, btype="lowpass", fs=sample_rate)
        return b, a

    async def downsample(self, samples):
        if not self._do_downsample:
            return samples, False

        samples_for_processing, self._unprocessed_samples_buffer = \
            _prepare_downsample(samples, self._unprocessed_samples_buffer, self._num_processed_samples)

        filtered_samples, self._z = signal.lfilter(self._b, self._a, samples_for_processing, zi=self._z)
        downsampled_samples = _perform_dropout(filtered_samples, self._dec_rate)

        if len(self._unprocessed_samples_buffer) >= self._num_processed_samples:
            return downsampled_samples, True
        else:
            return downsampled_samples, False


class ContinuousFilter(object):
    def __init__(self, b, a):
        self._b = b
        self._a = a
        self._z = signal.lfilter_zi(b, a)

    async def filter(self, samples):
        filtered_samples, self._z = signal.lfilter(self._b, self._a, samples, zi=self._z)
        return filtered_samples


class DeemphasisFilter(ContinuousFilter):
    def __init__(self, sample_rate, deemphasis_decay_time):
        b, a = DeemphasisFilter._design_filter(sample_rate, deemphasis_decay_time)
        super().__init__(b, a)

    @staticmethod
    def _design_filter(sample_rate, deemphasis_decay_time):
        d = sample_rate * deemphasis_decay_time  # calculate the number of samples to hit the -3dB point
        x = np.exp(-1 / d)  # calculate the decay between each sample
        b = [1 - x]
        a = [1, -x]
        return b, a


class DCOffsetFilter(ContinuousFilter):
    def __init__(self, filter_order, cutoff_freq, sample_rate):
        b, a = DCOffsetFilter._design_filter(filter_order, cutoff_freq, sample_rate)
        super().__init__(b, a)

    @staticmethod
    def _design_filter(filter_order, cutoff_freq, sample_rate):
        b, a = signal.iirfilter(filter_order, cutoff_freq, btype="highpass", fs=sample_rate)
        return b, a


class Limiter(object):
    def __init__(self, attack_coeff, release_coeff, delay, dtype=np.float32):
        self._delay_index = 0
        self._envelope = 0
        self._gain = 1
        self._delay = delay
        self._delay_line = np.zeros(delay, dtype=dtype)
        self._release_coeff = release_coeff
        self._attack_coeff = attack_coeff

    async def limit(self, signal, threshold):
        limited_signal, self._delay_line, self._delay_index, self._delay, self._envelope, self._gain, \
        self._attack_coeff, self._release_coeff = \
            _limit_impl(signal, threshold, self._delay_line, self._delay_index, self._delay, self._envelope,
                        self._gain, self._attack_coeff, self._release_coeff)

        return limited_signal


class SdrFmReceiver(object):
    SAMPLE_FACTOR = 16384

    def __init__(self, station_freq, rel_distance_to_center_freq, fm_bandwidth, sdr_sample_rate,
                 desired_audio_sample_rate, deemphasis_decay_time, downsample_to_channel_filter_order,
                 downsample_to_audio_filter_order, dc_offset_filter_order, dc_offset_cutoff_freq, desired_sample_period,
                 attack_coeff, release_coeff, delay, threshold):

        dec_rate, dec_rate_audio, channel_sample_rate = \
            self._configure_frequency_settings(station_freq, rel_distance_to_center_freq, fm_bandwidth, sdr_sample_rate,
                                               desired_audio_sample_rate, desired_sample_period)
        self._configure_filters(sdr_sample_rate, dec_rate, downsample_to_channel_filter_order, dec_rate_audio,
                                downsample_to_audio_filter_order, channel_sample_rate, deemphasis_decay_time,
                                dc_offset_filter_order, dc_offset_cutoff_freq, attack_coeff, release_coeff, delay)

        self._last_bandwidth_sample = 0
        self._last_angle = 0
        self._do_receive = True
        self._audio_samples_queue = asyncio.Queue()
        self._limiter_threshold = threshold

    def _configure_frequency_settings(self, station_freq, rel_distance_to_center_freq, fm_bandwidth, sdr_sample_rate,
                                      desired_audio_sample_rate, desired_sample_period):
        self._sdr_sample_rate = sdr_sample_rate
        self._sdr_num_samples = int(SdrFmReceiver.SAMPLE_FACTOR * np.ceil(2 * desired_sample_period * sdr_sample_rate /
                                                                          SdrFmReceiver.SAMPLE_FACTOR))
        self._center_freq = station_freq + rel_distance_to_center_freq * sdr_sample_rate
        self._offset_freq = station_freq - self._center_freq

        dec_rate = int(sdr_sample_rate / fm_bandwidth)
        self._channel_blocksize = int(self._sdr_num_samples / dec_rate)
        channel_sample_rate = sdr_sample_rate / dec_rate
        self._delta_phase = 2.0 * np.pi * self._offset_freq / sdr_sample_rate

        dec_rate_audio = int(channel_sample_rate / desired_audio_sample_rate)
        self._audio_sample_rate = channel_sample_rate / dec_rate_audio
        self._audio_blocksize = int(self._channel_blocksize / dec_rate_audio)

        return dec_rate, dec_rate_audio, channel_sample_rate

    def _configure_filters(self, sdr_sample_rate, dec_rate, downsample_to_channel_filter_order, dec_rate_audio,
                           downsample_to_audio_filter_order, channel_sample_rate, deemphasis_decay_time,
                           dc_offset_filter_order, dc_offset_cutoff_freq, attack_coeff, release_coeff, delay):

        self._downsampler_to_channel = ContinuousDownsampler(dec_rate, sdr_sample_rate, self._channel_blocksize,
                                                             downsample_to_channel_filter_order, dtype=complex)
        self._downsampler_to_audio = ContinuousDownsampler(dec_rate_audio, channel_sample_rate, self._audio_blocksize,
                                                           downsample_to_audio_filter_order, dtype=float)
        self._de_emphasis = DeemphasisFilter(channel_sample_rate, deemphasis_decay_time)
        self._dc_offset = DCOffsetFilter(dc_offset_filter_order, dc_offset_cutoff_freq, self._audio_sample_rate)
        self._limiter = Limiter(attack_coeff, release_coeff, delay)

    def get_audio_device_config(self):
        return self._audio_sample_rate, self._audio_blocksize

    async def stream_to_audio(self):
        sdr = RtlSdr()
        sdr.sample_rate = self._sdr_sample_rate
        sdr.center_freq = self._center_freq

        try:
            async for iq_samples in sdr.stream(num_samples_or_bytes=self._sdr_num_samples):
                if not self._do_receive:
                    break
                audio_samples, first_has_still_data, second_has_still_data = await self._demodulate_to_audio(iq_samples)
                self._audio_samples_queue.put_nowait(audio_samples)

                if first_has_still_data or second_has_still_data:
                    # process the still buffered samples
                    if first_has_still_data:
                        audio_samples, __, __ = \
                            await self._perform_baseband_processing(np.empty(0, dtype=complex))
                    else:
                        audio_samples, __ = \
                            await self._perform_channel_processing(np.empty(0, dtype=float))
                    self._audio_samples_queue.put_nowait(audio_samples)
        finally:
            await sdr.stop()
            sdr.close()

    async def get_next_audio_batch(self):
        return await self._audio_samples_queue.get()

    def stop(self):
        self._do_receive = False

    async def _baseband_mixer(self, iq_samples):
        baseband_samples, self._last_angle = _baseband_mixer_impl(iq_samples, self._last_angle, self._delta_phase,
                                                                  self._offset_freq, self._sdr_sample_rate)
        return baseband_samples

    async def _polar_discriminator(self, samples_in_bandwidth):
        channel_samples, self._last_bandwidth_sample = _polar_discriminator_impl(samples_in_bandwidth,
                                                                                 self._last_bandwidth_sample)
        return channel_samples

    async def _demodulate_to_audio(self, iq_samples):
        baseband_samples = await self._baseband_mixer(iq_samples)
        audio_samples, first_has_more, second_has_more = await self._perform_baseband_processing(baseband_samples)
        return audio_samples, first_has_more, second_has_more

    async def _perform_baseband_processing(self, baseband_samples):
        samples_in_bandwidth, first_has_more = await self._downsampler_to_channel.downsample(baseband_samples)
        channel_samples = await self._polar_discriminator(samples_in_bandwidth)

        deemphasised_channel_samples = await self._de_emphasis.filter(channel_samples)

        audio_samples, second_has_more = await self._perform_channel_processing(deemphasised_channel_samples)
        return audio_samples, first_has_more, second_has_more

    async def _perform_channel_processing(self, deemphasised_channel_samples):
        audio_samples_with_dc, has_more = await self._downsampler_to_audio.downsample(deemphasised_channel_samples)
        unlimited_audio_samples = await self._dc_offset.filter(audio_samples_with_dc)
        audio_samples = await self._limiter.limit(unlimited_audio_samples, self._limiter_threshold)

        audio_samples = audio_samples.reshape((len(audio_samples), 1))
        return audio_samples, has_more


class AudioDevice(object):
    def __init__(self, sample_rate, block_size, dtype="float32", channels=1):
        self._stream = sd.OutputStream(samplerate=sample_rate, blocksize=block_size, callback=self._audio_callback,
                                       dtype=dtype, channels=channels, latency="low")
        self._audio_samples_queue = queue.Queue()
        self._do_play_stream = True

    def _audio_callback(self, outdata, frames, time, status):
        if not self._audio_samples_queue.empty():
            outdata[:] = self._audio_samples_queue.get_nowait()

    async def play_stream(self, stream_provider):
        with self._stream:
            while self._do_play_stream:
                audio_samples = await stream_provider.get_next_audio_batch()
                self._audio_samples_queue.put_nowait(audio_samples)

    def stop(self):
        self._do_play_stream = False


async def main():
    fm_radio = SdrFmReceiver(STATION_FREQ, REL_DISTANCE_TO_CENTER_FREQ, FM_BANDWIDTH, SDR_SAMPLE_RATE,
                             DESIRED_AUDIO_SAMPLE_RATE, DEEMPHASIS_DECAY_TIME, DOWNSAMPLE_TO_CHANNEL_FILTER_ORDER,
                             DOWNSAMPLE_TO_AUDIO_FILTER_ORDER, DC_OFFSET_FILTER_ORDER, DC_OFFSET_CUTOFF_FREQ,
                             DESIRED_SAMPLE_PERIOD, LIMITER_ATTACK_COEFF, LIMITER_RELEASE_COEFF, LIMITER_DELAY,
                             LIMITER_THRESHOLD)
    audio_device = AudioDevice(*fm_radio.get_audio_device_config())

    sdr_task = asyncio.create_task(fm_radio.stream_to_audio())
    audio_task = asyncio.create_task(audio_device.play_stream(fm_radio))

    await sdr_task
    await audio_task


if __name__ == '__main__':
    p = psutil.Process(os.getpid())
    p.nice(psutil.HIGH_PRIORITY_CLASS)

    asyncio.run(main())
