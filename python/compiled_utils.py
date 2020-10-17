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

import numpy as np
from numba import jit


@jit(nopython=True, nogil=True, fastmath=True, cache=True)
def _prepare_downsample(samples, unprocessed_samples_buffer, num_processed_samples):
    samples_for_processing = np.empty(len(unprocessed_samples_buffer) + len(samples), dtype=samples.dtype)
    samples_for_processing[0:len(unprocessed_samples_buffer)] = unprocessed_samples_buffer
    samples_for_processing[len(unprocessed_samples_buffer):] = samples

    unprocessed_samples_buffer = samples_for_processing[num_processed_samples:]
    samples_for_processing = samples_for_processing[0:num_processed_samples]

    return samples_for_processing, unprocessed_samples_buffer


@jit(nopython=True, nogil=True, fastmath=True, cache=True)
def _perform_dropout(filtered_samples, dec_rate):
    return filtered_samples[0::dec_rate]


@jit(nopython=True, nogil=True, fastmath=True, cache=True)
def _polar_discriminator_impl(samples_in_bandwidth, last_bandwidth_sample):
    channel_samples = np.empty(len(samples_in_bandwidth))
    if len(samples_in_bandwidth) > 0:
        first_channel_sample = np.angle(samples_in_bandwidth[0] * np.conj(last_bandwidth_sample))
        channel_samples[0] = first_channel_sample
        channel_samples[1:] = np.angle(samples_in_bandwidth[1:] * np.conj(samples_in_bandwidth[:-1]))
        last_bandwidth_sample = samples_in_bandwidth[-1]

    return channel_samples, last_bandwidth_sample


@jit(nopython=True, nogil=True, fastmath=True, cache=True)
def _baseband_mixer_impl(iq_samples, last_angle, delta_phase, offset_freq, sdr_sample_rate):
    factor = np.exp(-1.0j * (last_angle + delta_phase)) * \
             np.exp(-1.0j * 2.0 * np.pi * offset_freq / sdr_sample_rate * np.arange(len(iq_samples)))
    baseband_samples = iq_samples * factor
    last_angle = 2 * np.pi - np.angle(factor[-1])

    return baseband_samples, last_angle


@jit(nopython=True, nogil=True, fastmath=True, cache=True)
def _limit_impl(signal, threshold, delay_line, delay_index, delay, envelope, gain, attack_coeff, release_coeff):
    limited_signal = np.empty(len(signal))

    for i in np.arange(len(signal)):
        delay_line[delay_index] = signal[i]
        delay_index = (delay_index + 1) % delay

        # calculate an envelope of the signal
        envelope *= release_coeff
        envelope = max(abs(signal[i]), envelope)

        # have gain go towards a desired limiter gain
        if envelope > threshold:
            target_gain = 1 + threshold - envelope
        else:
            target_gain = 1.0
        gain = gain * attack_coeff + target_gain * (1 - attack_coeff)

        # limit the delayed signal
        limited_signal[i] = delay_line[delay_index] * gain

    return limited_signal, delay_line, delay_index, delay, envelope, gain, attack_coeff, release_coeff
