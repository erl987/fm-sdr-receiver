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

# this script is based on: https://witestlab.poly.edu/blog/capture-and-decode-fm-radio/
import time

from matplotlib import pyplot as plt
from rtlsdr import *
import sounddevice as sd
import numpy as np
from scipy import signal

stations = {"Radio P4 Stockholm": 103.3e6}  # in Hz
center_freq = 104.0e6  # in Hz
t_length = 15.0  # in s
fm_bandwidth = 200.0e3  # in Hz
sample_rate = 2.4e6  # in Hz
desired_audio_sample_rate = 44.1e3  # in Hz
nfft = 1024

if __name__ == '__main__':
    dec_rate = int(sample_rate / fm_bandwidth)
    channel_sample_rate = sample_rate / dec_rate

    dec_rate_audio = int(channel_sample_rate / desired_audio_sample_rate)
    audio_sample_rate = channel_sample_rate / dec_rate_audio

    # design the de-emphasis filter
    d = channel_sample_rate * 50e-6  # calculate the # of samples to hit the -3dB point
    x = np.exp(-1 / d)  # calculate the decay between each sample
    b = [1 - x]
    a = [1, -x]

    # obtain the data from the SDR device
    print("Listen to the radio ...")
    sdr = RtlSdr()
    sdr.sample_rate = sample_rate
    sdr.center_freq = center_freq
    num_samples = int((t_length * sample_rate // nfft) * nfft)
    iq_samples = sdr.read_samples(num_samples)
    sdr.close()
    print("... finalized")

    # digital signal processing
    for station_id, station_freq in stations.items():
        print("Digital signal processing ...")
        start_time = time.time()

        # mixing to the base band and downsampling to channel bandwidth
        offset_freq = station_freq - center_freq
        factor = np.exp(-1.0j * 2.0 * np.pi * offset_freq / sample_rate * np.arange(len(iq_samples)))
        baseband_samples = iq_samples * factor
        samples_in_bandwidth = signal.decimate(baseband_samples, dec_rate, ftype="fir", n=15)

        # polar discriminator
        channel_samples = np.angle(samples_in_bandwidth[1:] * np.conj(samples_in_bandwidth[:-1]))
        t_range = np.arange(0, len(channel_samples)) / channel_sample_rate

        # de-emphasis filter
        deemphasised_channel_samples = signal.lfilter(b, a, channel_samples)

        # downsampling to audio frequency
        audio_samples = signal.decimate(deemphasised_channel_samples, dec_rate_audio, ftype="fir", n=15)
        end_time = time.time()
        print("... took {} s".format((end_time - start_time) * 1))

        # playing the audio
        print("Playing audio signal for station '{}' ...".format(station_id))
        sd.play(audio_samples, audio_sample_rate)

        # plotting data
        plt.figure()
        plt.scatter(np.real(samples_in_bandwidth[0:50000]), np.imag(samples_in_bandwidth[0:50000]),
                    color="red", alpha=0.05)
        plt.title("Constellation plot at {} MHz '{}'".format(station_freq / 1.0e6, station_id))
        plt.xlabel("Real")
        plt.xlim(-1.1, 1.1)
        plt.ylabel("Imag")
        plt.ylim(-1.1, 1.1)
        plt.show()

        plt.figure()
        plt.psd(channel_samples, NFFT=nfft, Fs=channel_sample_rate)
        plt.title("Channel signal at {} MHz '{}'".format(station_freq / 1.0e6, station_id))
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Power (dB)")
        plt.show()

        plt.figure()
        plt.psd(iq_samples, NFFT=nfft, Fc=center_freq / 1.0e6, Fs=sample_rate / 1.0e6)
        plt.title("Recorded spectrum")
        plt.xlabel("Frequency (MHz)")
        plt.ylabel("Power (dB)")
        plt.show()
