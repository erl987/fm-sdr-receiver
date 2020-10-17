/*	  FM-radio - software defined radio using RTL-SDR
  Copyright (C) 2019-2020 Ralf Rettig

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
# pragma once

#include <vector>
#include <mutex>
#include <stdexcept>
#include "portaudio.h"


template <class T>
class CAudioPlayer
{
public:
	CAudioPlayer(const PaDeviceIndex & deviceIndex, const T & samplingFreq, const unsigned int& samplesPerBuf=1'000);
	~CAudioPlayer();
	void AddAudioData(std::vector<T>&& newData);
private:
	int StreamCallback(const void* input, void* output, unsigned long frameCount, const PaStreamCallbackTimeInfo* timeInfo, PaStreamCallbackFlags statusFlags);
	static int CallbackMapper(const void* input, void* output, unsigned long frameCount, const PaStreamCallbackTimeInfo* timeInfo, PaStreamCallbackFlags statusFlags, void* userData);
	PaStreamParameters DefineIOParams(const PaDeviceIndex & deviceIndex, const unsigned int& numChannels) const;
	void CheckError(PaError error) const;

	struct Buffer {
		std::vector<T> data;
		std::mutex mutex;
	} buffer;
	PaStream* activeStream;
};


template <typename T>
CAudioPlayer<T>::CAudioPlayer(const PaDeviceIndex& deviceIndex, const T& samplingFreq, const unsigned int& samplesPerBuf)
{
	using namespace std;
	const unsigned int numChannels = 1;

	CheckError(Pa_Initialize());

	PaStreamParameters currOutputStream = DefineIOParams(deviceIndex, numChannels);
	CheckError(Pa_OpenStream(&activeStream, nullptr, &currOutputStream, samplingFreq, samplesPerBuf, paNoFlag, &CAudioPlayer<T>::CallbackMapper, this));
	CheckError(Pa_StartStream(activeStream));

	buffer.data.resize(samplesPerBuf * numChannels);
}


template <class T>
CAudioPlayer<T>::~CAudioPlayer()
{
	Pa_Terminate(); // never called when `Pa_Initialize` has failed
}


template <typename T>
int CAudioPlayer<T>::StreamCallback(const void* input, void* output, unsigned long frameCount, const PaStreamCallbackTimeInfo* timeInfo, PaStreamCallbackFlags statusFlags)
{
	if (statusFlags == paOutputUnderflow) {
		std::cout << "Output underflow ..." << std::endl;
	}

	std::lock_guard<std::mutex> lock(buffer.mutex);
	if (buffer.data.size() < frameCount) {
		return paContinue;
	}

	float* outputArray = (float*)(output);
	for (size_t i = 0; i < frameCount; i++) {
		outputArray[i] = buffer.data[i];
	}

	buffer.data.erase(buffer.data.begin(), buffer.data.begin() + frameCount);

	return paContinue;
}

template <typename T>
void CAudioPlayer<T>::AddAudioData(std::vector<T>&& newData)
{
	std::lock_guard<std::mutex> lock(buffer.mutex);
	buffer.data.insert(buffer.data.end(), std::make_move_iterator(newData.begin()), std::make_move_iterator(newData.end()));
}

template <class T>
PaStreamParameters CAudioPlayer<T>::DefineIOParams(const PaDeviceIndex& deviceIndex, const unsigned int& numChannels) const
{
	PaStreamParameters ioParameters;

	if ((deviceIndex < 0) || (deviceIndex >= Pa_GetDeviceCount())) {
		throw std::runtime_error("The device index is not valid");
	}

	ioParameters.device = deviceIndex;
	ioParameters.channelCount = numChannels;
	ioParameters.sampleFormat = paFloat32;
	ioParameters.suggestedLatency = Pa_GetDeviceInfo(deviceIndex)->defaultLowOutputLatency;
	ioParameters.hostApiSpecificStreamInfo = nullptr;

	return ioParameters;
}

template <typename T>
int CAudioPlayer<T>::CallbackMapper(const void* input, void* output, unsigned long frameCount, const PaStreamCallbackTimeInfo* timeInfo, PaStreamCallbackFlags statusFlags, void* userData)
{
	assert(userData != nullptr);
	auto self = reinterpret_cast<CAudioPlayer<T>*>(userData);

	return self->StreamCallback(input, output, frameCount, timeInfo, statusFlags);
}

template <typename T>
void CAudioPlayer<T>::CheckError(PaError error) const
{
	if (error != paNoError) {
		throw std::runtime_error(Pa_GetErrorText(error));
	}
}
