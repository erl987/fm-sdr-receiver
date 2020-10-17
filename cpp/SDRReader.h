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
#pragma once
#include <memory>
#include <vector>
#include <deque>
#include <string>
#include <future>
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <cassert>
#include "ComplexVector.h"
#include "ThreadExceptionStorage.h"
#include "rtl-sdr.h"

namespace AudioSP {
	namespace Audio {
		const unsigned int samplingSizeFactor = 16'384;

		template <typename T>
		class CSDRReader {
		public:
			CSDRReader(unsigned int deviceIndex, double centerFreq, double samplingFreq, unsigned int readBufferLength=262'144, unsigned int maxAvailableSignalLength=262'144); // with default buffer length of RTL-SDR
			virtual ~CSDRReader();
			std::future<ComplexVector<T>> GetSignalData();
		private:
			CSDRReader(const CSDRReader&) = delete;
			CSDRReader & operator= (const CSDRReader&) = delete;
			void StreamCallback(unsigned char* buffer, uint32_t numReadSamples);
			static void CallbackMapper(unsigned char* buffer, uint32_t numReadSamples, void* userData);
			void GetSignalDataThread();
			void CheckError(int returnValue) const;
			void WaitForNewData() const;
			void NotifyForNewData() const;
			void SetSignalDataPromise(size_t currentSampleStreamSize);

			std::unique_ptr<std::thread> samplingThread;
			std::unique_ptr<std::thread> getSignalThread;
			struct SampleStream {
				std::vector<unsigned char> data;
				std::mutex mutex;
			} sampleStream;
			struct PendingPromises {
				std::deque<std::promise<ComplexVector<T>>> data;
				std::mutex mutex;
			} pendingPromises;
			mutable struct NewDataCondition {
				std::mutex mutex;
				std::condition_variable condition;
				bool hasNew = false;
			} newDataCondition;
			unsigned int maxAvailableSignalLength;
			std::atomic<bool> doTerminate;
			rtlsdr_dev_t* dev;
		};
	}
}

template <typename T>
AudioSP::Audio::CSDRReader<T>::CSDRReader(unsigned int deviceIndex, double centerFreq, double samplingFreq, unsigned int readBufferLength, unsigned int maxAvailableSignalLength)
	: doTerminate(false),
	  maxAvailableSignalLength(maxAvailableSignalLength)
{
	using namespace std;
	newDataCondition.hasNew = false;

	assert(maxAvailableSignalLength >= readBufferLength);
	assert(readBufferLength % samplingSizeFactor == 0);

	// initialize the SDR-device
	CheckError(rtlsdr_open(&dev, deviceIndex));
	try
	{
		CheckError(rtlsdr_set_center_freq(dev, static_cast<unsigned int>(centerFreq)));
		CheckError(rtlsdr_set_sample_rate(dev, static_cast<unsigned int>(samplingFreq)));

		// activating automatic gain control
		CheckError(rtlsdr_set_tuner_gain_mode(dev, false));
		CheckError(rtlsdr_set_agc_mode(dev, true));

		// without resetting the USB-dongle will not start to capture
		CheckError(rtlsdr_reset_buffer(dev));

		samplingThread = unique_ptr<thread>(new thread([&]()
		{
			try
			{
				CheckError(rtlsdr_read_async(dev, &CSDRReader<T>::CallbackMapper, this, 0, readBufferLength));
			}
			catch (...)
			{
				ThreadExceptionStorage::getInstance().storeException(current_exception());
			}
		}));
		getSignalThread = unique_ptr<thread>(new thread(&CSDRReader<T>::GetSignalDataThread, this));
	}
	catch (const exception& e)
	{
		rtlsdr_close(dev); // ensure correct deinitialization in any situation
		throw;
	}
}

template <typename T>
AudioSP::Audio::CSDRReader<T>::~CSDRReader()
{
	doTerminate = true;
	NotifyForNewData();
	rtlsdr_cancel_async(dev); // TODO: for some reason this does NOT cancel the capture ...
	samplingThread->join();
	getSignalThread->join();
	rtlsdr_close(dev);
}

template <typename T>
void AudioSP::Audio::CSDRReader<T>::WaitForNewData() const
{
	using namespace std;
	unique_lock<mutex> lock(newDataCondition.mutex);
	newDataCondition.condition.wait(lock, [=](){ return newDataCondition.hasNew || doTerminate; });
	newDataCondition.hasNew = false;
}

template <typename T>
void AudioSP::Audio::CSDRReader<T>::NotifyForNewData() const
{
	using namespace std;
	{
		lock_guard<mutex> lock(newDataCondition.mutex);
		newDataCondition.hasNew = true;
	}
	newDataCondition.condition.notify_all();
}

template <typename T>
void AudioSP::Audio::CSDRReader<T>::GetSignalDataThread()
{
	using namespace std;
	size_t currentSampleStreamSize;
	try
	{
		while (!doTerminate)
		{
			WaitForNewData();

			while (true)
			{
				{
					lock_guard<mutex> lock(pendingPromises.mutex);
					if (pendingPromises.data.empty())
					{
						break;
					}
				}
				{
					lock_guard<mutex> sampleStreamLock(sampleStream.mutex);
					currentSampleStreamSize = sampleStream.data.size();
					if (currentSampleStreamSize > maxAvailableSignalLength)
					{
						currentSampleStreamSize = maxAvailableSignalLength;
					}
				}
				if (currentSampleStreamSize == 0)
				{
					break;
				}

				SetSignalDataPromise(currentSampleStreamSize);
			}
		}
	}
	catch (...)
	{
		ThreadExceptionStorage::getInstance().storeException(current_exception());
	}
}

template <typename T>
void AudioSP::Audio::CSDRReader<T>::SetSignalDataPromise(size_t currentSampleStreamSize)
{
	using namespace std;
	vector<T> real(currentSampleStreamSize / 2);
	vector<T> imag(currentSampleStreamSize / 2);

	{
		// convert the IQ-datastream to an analytical signal
		lock_guard<mutex> sampleStreamLock(sampleStream.mutex);
		T reciprocalVal = static_cast<T>(2.0 / 255);
		for (size_t i = 0; i < currentSampleStreamSize / 2; i++)
		{
			real[i] = static_cast<T>(sampleStream.data[2 * i])* reciprocalVal - 1.0;
			imag[i] = static_cast<T>(sampleStream.data[2 * i + 1])* reciprocalVal - 1.0;
		}

		if (sampleStream.data.size() <= currentSampleStreamSize) {
			sampleStream.data.clear();
		}
		else
		{
			auto beginIt = sampleStream.data.begin();
			sampleStream.data.erase(beginIt, beginIt + currentSampleStreamSize);
		}
	}

	{
		lock_guard<mutex> lock(pendingPromises.mutex);
		ComplexVector<T> newSignal = ComplexVector<T>(move(real), move(imag));
		pendingPromises.data.front().set_value(move(newSignal));
		pendingPromises.data.pop_front();
	}
}

template <typename T>
std::future<ComplexVector<T>> AudioSP::Audio::CSDRReader<T>::GetSignalData()
{
	using namespace std;
	promise<ComplexVector<T>> promise;
	auto future = promise.get_future();

	{
		lock_guard<mutex> lock(pendingPromises.mutex);
		pendingPromises.data.push_back(move(promise));
	}

	return future;
}

template <typename T>
void AudioSP::Audio::CSDRReader<T>::CallbackMapper(unsigned char* buffer, uint32_t numReadSamples, void* userData)
{
	assert(userData != nullptr);
	auto self = reinterpret_cast<CSDRReader<T>*>(userData);

	self->StreamCallback(buffer, numReadSamples);
}

template <typename T>
void AudioSP::Audio::CSDRReader<T>::StreamCallback(unsigned char* buffer, uint32_t numReadSamples)
{
	using namespace std;

	if ((numReadSamples % 2) != 0)
	{
		throw logic_error("The number of samples received from the SDR-receiver is not consistent (i.e. uneven).");
	}

	{
		lock_guard<mutex> samplesStreamLock(sampleStream.mutex);
		size_t prevLength = sampleStream.data.size();

		sampleStream.data.resize(sampleStream.data.size() + numReadSamples);
		for (uint32_t i = 0; i < numReadSamples; i++)
		{
			sampleStream.data[i + prevLength] = buffer[i];
		}
	}

	NotifyForNewData();
}

template <typename T>
void AudioSP::Audio::CSDRReader<T>::CheckError(int returnValue) const
{
	if (returnValue != 0)
	{
		throw std::runtime_error("An error has occurred in the RTL-SDR reader, error code: " + std::to_string(returnValue));
	}
}