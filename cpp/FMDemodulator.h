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
#include <vector>
#include <complex>
#include <chrono>
#include <thread>
#include <future>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <boost/math/constants/constants.hpp>
#include "ComplexFIRfilter.h"
#include "FIRfilter.h"
#include "IIRfilter.h"
#include "ComplexVector.h"
#include "Limiter.h"
#include "ThreadExceptionStorage.h"

struct FMDemodulatorSettings {
	double relDistanceToCenterFreq = 0.2;
	double desiredSampleTimePeriod = 0.05; // in s
	unsigned int factorForSignalSampleBackup = 2;
	double desiredFmBandwidth = 200.0e3; // in Hz
	double desiredAudioSampleRate = 48.0e3; // in Hz
	unsigned int filterOrderDownsamplingToChannel = 20;
	unsigned int filterOrderDownsamplingToAudio = 20;
	unsigned int filterOrderDc = 3;
	unsigned int cutoffFreqDc = 300; // in Hz
	double deemphasisDecayTime = 50.0e-6; // in s - valid for Europe
	double limiterAttackCoeff = 0.9;
	double limiterReleaseCoeff = 0.9999;
	double limiterThreshold = 0.9;
	unsigned int limiterNumSamplesDelay = 40;
};


template <typename T>
class FMDemodulator
{
public:
	FMDemodulator(double stationFreq, double samplingFreq, const FMDemodulatorSettings& settings);
	~FMDemodulator();
	std::future<std::vector<T>> Perform(ComplexVector<T>&& iqSamples);
	double GetAudioSamplingFreq() const;
	double GetCenterFreq() const;
	unsigned int GetNumSignalSamplesForSDR() const;
	unsigned int GetMaxNumSignalSamplesForProcessing() const;
	std::chrono::microseconds GetSampleTimePeriod() const;
private:
	ComplexVector<T> PrecalcSine(size_t length, double offsetFreq, double samplingFreq) const;
	std::vector<T> ProcessDataset(ComplexVector<T>&& iqSamples);
	void ProcessingThread();
	void WaitForNewData() const;
	void NotifyForNewData() const;

	std::unique_ptr<std::thread> processingThread;
	std::unique_ptr<Core::Processing::Filter::ComplexFIRfilter<T>> downsamplerToChannel;
	std::unique_ptr<Core::Processing::Filter::CFIRfilter<T, T>> downsamplerToAudio;
	std::unique_ptr<Core::Processing::Filter::CIIRfilter<T, T>> deemphasisFilter;
	std::unique_ptr<Core::Processing::Filter::CIIRfilter<T, T>> dcFilter;
	std::unique_ptr<Limiter<T>> limiter;
	FMDemodulatorSettings settings;
	double audioSamplingFreq;
	double centerFreq;
	unsigned int numSdrSignalSamples;
	unsigned int maxNumSamplesForProcessing;
	std::chrono::microseconds sampleTimePeriod;
	T phi;
	T deltaPhase;
	ComplexVector<T> sine;
	std::complex<T> lastBandwidthSample;
	struct PendingData {
		std::deque<std::promise<std::vector<T>>> promises;
		std::deque<ComplexVector<T>> iqSamples;
		std::mutex mutex;
	} pendingData;
	mutable struct NewDataCondition {
		std::condition_variable condition;
		bool hasNew = false;
		std::mutex mutex;
	} newDataCondition;
	std::atomic<bool> doTerminate;
};

template <typename T>
FMDemodulator<T>::FMDemodulator(double stationFreq, double samplingFreq, const FMDemodulatorSettings& settings)
	: settings(settings),
	  doTerminate(false)
{
	using namespace std;
	using namespace boost::math::constants;
	using namespace Core::Processing::Filter;

	double channelSamplingFreq;
	unsigned int downsamplingFactorToChannel, downsamplingFactorToAudio;

	downsamplingFactorToChannel = static_cast<unsigned int>(round(samplingFreq / settings.desiredFmBandwidth));
	channelSamplingFreq = samplingFreq / downsamplingFactorToChannel;
	downsamplingFactorToAudio = static_cast<unsigned int>(round(channelSamplingFreq / settings.desiredAudioSampleRate));
	audioSamplingFreq = channelSamplingFreq / downsamplingFactorToAudio;
	downsamplerToChannel = make_unique<ComplexFIRfilter<T>>(downsamplingFactorToChannel, 1, settings.filterOrderDownsamplingToChannel);
	downsamplerToAudio = make_unique<CFIRfilter<T, T>>(downsamplingFactorToAudio, 1, settings.filterOrderDownsamplingToAudio, 1e-5);

	vector<T> aDeemphasis = { 1, -exp(-1 / static_cast<T>(channelSamplingFreq * settings.deemphasisDecayTime)) };
	vector<T> bDeemphasis = { 1 - exp(-1 / static_cast<T>(channelSamplingFreq * settings.deemphasisDecayTime)) };
	deemphasisFilter = make_unique<CIIRfilter<T, T>>(aDeemphasis.cbegin(), aDeemphasis.cend(), bDeemphasis.cbegin(), bDeemphasis.cend());

	auto dcFilterParams = CIIRfilter<T, T>::DesignHighPassFilter(settings.filterOrderDc, settings.cutoffFreqDc, audioSamplingFreq);
	vector<T> aDC = get<1>(dcFilterParams);
	vector<T> bDC = get<0>(dcFilterParams);
	dcFilter = make_unique<CIIRfilter<T, T>>(aDC.cbegin(), aDC.cend(), bDC.cbegin(), bDC.cend());

	limiter = make_unique<Limiter<T>>(settings.limiterAttackCoeff, settings.limiterReleaseCoeff, settings.limiterNumSamplesDelay, settings.limiterThreshold);
	
	double freqOffset = -settings.relDistanceToCenterFreq * samplingFreq;
	centerFreq = stationFreq - freqOffset;

	numSdrSignalSamples = AudioSP::Audio::samplingSizeFactor * ceil((2 * samplingFreq * settings.desiredSampleTimePeriod) / AudioSP::Audio::samplingSizeFactor);
	sampleTimePeriod = chrono::microseconds(static_cast<long>(1.0e6 * numSdrSignalSamples / samplingFreq / 2)); // per sampling two samples (real + imaginary) are acquired

	// the data is processed as complex numbers, i.e. having half the length of the signal stream
	maxNumSamplesForProcessing = settings.factorForSignalSampleBackup * numSdrSignalSamples / 2;
	sine = PrecalcSine(maxNumSamplesForProcessing, freqOffset, samplingFreq);
	deltaPhase = static_cast<T>(2.0 * pi<double>() * freqOffset / samplingFreq);
	phi = -deltaPhase;

	processingThread = unique_ptr<thread>(new thread(&FMDemodulator<T>::ProcessingThread, this));
}

template <typename T>
FMDemodulator<T>::~FMDemodulator()
{
	doTerminate = true;
	NotifyForNewData();
	processingThread->join();
}

template <typename T>
std::vector<T> FMDemodulator<T>::ProcessDataset(ComplexVector<T>&& iqSamples)
{
	using namespace std;
	using namespace boost::math::constants;

	assert(iqSamples.size() <= maxNumSamplesForProcessing);
	if (iqSamples.size() == 0)
	{
		return vector<T>();
	}

	// mixing to the base band
	complex<T> multiplicator(cos( -phi - deltaPhase ), sin( -phi - deltaPhase ));
	ComplexVector<T> thisSine(sine, iqSamples.size());
	ComplexVector<T> factor = multiplicator * thisSine;
	iqSamples *= factor;
	phi = 2 * pi<T>() - arg(factor.at(factor.size() - 1));
	ComplexVector<T> basebandSamples = move(iqSamples);
	ComplexVector<T> samplesInBandwidth = downsamplerToChannel->Processing(move(basebandSamples));

	if (samplesInBandwidth.size() == 0)
	{
		return vector<T>();
	}

	// polar discriminator
	complex<T> nextLastBandwidthSample = complex<T>(samplesInBandwidth.at(samplesInBandwidth.size() - 1));
	vector<T> channelSamples = polarDiscriminator(move(samplesInBandwidth), move(lastBandwidthSample));
	lastBandwidthSample = nextLastBandwidthSample;

	// de-emphasis filter
	vector<T> deemphasisedChannelSamples = deemphasisFilter->Processing(move(channelSamples));

	// downsampling to audio frequency
	vector<T> audioSamplesWithDc = downsamplerToAudio->Processing(move(deemphasisedChannelSamples));

	// DC-offset filter
	vector<T> unlimitedAudioSamples = dcFilter->Processing(move(audioSamplesWithDc));

	// limiter
	return limiter->Processing(move(unlimitedAudioSamples));
}

template <typename T>
void FMDemodulator<T>::WaitForNewData() const
{
	using namespace std;
	unique_lock<mutex> lock(newDataCondition.mutex);
	newDataCondition.condition.wait(lock, [=](){return newDataCondition.hasNew || doTerminate;});
	newDataCondition.hasNew = false;
}

template <typename T>
void FMDemodulator<T>::NotifyForNewData() const
{
	using namespace std;
	{
		lock_guard<mutex> lock(newDataCondition.mutex);
		newDataCondition.hasNew = true;
	}
	newDataCondition.condition.notify_all();
}


template <typename T>
void FMDemodulator<T>::ProcessingThread()
{
	using namespace std;

	try
	{
		while (!doTerminate)
		{
			WaitForNewData();

			while (true)
			{
				unique_lock<mutex> lock(pendingData.mutex);
				if (pendingData.promises.empty())
				{
					break;
				}

				ComplexVector<T> iqSamples = move(pendingData.iqSamples.front());
				promise<vector<T>> promise = move(pendingData.promises.front());
				pendingData.iqSamples.pop_front();
				pendingData.promises.pop_front();
				lock.unlock();

				promise.set_value(ProcessDataset(move(iqSamples)));	
			}
		}
	} 
	catch (...)
	{
		ThreadExceptionStorage::getInstance().storeException(current_exception());
	}
}

template <typename T>
std::future<std::vector<T>> FMDemodulator<T>::Perform(ComplexVector<T>&& iqSamples)
{
	using namespace std;
	promise<vector<T>> promise;
	auto future = promise.get_future();

	{
		lock_guard<mutex> lock(pendingData.mutex);
		pendingData.promises.push_back(move(promise));
		pendingData.iqSamples.push_back(move(iqSamples));
	}
	NotifyForNewData();

	return future;
}

template <typename T>
ComplexVector<T> FMDemodulator<T>::PrecalcSine(size_t length, double offsetFreq, double samplingFreq) const
{
	using namespace std;
	using namespace boost::math::constants;

	vector<T> real(length);
	vector<T> imag(length);

	double constVal = -2.0 * pi<double>() * offsetFreq / samplingFreq;
	for (size_t n = 0; n < length; n++)
	{
		double phi = constVal * n;
		real[n] = static_cast<T>(cos(phi));
		imag[n] = static_cast<T>(sin(phi));
	}

	return ComplexVector<T>(move(real), move(imag));
}

template <typename T>
double FMDemodulator<T>::GetAudioSamplingFreq() const
{
	return audioSamplingFreq;
}

template <typename T>
double FMDemodulator<T>::GetCenterFreq() const
{
	return centerFreq;
}

template <typename T>
unsigned int FMDemodulator<T>::GetNumSignalSamplesForSDR() const
{
	return numSdrSignalSamples;
}

template <typename T>
unsigned int FMDemodulator<T>::GetMaxNumSignalSamplesForProcessing() const
{
	return 2 * maxNumSamplesForProcessing; // valid for the SDR-data stream that contains real + imaginary part
}

template <typename T>
std::chrono::microseconds FMDemodulator<T>::GetSampleTimePeriod() const
{
	return sampleTimePeriod;
}