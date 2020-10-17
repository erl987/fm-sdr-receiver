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
#include <utility>
#include "DataProcessing.h"

/*@{*/
/** \ingroup Core
*/
namespace Core {
	namespace Processing {
		namespace Filter {
			/** \ingroup Core
			*	Support class for creating digital filters.
			*/
			template <typename U, typename T>
			class CFilterUtils
			{
			public:
				CFilterUtils(void);
				~CFilterUtils(void) {};
				template <typename InIt> CFilterUtils(InIt bFirst, InIt bLast, int downsamplingFactor, int upsamplingFactor);
				std::vector<T> Processing(std::vector<T>&& signal, std::function<std::vector<T>(std::vector<T>&&)> filter);
				std::vector<T> Processing(std::vector<T>&& signal, std::function<std::vector<T>(std::vector<T>&&)> filter, const bool& resetFirstDatapoint);
				std::pair<std::vector<U>, std::vector<T>> Processing(std::vector<U>&& time, std::vector<T>&& signal, std::function<std::vector<T>(std::vector<T>&&)> filter);
				template <typename InIt> void SetParams(InIt bFirst, InIt bLast, int downsamplingFactor, int upsamplingFactor);
				int GetDownsamplingFactor();
				int GetUpsamplingFactor();
				std::vector<T>& GetPreviousSignalRef();
				std::vector<T>& GetBRef();
				int GetFirstDatapoint();
				int DownsamplingLength(int dataLength);
				int UpsamplingLength(int dataLength, const bool& isForTimeData);
				std::vector<U> UpsamplingForTime(std::vector<U>&& data);
				std::vector<T> UpsamplingForSignal(std::vector<T>&& data);
				template <typename V> std::vector<V> Downsampling(std::vector<V>&& data);
			private:
				void SetFirstDatapoint(const size_t& dataLength);
				template <typename V> std::vector<V> UpsamplingInternal(std::vector<V>&& data, std::vector<V>&& upsampledData);
				CFilterUtils(const CFilterUtils&) = delete;
				CFilterUtils& operator= (const CFilterUtils&) = delete;

				std::vector<T> b;
				std::vector<T> previousSignal;
				std::vector<U> previousDatapointUpsamplingTime;
				std::vector<T> previousDatapointUpsamplingSignal;
				int downsamplingFactor;
				int upsamplingFactor;
				int firstDatapoint;
				bool isInit;
			};
		}
	}
}

/*@}*/


/**
* 	@brief		Default constructor.
*/
template <typename U, typename T>
Core::Processing::Filter::CFilterUtils<U, T>::CFilterUtils(void)
	: isInit(false)
{
}


/**	@brief 		Constructor
*	@param		bFirst					Iterator to the beginning of the container storing the b-filter parameters of the FIR-filter
*	@param		bLast					Iterator to one element after the end of the container storing the b-filter parameters
*	@param		downsamplingFactor		Downsampling factor: every downsamplingFactor-th datapoint is used
*	@param		upsamplingFactor		Upsampling factor: diving the section between two datapoints in upsamplingFactor sections
*	@exception	std::runtime_error		Thrown if the downsampling or the upsampling factor is smaller than 1
*	@remarks							Upsampling and downsampling factors are automatically minimized by the greatest common divisor.
*/
template <typename U, typename T>
template <typename InIt>
Core::Processing::Filter::CFilterUtils<U, T>::CFilterUtils(InIt bFirst, InIt bLast, int downsamplingFactor, int upsamplingFactor)
	: isInit(false)
{
	SetParams(bFirst, bLast, downsamplingFactor, upsamplingFactor);
}


/**	@brief 		Safely resetting filter parameters (also possible for an already used object)
*	@param		bFirst					Iterator to the beginning of the container storing the b-filter parameters of the FIR-filter
*	@param		bLast					Iterator to one element after the end of the container storing the b-filter parameters
*	@param		downsamplingFactor		Downsampling factor: every downsamplingFactor-th datapoint is used
*	@param		upsamplingFactor		Upsampling factor: diving the section between two datapoints in upsamplingFactor sections
*	@return								None
*	@exception	std::length_error		Thrown if the filter parameters are empty
*	@exception	std::runtime_error		Thrown if the downsampling or the upsampling factor is smaller than 1
*	@remarks							Upsampling and downsampling factors are automatically minimized by the greatest common divisor.
*/
template <typename U, typename T>
template <typename InIt>
void Core::Processing::Filter::CFilterUtils<U, T>::SetParams(InIt bFirst, InIt bLast, int downsamplingFactor, int upsamplingFactor)
{
	using namespace std;
	int gcd;

	//minimization of factors
	gcd = CDataProcessing<T>::GreatestCommonDivisor(downsamplingFactor, upsamplingFactor);
	downsamplingFactor = downsamplingFactor / gcd;
	upsamplingFactor = upsamplingFactor / gcd;

	// check factors
	if ((downsamplingFactor < 1) || (upsamplingFactor < 1)) {
		throw std::runtime_error("Downsampling and upsampling factors must be at least 1.");
	}

	// assign input values
	firstDatapoint = 0;
	b.assign(bFirst, bLast);
	if (b.empty()) {
		throw std::length_error("Filter parameters must not be empty.");
	}

	CFilterUtils<U, T>::downsamplingFactor = downsamplingFactor;
	CFilterUtils<U, T>::upsamplingFactor = upsamplingFactor;

	// set storage container for recently processed data (required for restarting the filtering without losses)
	previousSignal.assign(b.size() - 1, 0);
	previousDatapointUpsamplingTime.clear();
	previousDatapointUpsamplingSignal.clear();

	isInit = true;
}


/** @brief		Upsampling of time data
*	@param		data					Time data to be processed
*	@return								Processed time data
*	@exception	std::rutime_error		Thrown if the class has not been initialized before use with the function SetParams() or the constructor
*	@remarks							The size of the output container can be determined by CIIRFilterBase<T>::UpsamplingLength.
*/
template <typename U, typename T>
std::vector<U> Core::Processing::Filter::CFilterUtils<U, T>::UpsamplingForTime(std::vector<U>&& data)
{
	std::vector<U> upsampledData;

	// check if filter parameters have been initialized
	if (!isInit) {
		throw std::runtime_error("Object has not been initialized before use.");
	}

	upsampledData.resize(UpsamplingLength(static_cast<int>(data.size()), true));

	// add last previous datapoint for interpolation in front of the current data
	if ((!previousDatapointUpsamplingTime.empty()) && (!data.empty())) { // it will not be used for the first execution of the function
		upsampledData.resize(upsampledData.size() + 1);
		data.insert(data.begin(), previousDatapointUpsamplingTime.begin(), previousDatapointUpsamplingTime.end());
	}

	UpsamplingInternal(move(data), move(upsampledData));

	// delete last previous datapoint
	if ((!previousDatapointUpsamplingTime.empty()) && (!upsampledData.empty())) {
		upsampledData.erase(upsampledData.begin());
	}

	// save last datapoint for further processing steps
	if (data.size() > 0) {
		previousDatapointUpsamplingTime.assign(1, data.back());
	}

	return upsampledData;
}


/** @brief		Upsampling of signal data
*	@param		data					Signal data to be processed
*	@return								Processed signal data
*	@exception	std::rutime_error		Thrown if the class has not been initialized before use with the function SetParams() or the constructor
*	@remarks							The size of the output container can be determined by CIIRFilterBase<T>::UpsamplingLength.
*/
template <typename U, typename T>
std::vector<T> Core::Processing::Filter::CFilterUtils<U, T>::UpsamplingForSignal(std::vector<T>&& data)
{
	std::vector<T> upsampledData;

	// check if filter parameters have been initialized
	if (!isInit) {
		throw std::runtime_error("Object has not been initialized before use.");
	}

	upsampledData.resize(UpsamplingLength(static_cast<int>(data.size()), false));

	// add last previous datapoint for interpolation in front of the current data
	if ((!previousDatapointUpsamplingSignal.empty()) && (!data.empty())) { // it will not be used for the first execution of the function
		upsampledData.resize(upsampledData.size() + 1);
		data.insert(data.begin(), previousDatapointUpsamplingSignal.begin(), previousDatapointUpsamplingSignal.end());
	}

	UpsamplingInternal(move(data), move(upsampledData));

	// delete last previous datapoint
	if ((!previousDatapointUpsamplingSignal.empty()) && (!upsampledData.empty())) {
		upsampledData.erase(upsampledData.begin());
	}

	// save last datapoint for further processing steps
	if (data.size() > 0) {
		previousDatapointUpsamplingSignal.assign(1, data.back());
	}

	return upsampledData;
}


/** @brief		Actual implementation of the upsampling
*	@param		data					Data to be processed
*	@param		upsampledData			Preallocated container that will contain the processed data within the method
*	@return								Processed data
*	@exception							None
*	@remarks							The preallocated container needs to have the correct size to hold all the data.
*/
template <typename U, typename T>
template <typename V>
std::vector<V> Core::Processing::Filter::CFilterUtils<U, T>::UpsamplingInternal(std::vector<V>&& data, std::vector<V>&& upsampledData)
{
	double invUpsamplingFactor;

	// prepare upsampled data
	for (size_t i = 0; i < data.size(); i++) {
		upsampledData[upsamplingFactor * i] = data[i];
	}

	// calculate upsampled data points by linear interpolation
	if (data.size() > 1) {
		invUpsamplingFactor = 1.0 / upsamplingFactor;
		for (size_t i = 0; i < (data.size() - 1); i++) {
			auto slope = (data[i + 1] - data[i]) * invUpsamplingFactor;
			for (int j = 1; j < upsamplingFactor; j++) {
				upsampledData[upsamplingFactor * i + j] = data[i] + slope * j;
			}
		}
	}

	return upsampledData;
}


/** @brief		Performs the filter processing.
*	@param		signal					Signal data to be processed
*	@param		filter					Function that will perform the actual filtering
*	@param		resetFirstDatapoint		If this method should reset the first datapoint to be used in the next call to this method
*	@return								Processed data
*	@exception							None
*	@remarks							This is an internal method.
*/
template <typename U, typename T>
std::vector<T> Core::Processing::Filter::CFilterUtils<U, T>::Processing(std::vector<T>&& signal, std::function<std::vector<T>(std::vector<T>&&)> filter, const bool& resetFirstDatapoint)
{
	using namespace std;
	size_t signalSize;
	vector<T> filteredSignal;

	// check if filter parameters have been initialized
	if (!isInit) {
		throw runtime_error("Object has not been initialized before use.");
	}

	signalSize = signal.size();

	// upsampling
	if (upsamplingFactor > 1) {
		signal = UpsamplingForSignal(move(signal));
	}

	// downsampling
	if ((downsamplingFactor > 1) || ((downsamplingFactor == 1) && (upsamplingFactor == 1))) {
		filteredSignal = filter(move(signal));
	}
	else {
		filteredSignal = signal;
	}

	// calculate index of next datapoint for continuous downsampling
	if (resetFirstDatapoint) {
		SetFirstDatapoint(signalSize);
	}

	return filteredSignal;
}


/** @brief		Performs the filter processing.
*	@param		signal					Signal data to be processed
*	@param		filter					Function that will perform the actual filtering
*	@return								Processed data
*	@exception							None
*	@remarks							None
*/
template <typename U, typename T>
std::vector<T> Core::Processing::Filter::CFilterUtils<U, T>::Processing(std::vector<T>&& signal, std::function<std::vector<T>(std::vector<T>&&)> filter)
{
	return Processing(move(signal), filter, true);
}


/** @brief		Performs the filter processing if also time-information belongs to the data.
*	@param		time					Time data to be processed
*	@param		signal					Signal data to be processed
*	@param		filter					Function that will perform the actual filtering
*	@return								Pair of the processed timepoint and signal data containers
*	@exception							None
*	@remarks							None
*/
template <typename U, typename T>
std::pair<std::vector<U>, std::vector<T>> Core::Processing::Filter::CFilterUtils<U, T>::Processing(std::vector<U>&& time, std::vector<T>&& signal, std::function<std::vector<T>(std::vector<T>&&)> filter)
{
	using namespace std;
	size_t signalSize;
	vector<T> filteredSignal;
	vector<U> filteredTime;

	if (time.size() != signal.size()) {
		throw logic_error("Time and signal data provided to the filter are differing in size");
	}

	signalSize = signal.size();

	// resampling of the signal (including required filtering)
	filteredSignal = Processing(move(signal), filter, false);

	// resampling of the time (without filtering)
	if (upsamplingFactor > 1) {
		time = UpsamplingForTime(move(time));
	}
	filteredTime = Downsampling(move(time));

	// calculate index of next datapoint for continuous downsampling
	SetFirstDatapoint(signalSize);

	if (filteredTime.size() != filteredSignal.size()) {
		throw std::logic_error("Filtered time and signal container differ in size");
	}

	return make_pair<>(filteredTime, filteredSignal);
}


/** @brief		Obtaining the length of a downsampled dataset
*	@param		dataLength				Length of dataset to be downsampled
*	@return								Length of downsampled dataset
*	@exception							None
*	@remarks							None
*/
template <typename U, typename T>
int Core::Processing::Filter::CFilterUtils<U, T>::DownsamplingLength(int dataLength)
{
	int downsamplingLength;

	if (dataLength < firstDatapoint) {
		downsamplingLength = 0;
	}
	else {
		downsamplingLength = static_cast<int>(ceil((dataLength - firstDatapoint) / static_cast<double>(downsamplingFactor)));
	}

	return downsamplingLength;
}


/** @brief		Obtaining the length of an upsampled dataset
*	@param		dataLength				Length of dataset to be upsampled
*	@param		isForTimeData			If the returned length is valid for time data (true) or signal data (false)
*	@return								Length of upsampled dataset
*	@exception							None
*	@remarks							None
*/
template <typename U, typename T>
int Core::Processing::Filter::CFilterUtils<U, T>::UpsamplingLength(int dataLength, const bool& isForTimeData)
{
	int upsamplingLength;

	if (dataLength < 1) {
		upsamplingLength = 0;
	}
	else {
		upsamplingLength = (dataLength - 1) * upsamplingFactor + 1;
		if ((isForTimeData && !previousDatapointUpsamplingTime.empty()) || (!isForTimeData && !previousDatapointUpsamplingSignal.empty())) {
			// if this is not the first processing of the dataset, interpolation in front of the current data is required
			upsamplingLength += upsamplingFactor - 1;
		}
	}

	return upsamplingLength;
}


/** @brief		Returns the upsampling factor
*	@return								The upsampling factor
*	@exception							None
*	@remarks							None
*/
template <typename U, typename T>
int Core::Processing::Filter::CFilterUtils<U, T>::GetUpsamplingFactor()
{
	return upsamplingFactor;
}


/** @brief		Returns the downsampling factor
*	@return								The downsampling factor
*	@exception							None
*	@remarks							None
*/
template <typename U, typename T>
int Core::Processing::Filter::CFilterUtils<U, T>::GetDownsamplingFactor()
{
	return downsamplingFactor;
}


/** @brief		Performs the downsampling of data without any filtering
*   @param		data					The data to be downsampled
*	@return								The downsampled data
*	@exception							None
*	@remarks							None
*/
template <typename U, typename T>
template <typename V>
std::vector<V> Core::Processing::Filter::CFilterUtils<U, T>::Downsampling(std::vector<V>&& data)
{
	std::vector<V> downsampledData;

	// downsampling
	downsampledData.resize(DownsamplingLength(static_cast<int>(data.size())));
	for (size_t i = 0; i < downsampledData.size(); i++) {
		downsampledData[i] = data[firstDatapoint + i * downsamplingFactor];
	}

	return downsampledData;
}


/** @brief		Obtains a reference to the previous signal data stored for usage in the next processing step
*	@return								Reference to the previous signal data stored
*	@exception							None
*	@remarks							None
*/
template <typename U, typename T>
std::vector<T>& Core::Processing::Filter::CFilterUtils<U, T>::GetPreviousSignalRef()
{
	return previousSignal;
}


/** @brief		Returns a reference to the b-filter coefficients
*	@return								Reference to the b-filter coefficients
*	@exception							None
*	@remarks							None
*/
template <typename U, typename T>
std::vector<T>& Core::Processing::Filter::CFilterUtils<U, T>::GetBRef()
{
	return b;
}


/** @brief		Returns the index of the first datapoint to be used in the next processing step
*	@return								The next first datapoint to be used
*	@exception							None
*	@remarks							None
*/
template <typename U, typename T>
int Core::Processing::Filter::CFilterUtils<U, T>::GetFirstDatapoint()
{
	return firstDatapoint;
}


/** @brief		Calculates and sets the index of the first datapoint to be used in the next processing step
*	@param		dataLength				The data length
*	@exception							None
*	@remarks							None
*/
template <typename U, typename T>
void Core::Processing::Filter::CFilterUtils<U, T>::SetFirstDatapoint(const size_t& dataLength)
{
	int maxIndex;

	maxIndex = firstDatapoint + (DownsamplingLength(dataLength) - 1) * downsamplingFactor;
	firstDatapoint = maxIndex + downsamplingFactor - dataLength;
}