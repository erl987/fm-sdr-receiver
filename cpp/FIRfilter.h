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
#include <limits>
#include <numeric>
#include <algorithm>
#include <math.h>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/pow.hpp>
#include "interpolation.h" // Alglib
#include "FilterUtils.h"
#include "DataProcessing.h"
#include "FFT.h"
#include "SearchTransferFunc.h"

/*@{*/
/** \ingroup Core
*/
namespace Core {
	namespace Processing {
		namespace Filter {
			/** \ingroup Core
			*	Class for filtering data with digital FIR-filters in the time domain. The class is especially well-suited for processing arbitrary length input signals.
			*/
			template <class U, class T>
			class CFIRfilter
			{
			public:
				CFIRfilter(void);
				~CFIRfilter(void) {};
				template <class InIt> CFIRfilter(InIt bFirst, InIt bLast, int downsamplingFactor = 1, int upsamplingFactor = 1, T maxSymmetryDevAllowed = 1e-12);
				CFIRfilter(int downsamplingFactor, int upsamplingFactor, int filterOrder = 30, T maxSymmetryDevAllowed = 1e-12);
				template <class InIt> void SetParams(InIt bFirst, InIt bLast, int downsamplingFactor = 1, int upsamplingFactor = 1, T maxSymmetryDevAllowed = 1e-12);
				void SetParams(int downsamplingFactor, int upsamplingFactor, int filterOrder = 30, T maxSymmetryDevAllowed = 1e-12);
				static std::vector<T> DesignLowPassFilter(int filterOrder, T normCutoffFreq);
				static std::vector<T> DesignLowPassFilter(const T& transitionWidth, const T& cutoffFreq, const T& samplingFreq, const int& startFilterOrder = 500, const int& maxFilterOrder = 1000);
				std::vector<T> Processing(std::vector<T>&& signal);
				std::pair<std::vector<U>, std::vector<T>> Processing(std::vector<U>&& time, std::vector<T>&& signal);
				static T CalcTransitionWidth(const int& filterOrder, const T& normCutoffFreq, const T& samplingFreq, const T& hLimit = 0.975, const T& deltaF = 20, const T& maxResidual = 1e-5, const int& maxNumIt = 500);
			private:
				template <class InIt, class OutIt> static OutIt DigitalFilterGain(const long& pointNumber, const T& samplingFreq, InIt bFirst, InIt bLast, OutIt gainFirst);
				template <class InIt> static void GenerateFreqResponseModel(InIt bFirst, InIt bLast, const T& deltaF, const T& samplingFreq, alglib::spline1dinterpolant& model);
				static T GetFreqResponse(const T& freq, const alglib::spline1dinterpolant& model);

				std::vector<T> DownsamplingFilter(std::vector<T>&& inputSignal);
				CFIRfilter(const CFIRfilter&) = delete;
				CFIRfilter& operator= (const CFIRfilter&) = delete;

				CFilterUtils<U, T> utils;
				bool isInit;
			};
		}
	}
}

/*@}*/


/**	@brief 		Filtering / resampling function depending on the set upsampling and downsampling factors
*	@param		signal					Container storing the data to be resampled
*	@return								Container storing the downsampled data
*	@exception	std::rutime_error		Thrown if the class has not been initialized before use with the function SetParams() or the constructor
*	@remarks							The algorithm follows the IEEE Programs for Digital Signal Processing. Wiley & Sons 1979, program 8.2.
*/
template <class U, class T>
std::vector<T> Core::Processing::Filter::CFIRfilter<U, T>::Processing(std::vector<T>&& signal)
{
	using namespace std::placeholders;
	return utils.Processing(move(signal), std::bind(&CFIRfilter<U, T>::DownsamplingFilter, this, _1));
}


/**	@brief		Default constructor
*/
template <class U, class T>
Core::Processing::Filter::CFIRfilter<U, T>::CFIRfilter(void)
{
}


/**	@brief 		Resampling function including the corresponding time data container (without filtering the time data)
*	@param		time					Container storing the time corresponding to the data to be resampled
*	@param		signal					Container storing the data to be resampled. The length has to be the same as for the time container.
*	@return								Pair of containers storing the downsampled time and data
*	@exception	std::rutime_error		Thrown if the class has not been initialized before use with the function SetParams() or the constructor
*	@remarks							The algorithm follows the IEEE Programs for Digital Signal Processing. Wiley & Sons 1979, program 8.2.
*/
template <class U, class T>
std::pair<std::vector<U>, std::vector<T>> Core::Processing::Filter::CFIRfilter<U, T>::Processing(std::vector<U>&& time, std::vector<T>&& signal)
{
	return utils.Processing(move(time), move(signal), std::bind(&CFIRfilter<U, T>::DownsamplingFilter, this, std::placeholders::_1));
}


/**	@brief 		Constructor
*	@param		bFirst					Iterator to the beginning of the container storing the filter parameters of the FIR-filter
*	@param		bLast					Iterator to the end of the container storing the filter parameters
*	@param		downsamplingFactor		Downsampling factor: every downsamplingFactor-th datapoint is used. Can be leaved out together with upsamplingFactor if only filtering is used.
*	@param		upsamplingFactor		Upsampling factor: diving the section between two datapoints in upsamplingFactor sections. Can be leaved out if only filtering is used.
*	@param		maxSymmetryDevAllowed	Maximum of the (symmetrical) b-parameters allowed to assume symmetry. This is to handle very small deviations which might occur due to storage of the parameters on file. Can be left out, default value is 1e-12.
*	@exception	std::length_error		Thrown if the filter order is uneven (i.e. the number of coefficients is even)
*	@exception	std::logic_error		Thrown if the filter is not symmetric (criterium is maxSymmetryDevAllowed)
*	@exception	std::runtime_error		Thrown if the downsampling or upsampling factor is smaller than 1.
*	@remarks							Upsampling and downsampling factors are automatically minimized by the greatest common divisor. In case of downsampling, the filter must be correctly designed regarding the cut-off frequency at fc = ( upsamplingFactor / downsamplingFactor * sampling rate / 2 ).
*										The algorithm follows the IEEE Programs for Digital Signal Processing. Wiley & Sons 1979, program 8.2. A 30-order FIR lowpass filter is recommended for downsampling.
*/
template <class U, class T>
template <class InIt> 
Core::Processing::Filter::CFIRfilter<U, T>::CFIRfilter(InIt bFirst, InIt bLast, int downsamplingFactor, int upsamplingFactor, T maxSymmetryDevAllowed)
{
	SetParams(bFirst, bLast, downsamplingFactor, upsamplingFactor, maxSymmetryDevAllowed);
}


/**	@brief 		Constructor for resampling with automatic design of the required low-pass filter
*	@param		downsamplingFactor			Downsampling factor: every downsamplingFactor-th datapoint is used. Can be leaved out together with upsamplingFactor if only filtering is used.
*	@param		upsamplingFactor			Upsampling factor: diving the section between two datapoints in upsamplingFactor sections. Can be leaved out if only filtering is used.
*	@param		filterOrder					Filter order, if omitted it is 30
*	@param		maxSymmetryDevAllowed		Maximum of the (symmetrical) b-parameters allowed to assume symmetry. This is to handle very small deviations which might occur due to storage of the parameters on file. Can be left out, default value is 1e-12.
*	@exception	std::runtime_error			Thrown if the downsampling or upsampling factor is smaller than 1.
*	@remarks								Upsampling and downsampling factors are automatically minimized by the greatest common divisor. In case of downsampling, the filter is automatically designed as 30-th order FIR-filter with the ideal cut-off frequency at fc = ( upsamplingFactor / downsamplingFactor * sampling rate / 2 ).
*											The algorithm follows the IEEE Programs for Digital Signal Processing. Wiley & Sons 1979, program 8.2.
*/
template <class U, class T>
Core::Processing::Filter::CFIRfilter<U, T>::CFIRfilter(int downsamplingFactor, int upsamplingFactor, int filterOrder, T maxSymmetryDevAllowed)
{
	SetParams(downsamplingFactor, upsamplingFactor, filterOrder, maxSymmetryDevAllowed);
}


/**	@brief 		Safely resetting with automatic design of the required low-pass filter (also possible for an already used object)
*	@param		downsamplingFactor			Downsampling factor: every downsamplingFactor-th datapoint is used. Can be leaved out together with upsamplingFactor if only filtering is used.
*	@param		upsamplingFactor			Upsampling factor: diving the section between two datapoints in upsamplingFactor sections. Can be leaved out if only filtering is used.
*	@param		filterOrder					Filter order, if omitted it is 30
*	@param		maxSymmetryDevAllowed		Maximum of the (symmetrical) b-parameters allowed to assume symmetry. This is to handle very small deviations which might occur due to storage of the parameters on file. Can be left out, default value is 1e-12.
*	@return									None
*	@exception	std::runtime_error			Thrown if the downsampling or upsampling factor is smaller than 1.
*	@remarks								Upsampling and downsampling factors are automatically minimized by the greatest common divisor. In case of downsampling, the filter is by default automatically designed as 30-th order FIR-filter with the ideal cut-off frequency at fc = ( upsamplingFactor / downsamplingFactor * sampling rate / 2 ).
*											The algorithm follows the IEEE Programs for Digital Signal Processing. Wiley & Sons 1979, program 8.2.
*/
template <class U, class T>
void Core::Processing::Filter::CFIRfilter<U, T>::SetParams(int downsamplingFactor, int upsamplingFactor, int filterOrder, T maxSymmetryDevAllowed)
{
	T normCutoffFreq;
	std::vector<T> b;

	normCutoffFreq = static_cast<T>(upsamplingFactor) / downsamplingFactor;

	// standard FIR low-pass filtering for downsampling according to IEEE program 8.2	
	if (normCutoffFreq < 1) {
		b = DesignLowPassFilter(filterOrder, normCutoffFreq + static_cast<T>(3.3) / (filterOrder + 1));
	} else {
		b = DesignLowPassFilter(filterOrder, static_cast<T>(1) + static_cast<T>(3.3) / (filterOrder + 1));
	}
	SetParams(b.begin(), b.end(), downsamplingFactor, upsamplingFactor, maxSymmetryDevAllowed);
}


/**	@brief 		Safely resetting filter parameters (also possible for an already used object)
*	@param		bFirst						Iterator to the beginning of the container storing the filter parameters of the FIR-filter
*	@param		bLast						Iterator to the end of the container storing the filter parameters
*	@param		downsamplingFactor			Downsampling factor: every downsamplingFactor-th datapoint is used. Can be leaved out together with upsamplingFactor if only filtering is used
*	@param		upsamplingFactor			Upsampling factor: diving the section between two datapoints in upsamplingFactor sections. Can be leaved out if only filtering is used.
*	@param		maxSymmetryDevAllowed		Maximum of the (symmetrical) b-parameters allowed to assume symmetry. This is to handle very small deviations which might occur due to storage of the parameters on file. Can be left out, default value is 1e-12.
*	@return									None
*	@exception	std::length_error			Thrown if the filter order is uneven (i.e. the number of coefficients is even) or the filter parameters are empty
*	@exception	std::logic_error			Thrown if the filter is not symmetric (criterium is maxSymmetryDevAllowed)
*	@exception	std::runtime_error			Thrown if the downsampling or upsampling factor is smaller than 1.
*	@remarks								Upsampling and downsampling factors are automatically minimized by the greatest common divisor. In case of downsampling, the filter must be correctly designed regarding the cut-off frequency at fc = ( upsamplingFactor / downsamplingFactor * sampling rate / 2 ).
*											The algorithm follows the IEEE Programs for Digital Signal Processing. Wiley & Sons 1979, program 8.2. A 30-order FIR lowpass filter is recommended for downsampling.
*/
template <class U, class T>
template <class InIt> 
void Core::Processing::Filter::CFIRfilter<U, T>::SetParams(InIt bFirst, InIt bLast, int downsamplingFactor, int upsamplingFactor, T maxSymmetryDevAllowed)
{
	using namespace std;
	vector<T> checkSymmetry;

	utils.SetParams(bFirst, bLast, downsamplingFactor, upsamplingFactor);

	// check filter conditions
	if ((utils.GetBRef().size() % 2) != 1) {
		throw std::length_error("Filter must have even order.");
	}

	transform(utils.GetBRef().begin(), utils.GetBRef().end(), utils.GetBRef().rbegin(), back_inserter(checkSymmetry), [](T a, T b) { return (abs(a - b)); });
	if (*max_element(checkSymmetry.begin(), checkSymmetry.end()) > maxSymmetryDevAllowed) {
		throw std::logic_error("Filter must be symmetric.");
	}

	this->isInit = true;
}


/** @brief		Downsampling and filtering of a dataset. This function is specifically implemented for each filter type.
*	@param		signal						Data to be processed
*	@return									Downsampled output data
*	@exception	std::logic_error			Thrown if the minimum filter order is smaller than 8.
*	@remarks								This function performs efficiently filtering and downsampling at the same time.
*/
template <class U, class T>
std::vector<T> Core::Processing::Filter::CFIRfilter<U, T>::DownsamplingFilter(std::vector<T>&& signal)
{
	std::vector<T> filteredSignal, b;
	int filterLengthB, leftHalfIndex, rightHalfIndex, signalIndex, filteredSignalZeroPads, signalZeroPads, requiredSignalLength;
	const int innerUnrollLength = 4;
	const int outerUnrollLength = 4;

	// obtain input data
	b.assign(utils.GetBRef().begin(), utils.GetBRef().end());
	reverse(b.begin(), b.end()); // improving efficiency in the loop
	size_t downsamplingLength = static_cast<size_t>(utils.DownsamplingLength(static_cast<int>(signal.size())));
	int downsamplingFactor = utils.GetDownsamplingFactor();

	// check if the filter is large enough for inner loop unrolling in any situation
	if (b.size() < (2 * innerUnrollLength + 1)) {
		throw std::logic_error("Filter orders below 8 are not supported.");
	}

	filterLengthB = static_cast<int>(b.size());
	signalIndex = filterLengthB - 1 + utils.GetFirstDatapoint();

	// add previous data in front of the current data
	signal.insert(signal.begin(), utils.GetPreviousSignalRef().begin(), utils.GetPreviousSignalRef().end());

	// preconditioning of the data for loop unrolling
	filteredSignalZeroPads = outerUnrollLength - downsamplingLength % outerUnrollLength;
	if (filteredSignalZeroPads == outerUnrollLength) {
		filteredSignalZeroPads = 0;
	}
	filteredSignal.resize(downsamplingLength + filteredSignalZeroPads);
	requiredSignalLength = signalIndex + downsamplingFactor * filteredSignal.size(); // index of first processed datapoint + first iterator distance beyond the processed signal data points
	signalZeroPads = requiredSignalLength - signal.size();
	signal.resize(signal.size() + signalZeroPads);

	// delete the second half of the filter (not required because it is symmetric) and zero-pad it
	T bCenter = b[filterLengthB / 2];
	b.erase(b.begin() + static_cast<size_t>(filterLengthB / 2), b.end());
	b.insert(b.end(), innerUnrollLength - 1, static_cast<T>(0)); // this is the upper limit of the required length, the exact length is not important

	// perform filtering with outer and inner loop explicitly unrolled for better performance (SIMD-vectorization, literature: A. Shahbahrami, B.H.H. Juurlink, and S. Vassiliadis. Efficient vectorization of the fir filter. In ProRisc 2005, pages 432--437, November 2005.)
	for (size_t counter = 0; counter < filteredSignal.size(); counter += outerUnrollLength) {
		// symmetrical filter additions
		for (int k = -(filterLengthB - 1); k < -(filterLengthB / 2); k += innerUnrollLength) {
			leftHalfIndex = signalIndex + k;
			rightHalfIndex = signalIndex - (filterLengthB - 1) - k;
			filteredSignal[counter] += b[(filterLengthB - 1) + k] * (signal[leftHalfIndex] + signal[rightHalfIndex])
				+ b[(filterLengthB - 1) + k + 1] * (signal[leftHalfIndex + 1] + signal[rightHalfIndex - 1])
				+ b[(filterLengthB - 1) + k + 2] * (signal[leftHalfIndex + 2] + signal[rightHalfIndex - 2])
				+ b[(filterLengthB - 1) + k + 3] * (signal[leftHalfIndex + 3] + signal[rightHalfIndex - 3]);

			filteredSignal[counter + 1] += b[(filterLengthB - 1) + k] * (signal[leftHalfIndex + downsamplingFactor] + signal[rightHalfIndex + downsamplingFactor])
				+ b[(filterLengthB - 1) + k + 1] * (signal[leftHalfIndex + 1 + downsamplingFactor] + signal[rightHalfIndex - 1 + downsamplingFactor])
				+ b[(filterLengthB - 1) + k + 2] * (signal[leftHalfIndex + 2 + downsamplingFactor] + signal[rightHalfIndex - 2 + downsamplingFactor])
				+ b[(filterLengthB - 1) + k + 3] * (signal[leftHalfIndex + 3 + downsamplingFactor] + signal[rightHalfIndex - 3 + downsamplingFactor]);

			filteredSignal[counter + 2] += b[(filterLengthB - 1) + k] * (signal[leftHalfIndex + 2 * downsamplingFactor] + signal[rightHalfIndex + 2 * downsamplingFactor])
				+ b[(filterLengthB - 1) + k + 1] * (signal[leftHalfIndex + 1 + 2 * downsamplingFactor] + signal[rightHalfIndex - 1 + 2 * downsamplingFactor])
				+ b[(filterLengthB - 1) + k + 2] * (signal[leftHalfIndex + 2 + 2 * downsamplingFactor] + signal[rightHalfIndex - 2 + 2 * downsamplingFactor])
				+ b[(filterLengthB - 1) + k + 3] * (signal[leftHalfIndex + 3 + 2 * downsamplingFactor] + signal[rightHalfIndex - 3 + 2 * downsamplingFactor]);

			filteredSignal[counter + 3] += b[(filterLengthB - 1) + k] * (signal[leftHalfIndex + 3 * downsamplingFactor] + signal[rightHalfIndex + 3 * downsamplingFactor])
				+ b[(filterLengthB - 1) + k + 1] * (signal[leftHalfIndex + 1 + 3 * downsamplingFactor] + signal[rightHalfIndex - 1 + 3 * downsamplingFactor])
				+ b[(filterLengthB - 1) + k + 2] * (signal[leftHalfIndex + 2 + 3 * downsamplingFactor] + signal[rightHalfIndex - 2 + 3 * downsamplingFactor])
				+ b[(filterLengthB - 1) + k + 3] * (signal[leftHalfIndex + 3 + 3 * downsamplingFactor] + signal[rightHalfIndex - 3 + 3 * downsamplingFactor]);
		}
		// central data point addition
		for (size_t i = 0; i < outerUnrollLength; i++) {
			filteredSignal[counter + i] += bCenter * signal[signalIndex + -filterLengthB / 2 + i * downsamplingFactor];
		}

		signalIndex += outerUnrollLength * downsamplingFactor;
	}

	// remove preconditioning of the data for loop unrolling
	filteredSignal.erase(filteredSignal.end() - filteredSignalZeroPads, filteredSignal.end());
	signal.erase(signal.end() - signalZeroPads, signal.end());

	// delete temporary previous data
	signal.erase(signal.begin(), signal.begin() + utils.GetPreviousSignalRef().size());

	// preserve the latest data corresponding to the length of the filter (for restarting the filter)
	if (signal.size() >= static_cast<size_t>(filterLengthB)) {
		utils.GetPreviousSignalRef().insert(utils.GetPreviousSignalRef().end(), signal.begin() + signal.size() - (filterLengthB - 1), signal.end());
	} else {
		utils.GetPreviousSignalRef().insert(utils.GetPreviousSignalRef().end(), signal.begin(), signal.end());
	}
	utils.GetPreviousSignalRef().erase(utils.GetPreviousSignalRef().begin(), utils.GetPreviousSignalRef().begin() + utils.GetPreviousSignalRef().size() - (filterLengthB - 1));

	return filteredSignal;
}


/**	@brief		FIR-filter Hamming window based design for linearized phase low-pass filters
*	@param		filterOrder					Filter order
*	@param		normCutoffFreq				Normalized cutoff frequency of the (real) low-pass filter. It is the frequency where the main leaf of the filter has the minimum attenuation of the stop band (i.e. that of first side leaf). The frequency is normalized with the Nyquist frequency, i.e. 0.5 * sampling frequency.
*	@return									The b-filter coeffficient container.
*	@exception	std::length_error			Thrown if the filter order is not even. Only even filters are implemented.
*	@exception	std::logic_error			Thrown if the normalized ideal cutoff frequency is not in the range of (0..1], relation ideal - real low pass filter: fcNormIdeal = fcNormReal - 3.3 / ( filterOrder + 1 )
*	@remarks								The algorithm follows the IEEE Programs for Digital Signal Processing. Wiley & Sons 1979, program 5.2.
*/
template <class U, class T>
std::vector<T> Core::Processing::Filter::CFIRfilter<U, T>::DesignLowPassFilter(int filterOrder, T normCutoffFreq)
{
	using namespace std;
	using namespace boost::math::constants;
	T sum, idealNormCutoffFreq;

	vector<T> impulseAnswer(filterOrder + 1);

	// convert from real low-pass cutoff frequency to the ideal one
	idealNormCutoffFreq = normCutoffFreq - static_cast<T>(3.3) / (filterOrder + 1);

	// input checks
	static_assert(((is_same< float, T >::value) || (is_same< double, T >::value)), "normCutoffFreq must be either of type float or double.");
	if ((filterOrder % 2) != 0) {
		throw std::length_error("Filter order must be even.");
	}
	if ((idealNormCutoffFreq <= 0) || (idealNormCutoffFreq > 1)) {
		throw std::logic_error("Normalized cutoff frequency must have a value of (0..1].");
	}

	// ideal low-pass filter
	int counter = 0;
	for (int n = -(filterOrder / 2); n <= (filterOrder / 2); n++) {
		if (n != 0) {
			impulseAnswer[counter] = sin(n * pi<T>() * idealNormCutoffFreq) / (n * pi<T>());
		} else {
			impulseAnswer[counter] = idealNormCutoffFreq;
		}
		counter++;
	}

	// apply Hamming window
	CDataProcessing<T>::HammingWindow(impulseAnswer.begin(), impulseAnswer.end(), impulseAnswer.begin());

	// prepare normalization for the center frequency of the passband (i.e. 0 Hz) being 0 dB
	sum = accumulate(impulseAnswer.begin(), impulseAnswer.end(), static_cast<T>(0));
	transform(impulseAnswer.begin(), impulseAnswer.end(), impulseAnswer.begin(), [=](T value) { return (value / sum); });

	// set output data
	return impulseAnswer;
}


/**	@brief		Calculates the filter response value for a certain frequency
*	@param		freq						Frequency for which the filter response value is calculated [Hz]. For frequencies larger than samplingFreq / 2.0, the value for samplingFreq / 2.0 is returned.
*	@param		model						Alglib-library interpolation model containing the filter response function
*	@return									Absolute value of the filter response function for the given frequency
*	@exception	std::range_error			Thrown if the input frequency is negative
*	@remarks								The function requires that an interpolation model containing the filter response function is existing. This allows a high performance.
*/
template <class U, class T>
T Core::Processing::Filter::CFIRfilter<U, T>::GetFreqResponse(const T& freq, const alglib::spline1dinterpolant& model)
{
	T currResponse;

	if (freq < 0) {
		throw std::range_error("Frequency is negative.");
	}

	// determine the required reponse for a certain frequency using the cubic Hermite interpolant
	currResponse = alglib::spline1dcalc(model, freq);

	return currResponse;
}


/**	@brief		Constructs an interpolation model for fast evaluation of the filter response function
*	@param		bFirst						Iterator to the beginning of the container with the filter coefficients of the FIR-filter
*	@param		bLast						Iterator past the end of the container with the filter coefficients of the FIR-filter
*	@param		deltaF						Required stepping for the evaluation of the filter reponse function [Hz]. Larger values increase speed.
*	@param		samplingFreq				Sampling frequency required for the filter [Hz]
*	@param		model						Contains the Alglib-library interpolation model after the function call
*	@return									None
*	@exception	std::range_error			Thrown if either the frequency stepping or the sampling frequency are negative
*	@exception	std::logic_error			Thrown if the frequency stepping is larger than half the sampling frequency
*	@exception	std::runtime_error			Thrown if the filter coefficient container is empty
*	@remarks								The function depends on the Alglib-library interpolation model
*/
template <class U, class T>
template <class InIt>
void Core::Processing::Filter::CFIRfilter<U, T>::GenerateFreqResponseModel(InIt bFirst, InIt bLast, const T& deltaF, const T& samplingFreq, alglib::spline1dinterpolant& model)
{
	using namespace std;
	long pointNum;
	vector<double> b;
	vector< std::complex<T> > gain;
	alglib::real_1d_array f, gainAbs;

	// obtain input values
	b.assign(bFirst, bLast);

	// check for wrong inputs
	if ((deltaF <= 0) || (samplingFreq <= 0)) {
		throw range_error("Frequency stepping or sampling frequency are not positive.");
	}
	if (deltaF > samplingFreq / 2.0) {
		throw logic_error("The frequency stepping is larger than half the sampling frequency.");
	}

	// calculate filter response function
	pointNum = static_cast<long>((samplingFreq / 2.0) / deltaF + 1);
	DigitalFilterGain(pointNum, samplingFreq, b.begin(), b.end(), back_inserter(gain));

	// generate model for interpolation	
	f.setlength(gain.size());
	gainAbs.setlength(gain.size());
	for (int i = 0; i < f.length(); i++) {
		f(i) = 0.0 + deltaF * i;
		gainAbs(i) = abs(gain[i]);
	}
	alglib::spline1dbuildmonotone(f, gainAbs, model);
}


/**	@brief		FIR-filter Hamming window based design for linearized phase low-pass filters finding the required filter order for a defined transition width of the filter
*	@param		transitionWidth				Width of the transition zone of the filter [Hz]. The transition zone is defined as the distance between the frequency where the response function goes below 0.975 and the 'cutoffFreq'.
*	@param		cutoffFreq					Cutoff frequency of the (real) low-pass filter [Hz]. It is the frequency where the main leaf of the filter has the minimum attenuation of the stop band (i.e. that of first side leaf).
*	@param		samplingFreq				Sampling frequency for which the filter has to be valid [Hz]
*	@param		startFilterOrder			Filter order used for starting the algorithm. If the start value is inappropriate, the resulting filter might be in a local minimum and thus wrong. The default value is 500.
*	@param		maxFilterOrder				Maximum filter order allowed. The default value is 1000.
*	@return									The b-filter container
*	@exception	std::length_error			Thrown if either the starting value for the filter order or the maximum filter order are negative or if the starting value for the filter order is larger than the maximum filter order
*	@exception	std::range_error			Thrown if the sampling frequency, the transition width or the cutoff frequency are negative, the transition width is larger than the cutoff frequency or if the cutoff frequency is larger than half of the sampling frequency
*	@exception	std::domain_error			Thrown if the maximum filter order is smaller than the minimum filter order required for realizing the given cutoff frequency
*	@remarks								The algorithm follows the IEEE Programs for Digital Signal Processing. Wiley & Sons 1979, program 5.2.
*/
template <class U, class T>
std::vector<T> Core::Processing::Filter::CFIRfilter<U, T>::DesignLowPassFilter(const T& transitionWidth, const T& cutoffFreq, const T& samplingFreq, const int& startFilterOrder, const int& maxFilterOrder)
{
	using namespace std;
	int finalFilterOrder;

	// check input
	if (samplingFreq <= 0) {
		throw std::range_error("Sampling frequency is not positive.");
	}
	if ((transitionWidth <= 0) || (cutoffFreq <= 0) || (transitionWidth > cutoffFreq)) {
		throw std::range_error("Transition width or cutoff frequency are negative of the transition width is larger than the cutoff frequency.");
	}
	if ((startFilterOrder <= 0) || (maxFilterOrder <= 0)) {
		throw std::length_error("The starting filter order or the maximum filter order are not positive.");
	}
	if (startFilterOrder > maxFilterOrder) {
		throw std::length_error("The starting filter order is larger than the maximum filter order");
	}

	// find the required transition width by solution of the non-linear system of transition width vs. filter order
	finalFilterOrder = Solver::FindTransitionWidth(startFilterOrder, transitionWidth, cutoffFreq, samplingFreq, maxFilterOrder);

	// generate the required filter
	return DesignLowPassFilter(finalFilterOrder, cutoffFreq / samplingFreq * static_cast<T>(2.0));
}


/**	@brief	Calculation of the response function of an arbitrary digital filter
*	@param		pointNumber					Number of points used for response function calculation. This refering to the number of points between [0, samplingFreq / 2.0]. The next larger power of two of (2 * pointNumber) must be at least as large as the filter order + 1.
*	@param		samplingFreq				Sampling frequency required for the filter [Hz]
*	@param		bFirst						First iterator of the container with digital filter parameters b
*	@param		bLast						Iterator past the last element of the container with the digital filter parameters b
*	@param		gainFirst					Iterator to the first element of the container, which will contain the complex response function between [0, samplingFreq / 2.0] after the function call. It must be of the required size ('pointNumber') or std::back_inserter must be used.
*	@return									Iterator past the last element of the response function container
*	@exception	std::runtime_error			Thrown if the filter coefficient container is empty
*	@exception	std::range_error			Thrown if the next larger power number of two of (2 * 'pointNumber') is smaller than the size of b (i.e. filter size + 1) or the sampling frequency 'samplingFreq' is negative or zero
*	@remarks								This function may not work correctly for filters which show strong changes near the frequency limits. This problem can be reduced by using datatype T = double.
*											The filter parameters can be obtained with MATLAB.
*/
template <class U, class T>
template <class InIt, class OutIt>
OutIt Core::Processing::Filter::CFIRfilter<U, T>::DigitalFilterGain(const long& pointNumber, const T& samplingFreq, InIt bFirst, InIt bLast, OutIt gainFirst)
{
	using namespace std;
	using namespace alglib;

	vector<T> a, b, f, amplitude, argument;
	vector< std::complex<T> > numerator, denominator, filterGain;
	long length;
	real_1d_array abscissa, pow2Amplitude, pow2Argument;
	spline1dinterpolant interpolationModelAmpl, interpolationModelArg;

	// get FIR-filter coefficients
	a.assign(1, 1);
	b.assign(bFirst, bLast);
	if (b.size() == 0) {
		throw runtime_error("The filter coefficients are empty.");
	}
	length = static_cast<long>(pow(2, ceil(log(2.0 * pointNumber) / log(2.0)))); // returns next larger power of two

	// check input data
	if (length < static_cast<int>(b.size())) {
		throw range_error("The filter length + 1 is smaller than the next larger power of two of (2 * 'pointNumber')");
	}
	if (samplingFreq <= 0) {
		throw std::range_error("Sampling frequency is not positive.");
	}

	// perform Fourier transformations for a and b filter coefficients (internally handled with powers of two for numerical stability)
	Core::Processing::CFFT<T> fft(length);
	f.resize(length);
	numerator.resize(length);
	denominator.resize(length / 2);
	denominator.insert(denominator.begin(), length / 2, 1); // simplification due to FIR-filter type
	fft.ComplexFFT(f.begin(), numerator.begin(), b.begin(), b.end(), samplingFreq);

	// adjust size (removes the unnecessary identical data in the second half)
	numerator.erase(numerator.begin() + numerator.size() / 2, numerator.end());
	denominator.erase(denominator.begin() + denominator.size() / 2, denominator.end());

	// calculate filter gain
	filterGain.resize(length / 2);
	transform(numerator.begin(), numerator.end(), denominator.begin(), filterGain.begin(), [](std::complex<T> val1, std::complex<T> val2) { return (val1 / val2); });

	// transform to amplitude and angle
	amplitude.resize(length / 2);
	argument.resize(length / 2);
	transform(filterGain.begin(), filterGain.end(), amplitude.begin(), [](std::complex<T> val) { return (std::abs<T>(val)); });
	transform(filterGain.begin(), filterGain.end(), argument.begin(), [](std::complex<T> val) { return (std::arg<T>(val)); });

	// interpolate filter gain to original length
	abscissa.setlength(filterGain.size());
	pow2Amplitude.setlength(filterGain.size());
	pow2Argument.setlength(filterGain.size());
	for (size_t i = 0; i < amplitude.size(); i++) {
		abscissa(i) = static_cast<double>(i);						// converted to double due to internal purposes
		pow2Amplitude(i) = static_cast<double>(amplitude[i]);		// converted to double due to internal purposes
		pow2Argument(i) = static_cast<double>(argument[i]);		// converted to double due to internal purposes
	}
	spline1dbuildcubic(abscissa, pow2Amplitude, interpolationModelAmpl);
	spline1dbuildcubic(abscissa, pow2Argument, interpolationModelArg);
	amplitude.resize(pointNumber);
	argument.resize(pointNumber);
	for (int i = 0; i < pointNumber; i++) { // conversion back to original datatype
		amplitude[i] = static_cast<T>(spline1dcalc(interpolationModelAmpl, i / static_cast<double>(pointNumber - 1)* (length / 2 - 1)));
		argument[i] = static_cast<T>(spline1dcalc(interpolationModelArg, i / static_cast<double>(pointNumber - 1)* (length / 2 - 1)));
	}

	// transform back to complex numbers
	filterGain.resize(amplitude.size());
	transform(amplitude.begin(), amplitude.end(), argument.begin(), filterGain.begin(), [](T ampl, T arg) { return (std::polar(ampl, arg)); });

	// store output data
	return std::move(filterGain.begin(), filterGain.end(), gainFirst);
}


/**	@brief	Calculation of the transition width of a certain filter
*	@param		filterOrder					Order of the filter (must be even)
*	@param		normCutoffFreq				Normalized cutoff frequency of the (real) low-pass filter. It is the frequency where the main leaf of the filter has the minimum attenuation of the stop band (i.e. that of first side leaf). The frequency is normalized with the Nyquist frequency, i.e. 0.5 * sampling frequency.
*	@param		samplingFreq				Sampling frequency [Hz]
*	@param		hLimit						Limit of the frequency response function defining the lower frequency border of the transition region. The default value is 0.975.
*	@param		deltaF						Frequency stepping for evaluating the frequency response function exactly, all other frequency are internally cubically interpolated. This evaluation is required for searching the transition width. Larger values reduce the exactness of the result. The default value is 20 Hz.
*	@param		maxResidual					Maximum residual accepted for searching the transition width. The default value is 1e-5.
*	@param		maxNumIt					Maximum number of iterations allowed for searching the transition width. The default value is 500.
*	@return									Width of the transition zone of the filter normalized with the sampling frequency. The transition zone is defined as the distance between the frequency where the response function goes below 'hLimit' to the real cutoff frequency'
*	@exception	std::range_error			Thrown if the sampling frequency is not positive or the normalized ideal cutoff frequency is not in the range of (0..1], the relation between the ideal and the real low pass filter is: fcNormIdeal = fcNormReal - 3.3 / ( filterOrder + 1 )
*	@exception	std::length_error			Thrown if the filter order is not even
*	@exception	std::logic_error			Thrown if 'deltaF' is larger than half the sampling frequency
*	@exception	std::runtime_error			Thrown if the sampling frequency is smaller than allowed by the frequency stepping 'deltaF'. The minimum sampling frequency can be achieved as follows: minPow2 = Next larger power of two of ( 2 * filterOrder ) / 4, then follows: minimum sampling frequency = 2 * deltaF * minPow2.
*	@remarks								This function may not work correctly for filters which show strong changes near the frequency limits. This problem can be reduced by using datatype T = double.
*/
template <class U, class T>
T Core::Processing::Filter::CFIRfilter<U, T>::CalcTransitionWidth(const int& filterOrder, const T& normCutoffFreq, const T& samplingFreq, const T& hLimit, const T& deltaF, const T& maxResidual, const int& maxNumIt)
{
	using namespace std;
	int numIt, minPow2;
	T freq, lowerFreq, upperFreq, response, transWidth, idealNormCutoffFreq, samplingFreqSmallest;
	alglib::spline1dinterpolant freqResponseModel;
	vector<T> b;

	// checking if the real cutoff frequency is valid (i.e. if the ideal cutoff frequency is positive)
	idealNormCutoffFreq = normCutoffFreq - 3.3 / (filterOrder + 1);
	if ((idealNormCutoffFreq <= 0) || (idealNormCutoffFreq > 1)) {
		throw std::range_error("The ideal cutoff frequency is not in the range [0,1).");
	}

	// checking for the minimum sampling frequency allowed
	minPow2 = static_cast<int>(pow(2, ceil(log(filterOrder + 1.0) / log(2.0))));
	minPow2 /= 2 * 2;
	samplingFreqSmallest = 2 * minPow2 * deltaF;
	if (abs(samplingFreq) < samplingFreqSmallest) {
		throw std::runtime_error("Sampling frequency is smaller than allowed by the frequency stepping.");
	}

	// determine the filter response function
	b = Core::Processing::Filter::CFIRfilter<U, T>::DesignLowPassFilter(filterOrder, normCutoffFreq);
	Core::Processing::Filter::CFIRfilter<U, T>::GenerateFreqResponseModel(b.begin(), b.end(), deltaF, samplingFreq, freqResponseModel);

	// determine the transition width by using a divide and conquer algorithm
	lowerFreq = 0;
	upperFreq = normCutoffFreq * samplingFreq / 2;
	numIt = 0;
	while ((abs(lowerFreq - upperFreq) > maxResidual) && (numIt < maxNumIt)) {
		freq = (lowerFreq + upperFreq) / 2;
		response = Core::Processing::Filter::CFIRfilter<U, T>::GetFreqResponse(freq, freqResponseModel);
		if (response > hLimit) {
			lowerFreq = freq;
		} else {
			upperFreq = freq;
		}
		numIt++;
	}
	transWidth = (normCutoffFreq * samplingFreq / 2.0 - freq) / samplingFreq;

	return transWidth;
}