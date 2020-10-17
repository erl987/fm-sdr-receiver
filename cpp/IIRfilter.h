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
#include <algorithm>
#include <complex>
#include <boost/math/constants/constants.hpp>
#include "FilterUtils.h"
#include "Polynomial.h"


/*@{*/
/** \ingroup Core
*/
namespace Core {
	namespace Processing {
		namespace Filter {
			/** \ingroup Core
			*	Class for filtering data with digital IIR-filters in the time domain. The class is especially well-suited for processing arbitrary length input signals.
			*/
			template <typename U, typename T> class CIIRfilter
			{
			public:
				CIIRfilter(void);
				~CIIRfilter(void) {};
				template <typename InIt1, typename InIt2> CIIRfilter(InIt1 aFirst, InIt1 aLast, InIt2 bFirst, InIt2 bLast, int downsamplingFactor = 1, int upsamplingFactor = 1);
				template <typename InIt1, typename InIt2> void SetParams(InIt1 aFirst, InIt1 aLast, InIt2 bFirst, InIt2 bLast, int downsamplingFactor = 1, int upsamplingFactor = 1);
				std::vector<T> Processing(std::vector<T>&& signal);
				std::pair<std::vector<U>, std::vector<T>> Processing(std::vector<U>&& time, std::vector<T>&& signal);
				static std::pair<std::vector<T>, std::vector<T>> DesignHighPassFilter(const int& filterOrder, const double& cutoffFreq, const double& samplingFreq);
			private:
				std::vector<T> DownsamplingFilter(std::vector<T>&& inputSignal);
				CIIRfilter(const CIIRfilter&) = delete;
				CIIRfilter& operator= (const CIIRfilter&) = delete;
				template <typename V> static std::vector<V> poly(const std::vector<V>& roots);

				CFilterUtils<U, T> utils;
				std::vector<T> a;
				std::vector<T> previousFilteredSignal;
				bool isInit;
			};
		}
	}
}

/*@}*/


/**	@brief		Default constructor
*/
template <typename U, typename T>
Core::Processing::Filter::CIIRfilter<U, T>::CIIRfilter(void)
	: isInit(false)
{
}


/**	@brief 		Constructor
*	@param		aFirst					Iterator to the beginning of the container storing the a-filter parameters of the IIR-filter
*	@param		aLast					Iterator to the end of the container storing the a-filter parameters
*	@param		bFirst					Iterator to the beginning of the container storing the b-filter parameters of the IIR-filter
*	@param		bLast					Iterator to the end of the container storing the b-filter parameters
*	@param		downsamplingFactor		Downsampling factor: every downsamplingFactor-th datapoint is used. Can be leaved out together with upsamplingFactor if only filtering is used.
*	@param		upsamplingFactor		Upsampling factor: diving the section between two datapoints in upsamplingFactor sections. Can be leaved out if only filtering is used.
*	@exception	std::runtime_error		Thrown if the downsampling or upsampling factor is smaller than 1.
*	@remarks							Upsampling and downsampling factors are automatically minimized by the greatest common divisor. If it is used for downsampling, the filter must be correctly designed regarding the cut-off frequency at fc = 0.8 * ( upsamplingFactor / downsamplingFactor * sampling rate / 2 ).
*										The algorithm follows the IEEE Programs for Digital Signal Processing. Wiley & Sons 1979, program 8.2. A 8-order IIR lowpass Chebyshev type I filter with passband ripples smaller than 0.05 dB is recommended for downsampling.
*/
template <typename U, typename T>
template <typename InIt1, typename InIt2>
Core::Processing::Filter::CIIRfilter<U, T>::CIIRfilter(InIt1 aFirst, InIt1 aLast, InIt2 bFirst, InIt2 bLast, int downsamplingFactor, int upsamplingFactor)
	: isInit(false)
{
	SetParams(aFirst, aLast, bFirst, bLast, downsamplingFactor, upsamplingFactor);
}


/**	@brief 		Safely resetting filter parameters (also possible for an already used object)
*	@param		aFirst					Iterator to the beginning of the container storing the a-filter parameters of the IIR-filter
*	@param		aLast					Iterator to the end of the container storing the a-filter parameters
*	@param		bFirst					Iterator to the beginning of the container storing the b-filter parameters of the IIR-filter
*	@param		bLast					Iterator to the end of the container storing the b-filter parameters
*	@param		downsamplingFactor		Downsampling factor: every downsamplingFactor-th datapoint is used. Can be leaved out together with upsamplingFactor if only filtering is used.
*	@param		upsamplingFactor		Upsampling factor: diving the section between two datapoints in upsamplingFactor sections. Can be leaved out if only filtering is used.
*	@return								None
*	@exception	std::length_error		Thrown if the filter parameters are empty
*	@exception	std::runtime_error		Thrown if the downsampling or upsampling factor is smaller than 1.
*	@remarks							Upsampling and downsampling factors are automatically minimized by the greatest common divisor. If it is used for downsampling, the filter must be correctly designed regarding the cut-off frequency at fc = 0.8 * ( upsamplingFactor / downsamplingFactor * sampling rate / 2 ).
*										The algorithm follows the IEEE Programs for Digital Signal Processing. Wiley & Sons 1979, program 8.2. A 8-order IIR lowpass Chebyshev type I filter with passband ripples smaller than 0.05 dB is recommended for downsampling.
*/
template <typename U, typename T>
template <typename InIt1, typename InIt2>
void Core::Processing::Filter::CIIRfilter<U, T>::SetParams(InIt1 aFirst, InIt1 aLast, InIt2 bFirst, InIt2 bLast, int downsamplingFactor, int upsamplingFactor)
{
	using namespace std;

	utils.SetParams(bFirst, bLast, downsamplingFactor, upsamplingFactor);

	// assign input values
	a.assign(aFirst, aLast);
	if (a.empty()) {
		throw std::length_error("Filter parameters must not be empty.");
	}

	// set storage container for recently processed data (required for restarting the filtering without losses)
	previousFilteredSignal.assign(a.size() - 1, 0);

	isInit = true;
}


/**	@brief 		Performs the filtering of the data.
*	@param		signal					The data to be filtered
*	@return								The filtered data
*	@exception							None
*	@remarks							None
*/
template <typename U, typename T>
std::vector<T> Core::Processing::Filter::CIIRfilter<U, T>::Processing(std::vector<T>&& signal)
{
	using namespace std::placeholders;
	return utils.Processing(move(signal), std::bind(&CIIRfilter<U, T>::DownsamplingFilter, this, _1));
}


/**	@brief 		Performs the filtering of the data that also contains time-information.
*	@param		time					The timepoints of the data to be filtered
*	@param		signal					The data to be filtered
*	@return								A pair of the filtered time and data containers
*	@exception							None
*	@remarks							None
*/
template <typename U, typename T>
std::pair<std::vector<U>, std::vector<T>> Core::Processing::Filter::CIIRfilter<U, T>::Processing(std::vector<U>&& time, std::vector<T>&& signal)
{
	return utils.Processing(move(time), move(signal), std::bind(&CIIRfilter<T>::DownsamplingFilter, this, std::placeholders::_1));
}


/** @brief		Downsampling and filtering of a dataset. This function is specifically implemented for each filter type.
*	@param		signal					Data to be processed
*	@return								Downsampled output data
*	@exception							None
*	@remarks							This function performs efficiently filtering and downsampling at the same time
*/
template <typename U, typename T>
std::vector<T> Core::Processing::Filter::CIIRfilter<U, T>::DownsamplingFilter(std::vector<T>&& inputSignal)
{
	size_t n, nFiltered, signalSize;
	std::vector<T> downsampledSignal, a, b, filteredSignal;
	int startIndexForA, startIndexForB;
	T filtered;
	std::vector<T> signal;

	signalSize = inputSignal.size();
	signal.reserve(inputSignal.size() + utils.GetPreviousSignalRef().size());
	signal.assign(inputSignal.begin(), inputSignal.end());

	a.assign(CIIRfilter::a.begin(), CIIRfilter::a.end());
	b.assign(utils.GetBRef().begin(), utils.GetBRef().end());

	if (a[0] != 1) {
		throw std::runtime_error("Filter parameter a[0]=1 is required, instead it is a[0]=" + std::to_string(a[0]));
	}

	// for better loop performance
	std::reverse(a.begin(), a.end());
	std::reverse(b.begin(), b.end());

	// perform filtering
	filteredSignal.reserve(signalSize + this->previousFilteredSignal.size());
	filteredSignal.resize(signalSize);
	startIndexForA = static_cast<int>(a.size()) - 1;
	startIndexForB = static_cast<int>(b.size()) - 1;

	// add previous data in front of the current data
	signal.insert(signal.begin(), utils.GetPreviousSignalRef().begin(), utils.GetPreviousSignalRef().end());
	filteredSignal.insert(filteredSignal.begin(), this->previousFilteredSignal.begin(), this->previousFilteredSignal.end());

	for (size_t index = 0; index < signalSize; index++) {
		nFiltered = index + this->previousFilteredSignal.size();
		n = index + utils.GetPreviousSignalRef().size();

		filtered = 0;
		for (int k = startIndexForB; k >= 0; k--) {
			filtered += b[startIndexForB - k] * signal[n - k];
		}
		for (int l = startIndexForA; l >= 1; l--) {
			filtered -= a[startIndexForA - l] * filteredSignal[nFiltered - l];
		}
		filteredSignal[nFiltered] = filtered;
	}

	// delete temporary previous data
	signal.erase(signal.begin(), signal.begin() + utils.GetPreviousSignalRef().size());
	filteredSignal.erase(filteredSignal.begin(), filteredSignal.begin() + previousFilteredSignal.size());

	// preserve the latest data corresponding to the length of the filter (for restarting the filter)
	if (signal.size() > startIndexForB) {
		utils.GetPreviousSignalRef().assign(signal.end() - startIndexForB, signal.end());
	} else {
		utils.GetPreviousSignalRef().assign(signal.begin(), signal.end());
	}

	if (filteredSignal.size() > startIndexForA) {
		this->previousFilteredSignal.assign(filteredSignal.end() - startIndexForA, filteredSignal.end());
	} else {
		this->previousFilteredSignal.assign(filteredSignal.begin(), filteredSignal.end());
	}

	// downsampling
	if (utils.GetDownsamplingFactor() > 1)
	{
		filteredSignal = utils.Downsampling(move(filteredSignal));
	}

	return filteredSignal;
}


/** @brief		Determines the coefficients of a polynomial if its roots are given.
*	@param		root					The roots of the polynomial
*	@return								The coefficients of the polynomial, from the highest to the lowest (the highest is always 1)
*	@exception							None
*	@remarks							None
*/
template <typename U, typename T>
template <typename V>
std::vector<V> Core::Processing::Filter::CIIRfilter<U, T>::poly(const std::vector<V>& roots) {
	Polynomial<V> p = Polynomial<V>(std::vector<V>{1.0});

	for (const auto& root : roots) {
		Polynomial<V> pClone = p.clone();
		p.multiplyByX();
		pClone.multiply(root);
		p.subtract(pClone);
	}

	return p.getCoefficients();
}


/** @brief		Find the coefficients of an IIR Butterworth high-pass filter using bilinear transform.
*	@param		filterOrder				Filter order
*	@param		cutoffFreq				Cutoff frequency of the filter (-3dB frequency) [Hz]
*   @param		samplingFreq			Sampling frequency of the data to be filtered [Hz]
*	@return								Numerator, denominator coefficients of the digital filter (b, a)
*	@exception	std::logic_error		If the cutoff frequency is not less than 0.5 of the sampling frequency
*	@remarks							None
*/
template <typename U, typename T>
std::pair<std::vector<T>, std::vector<T>> Core::Processing::Filter::CIIRfilter<U, T>::DesignHighPassFilter(const int& filterOrder, const double& cutoffFreq, const double& samplingFreq)
{
	using namespace std;
	using namespace std::complex_literals;
	using namespace boost::math::constants;

	double F_c, kFac;
	vector<double> k(filterOrder);
	vector<double> theta(filterOrder);
	vector<complex<double>> pa(filterOrder);
	vector<complex<double>> p(filterOrder);
	vector<double> q(filterOrder);
	vector<double> m(filterOrder + 1);
	vector<double> term1(filterOrder + 1);
	vector<double> term2(filterOrder + 1);
	vector<complex<double>> aComplex(filterOrder + 1);
	vector<double> bUnscaled(filterOrder + 1);
	vector<T> a(filterOrder + 1);
	vector<T> b(filterOrder + 1);


	if (cutoffFreq / samplingFreq >= 0.5) {
		throw logic_error("The normalized cutoff frequency must be less than 0.5");
	}

	// find poles of the normalized analog filter
	iota(k.begin(), k.end(), 1.0);
	transform(k.begin(), k.end(), theta.begin(), [=](auto k) {return ((2 * k - 1) * pi<double>() / (2 * filterOrder)); });
	transform(theta.begin(), theta.end(), pa.begin(), [](auto theta) { return (-sin(theta) + 1i * cos(theta)); }); // poles of low-pass filter with cutoff = 1 rad/s

	// transform poles for high-pass filtering
	F_c = samplingFreq / pi<double>() * tan(pi<double>() * cutoffFreq / samplingFreq);  // continuous pre-warped frequency
	transform(pa.begin(), pa.end(), pa.begin(), [=](auto pa) { return (2 * pi<double>() * F_c / pa); }); // analog high-pass poles

	// find coefficients of the digital filter
	// poles and zeros in the z-plane
	transform(pa.begin(), pa.end(), p.begin(), [=](auto pa) { return (1.0 + pa / (2 * samplingFreq)) / (1.0 - pa / (2 * samplingFreq)); });
	fill(q.begin(), q.end(), 1.0); // zeros at z=1 (f=0)

	// convert poles and zeros to polynomial coefficients
	aComplex = poly(p);
	transform(aComplex.begin(), aComplex.end(), a.begin(), [](auto a) { return static_cast<T>(a.real()); });

	// amplitude scale factor for gain=1 at f=fs/2 (z=-1)
	bUnscaled = poly(q);
	iota(m.begin(), m.end(), 0.0);
	transform(a.begin(), a.end(), m.begin(), term1.begin(), [](auto a, auto m) { return (pow(-1.0, m) * a); });
	transform(bUnscaled.begin(), bUnscaled.end(), m.begin(), term2.begin(), [](auto b, auto m) { return (pow(-1.0, m) * b); });

	kFac = accumulate(term1.begin(), term1.end(), 0.0) / accumulate(term2.begin(), term2.end(), 0.0);
	transform(bUnscaled.begin(), bUnscaled.end(), b.begin(), [=](auto b) { return static_cast<T>(kFac * b); });

	return make_pair<>(b, a);
}