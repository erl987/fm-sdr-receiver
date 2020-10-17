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
#include <complex>
#include "FIRfilter.h"

namespace Core {
	namespace Processing {
		namespace Filter {
			template <typename T>
			class ComplexFIRfilter
			{
			public:
				ComplexFIRfilter(const unsigned int& downsamplingFactor, const unsigned int& upsamplingFactor, const unsigned int& filterOrder=30);
				ComplexVector<T> Processing(ComplexVector<T>&& signal);
				unsigned int ProcessedLength(int dataLength) const;
			private:
				CFIRfilter<T, T> realPartFilter;
				CFIRfilter<T, T> imagPartFilter;
			};
		}
	}
}

template <typename T>
Core::Processing::Filter::ComplexFIRfilter<T>::ComplexFIRfilter(const unsigned int& downsamplingFactor, const unsigned int& upsamplingFactor, const unsigned int& filterOrder)
	: realPartFilter(static_cast<int>(downsamplingFactor), static_cast<int>(upsamplingFactor), static_cast<int>(filterOrder), static_cast<T>(1e-5)),
	  imagPartFilter(static_cast<int>(downsamplingFactor), static_cast<int>(upsamplingFactor), static_cast<int>(filterOrder), static_cast<T>(1e-5))
{
}

template <typename T>
ComplexVector<T> Core::Processing::Filter::ComplexFIRfilter<T>::Processing(ComplexVector<T>&& signal)
{
	using namespace std;

	// perform the filtering separately for the real and imaginary part (assumes a linear filter)
	vector<T> realFilteredSignal= realPartFilter.Processing(move(signal.realRef()));
	vector<T> imagFilteredSignal = imagPartFilter.Processing(move(signal.imagRef()));

	return ComplexVector<T>(move(realFilteredSignal), move(imagFilteredSignal));
}

template <typename T>
unsigned int Core::Processing::Filter::ComplexFIRfilter<T>::ProcessedLength(int dataLength) const
{
	return realPartFilter.ProcessedLength(dataLength);
}
