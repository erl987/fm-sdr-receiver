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
#include <vector>

template <typename T>
class ComplexVector {
public:
	ComplexVector() = default;
	ComplexVector(size_t size);
	ComplexVector(std::vector<T>&& real, std::vector <T>&& imag);
	ComplexVector(const ComplexVector<T>& rhs);
	ComplexVector(const ComplexVector<T>& rhs, size_t length);
	ComplexVector(ComplexVector<T>&& rhs) noexcept;
	void resize(size_t newSize);
	size_t size() const;
	std::complex<T> at(size_t pos) const;
	void conj();
	ComplexVector<T>& multiplySubset(const ComplexVector<T>& rhs, size_t startIndexLhs, size_t endIndexLhs, size_t startIndexRhs);
	ComplexVector<T>& operator*=(const ComplexVector<T>& rhs);
	ComplexVector<T>& operator=(ComplexVector<T>&& rhs);
	template<typename U> friend const ComplexVector<T> operator*(const ComplexVector<U>& lhs, const ComplexVector<U>& rhs);
	template<typename U> friend const ComplexVector<T> operator*(const std::complex<U>& lhs, const ComplexVector<U>& rhs);

	const std::vector<T>& real() const;
	const std::vector<T>& imag() const;
	std::vector<T>& realRef();
	std::vector<T>& imagRef();
private:
	std::vector<T> realVal;
	std::vector<T> imagVal;
};

template <typename T>
ComplexVector<T>::ComplexVector(size_t size)
	: realVal(size),
      imagVal(size)
{
}

template <typename T>
ComplexVector<T>::ComplexVector(std::vector<T>&& real, std::vector<T>&& imag)
	: realVal(std::move(real)),
	  imagVal(std::move(imag))
{
}

template <typename T>
ComplexVector<T>::ComplexVector(const ComplexVector<T>& rhs)
	: realVal(rhs.realVal),
	  imagVal(rhs.imagVal)
{
}

template <typename T>
ComplexVector<T>::ComplexVector(const ComplexVector<T>& rhs, size_t length)
{
	assert(length <= rhs.realVal.size());
	assert(length <= rhs.imagVal.size());
	realVal.assign(rhs.realVal.begin(), rhs.realVal.begin() + length);
	imagVal.assign(rhs.imagVal.begin(), rhs.imagVal.begin() + length);
}

template <typename T>
ComplexVector<T>::ComplexVector(ComplexVector<T>&& rhs) noexcept
	: realVal(std::move(rhs.realVal)),
	  imagVal(std::move(rhs.imagVal))
{
}

template<typename T>
const std::vector<T>& ComplexVector<T>::real() const
{
	return realVal;
}

template<typename T>
const std::vector<T>& ComplexVector<T>::imag() const
{
	return imagVal;
}

template <typename T>
std::vector<T>& ComplexVector<T>::realRef()
{
	return realVal;
}

template <typename T>
std::vector<T>& ComplexVector<T>::imagRef()
{
	return imagVal;
}

template <typename T>
size_t ComplexVector<T>::size() const
{
	return realVal.size();
}

template<typename T>
void ComplexVector<T>::resize(size_t newSize)
{
	realVal.resize(newSize);
	imagVal.resize(newSize);
}

template <typename T>
ComplexVector<T>& ComplexVector<T>::operator*=(const ComplexVector<T>& rhs)
{
	using namespace std;

	for (int i = 0; i < rhs.size(); i++) 
	{
		T newReal = real()[i] * rhs.real()[i] - imag()[i] * rhs.imag()[i];
		T newImag = real()[i] * rhs.imag()[i] + imag()[i] * rhs.real()[i];
		realVal[i] = newReal;
		imagVal[i] = newImag;
	}

	return *this;
}

template <typename U>
const ComplexVector<U> operator*(const ComplexVector<U>& lhs, const ComplexVector<U>& rhs)
{
	using namespace std;

	vector<U> multipliedReal(lhs.size());
	vector<U> multipliedImag(lhs.size());

	for (int i = 0; i < lhs.size(); i++) 
	{
		multipliedReal[i] = lhs.real()[i] * rhs.real()[i] - lhs.imag()[i] * rhs.imag()[i];
		multipliedImag[i] = lhs.real()[i] * rhs.imag()[i] + lhs.imag()[i] * rhs.real()[i];
	}

	return ComplexVector<U>(move(multipliedReal), move(multipliedImag));
}

template<typename U>
const ComplexVector<U> operator*(const std::complex<U>& lhs, const ComplexVector<U>& rhs)
{
	using namespace std;

	vector<U> multipliedReal(rhs.size());
	vector<U> multipliedImag(rhs.size());

	for (int i = 0; i < rhs.size(); i++) 
	{
		multipliedReal[i] = lhs.real() * rhs.real()[i] - lhs.imag() * rhs.imag()[i];
		multipliedImag[i] = lhs.real() * rhs.imag()[i] + lhs.imag() * rhs.real()[i];
	}

	return ComplexVector<U>(move(multipliedReal), move(multipliedImag));
}

template <typename T>
ComplexVector<T>& ComplexVector<T>::operator=(ComplexVector<T>&& rhs)
{
	if (this != &rhs) 
	{
		realVal = std::move(rhs.real());
		imagVal = std::move(rhs.imag());
	}
	return *this;
}

template <typename T>
void ComplexVector<T>::conj()
{
	for (int i = 0; i < imagVal.size(); i++) 
	{
		imagVal[i] = -imagVal[i];
	}
}

template <typename T>
std::vector<T> arg(const ComplexVector<T>& numbers)
{
	std::vector<T> arg(numbers.size());
	for (int i = 0; i < numbers.size(); i++) 
	{
		arg[i] = std::atan2(numbers.imag()[i], numbers.real()[i]);
	}
	return arg;
}

template <typename T>
ComplexVector<T>& ComplexVector<T>::multiplySubset(const ComplexVector<T>& rhs, size_t startIndexLhs, size_t endIndexLhs, size_t startIndexRhs)
{
	size_t length = endIndexLhs - startIndexLhs;
	assert(startIndexLhs + length < size() && startIndexRhs + length < rhs.size());

	for (int i = 0; i < length; i++)
	{
		size_t iLhs = i + startIndexLhs;
		size_t iRhs = i + startIndexRhs;
		T newReal = real()[iLhs] * rhs.real()[iRhs] - imag()[iLhs] * rhs.imag()[iRhs];
		T newImag = real()[iLhs] * rhs.imag()[iRhs] + imag()[iLhs] * rhs.real()[iRhs];
		realVal[i] = newReal;
		imagVal[i] = newImag;
	}

	realVal.resize(length);
	imagVal.resize(length);

	return *this;
}

template <typename T>
std::vector<T> polarDiscriminator(ComplexVector<T>&& signal, std::complex<T>&& lastBandwidthSample)
{
	ComplexVector<T> conjSignal = signal;
	conjSignal.conj();
	
	T firstChannelSample = std::arg(signal.at(0) * std::conj(lastBandwidthSample));

	std::vector<T> result(signal.size());
	result[0] = firstChannelSample;
	std::vector<T> tempResult = arg(signal.multiplySubset(conjSignal, 1, conjSignal.size(), 0));
	for (int i=0; i < tempResult.size(); i++) {
		result[i+1] = tempResult[i];
	}
	return result;
}


template <typename T>
std::complex<T> ComplexVector<T>::at(size_t pos) const
{
	return std::complex<T>(realVal[pos], imagVal[pos]);
}