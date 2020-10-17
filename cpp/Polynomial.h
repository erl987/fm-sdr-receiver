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

/*@{*/
/** \ingroup Core
*/
namespace Core {
	namespace Processing {
		namespace Filter {
			/** \ingroup Core
			*	Representation of a polynomial.
			*/
			template <typename T>
			class Polynomial {
			public:
				Polynomial(const std::vector<T>& coefficientsFromLowest);
				std::vector<T> getCoefficients() const;
				Polynomial clone();
				void multiply(const T& factor);
				void multiplyByX();
				void subtract(const Polynomial& other);
				size_t getOrder() const;
			private:
				std::vector<T> getCoefficientsFromLowest() const;

				std::vector<T> coefficientsFromLowest;
			};
		}
	}
}

/*@}*/


/**	@brief 		Constructor
*	@param		coefficientsFromLowest	The coefficients of the polynomial
*	@exception							None
*	@remarks							The order of the coefficients is from the lowest to the highest!
*/
template <typename T>
Core::Processing::Filter::Polynomial<T>::Polynomial(const std::vector<T>& coefficientsFromLowest)
	: coefficientsFromLowest(coefficientsFromLowest)
{
}


/**	@brief 		Returns the coefficients of the polynomial from the lowest to the highest
*	@return								The coefficients of the polynomial
*	@exception							None
*	@remarks							None
*/
template <typename T>
std::vector<T> Core::Processing::Filter::Polynomial<T>::getCoefficientsFromLowest() const
{
	return coefficientsFromLowest;
}


/**	@brief 		Returns the order of the polynomial
*	@return								The order of the polynomial
*	@exception							None
*	@remarks							None
*/
template <typename T>
size_t Core::Processing::Filter::Polynomial<T>::getOrder() const
{
	return coefficientsFromLowest.size();
}


/**	@brief 		Returns the coefficients of the polynomial from the highest to the lowest
*	@return								The coefficients of the polynomial
*	@exception							None
*	@remarks							None
*/
template <typename T>
std::vector<T> Core::Processing::Filter::Polynomial<T>::getCoefficients() const
{
	std::vector<T> coefficientsFromHighest = coefficientsFromLowest;
	std::reverse(coefficientsFromHighest.begin(), coefficientsFromHighest.end());
	return coefficientsFromHighest;
}


/**	@brief 		Clones the polynomial object
*	@return								The cloned polynomial object
*	@exception							None
*	@remarks							None
*/
template <typename T>
Core::Processing::Filter::Polynomial<T> Core::Processing::Filter::Polynomial<T>::clone()
{
	return Polynomial(coefficientsFromLowest);
}


/**	@brief 		Multiplies the whole polynomial by a factor
*	@param		factor					The factor
*	@return								None
*	@exception							None
*	@remarks							The multiplication is equivalent to for example: `factor * (1 * x^2 + 5 * x - 3)`
*/
template <typename T>
void Core::Processing::Filter::Polynomial<T>::multiply(const T& factor)
{
	std::transform(coefficientsFromLowest.begin(), coefficientsFromLowest.end(), coefficientsFromLowest.begin(), [=](auto val) { return (factor * val); });
}


/**	@brief 		Multiplies the whole polynomial by `x`, i.e. this increases its order by 1
*	@return								None
*	@exception							None
*	@remarks							The multiplication is equivalent to for example: `x * (1 * x^2 + 5 * x - 3) = 1 * x^3 + 5 * x^2 - 3*x + 0`
*/
template <typename T>
void Core::Processing::Filter::Polynomial<T>::multiplyByX()
{
	coefficientsFromLowest.insert(coefficientsFromLowest.begin(), 0.0);
}


/**	@brief 		Subtracts another polynomial from this polynomial
*	@param		other					The other polynomial
*	@return								None
*	@exception							None
*	@remarks							The subtraction is equivalent to for example: `(1 * x^2 + 5 * x - 3) - (5*x - 9)`
*/
template <typename T>
void Core::Processing::Filter::Polynomial<T>::subtract(const Polynomial<T>& other)
{
	size_t minOrder = std::min(getOrder(), other.getOrder());
	for (int i = 0; i < minOrder; i++) {
		coefficientsFromLowest[i] -= other.getCoefficientsFromLowest()[i];
	}

	if (minOrder == getOrder()) {
		std::vector<T> coefficientsToAdd = other.getCoefficientsFromLowest();
		coefficientsToAdd.erase(coefficientsToAdd.begin(), coefficientsToAdd.begin() + minOrder + 1);
		coefficientsFromLowest.insert(coefficientsFromLowest.end(), coefficientsToAdd.begin(), coefficientsToAdd.end());
	}
}