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

template <typename T>
class Limiter
{
public:
	Limiter(T attackCoeff, T releaseCoeff, unsigned int numSamplesDelay, T threshold);
	std::vector<T> Processing(std::vector<T>&& signal);
private:
	size_t delayIndex;
	T envelope;
	T gain;
	unsigned int numSamplesDelay;
	std::vector<T> delayLine;
	T attackCoeff;
	T releaseCoeff;
	T threshold;
};

template <typename T>
Limiter<T>::Limiter(T attackCoeff, T releaseCoeff, unsigned int numSamplesDelay, T threshold)
	: delayIndex(0),
	  envelope(0),
	  gain(1),
	  numSamplesDelay(numSamplesDelay),
	  delayLine(numSamplesDelay),
	  releaseCoeff(releaseCoeff),
	  attackCoeff(attackCoeff),
	  threshold(threshold)
{
}

template <typename T>
std::vector<T> Limiter<T>::Processing(std::vector<T>&& signal)
{
	T targetGain;
	std::vector<T> limitedSignal(signal.size());

	for (size_t i = 0; i < signal.size(); i++)
	{
		delayLine[delayIndex] = signal[i];
		delayIndex = (delayIndex + 1) % numSamplesDelay;

		envelope *= releaseCoeff;
		envelope = std::max(abs(signal[i]), envelope);

		if (envelope > threshold)
		{
			targetGain = 1 + threshold - envelope;
		}
		else
		{
			targetGain = 1.0;
		}

		gain = gain * attackCoeff + targetGain * (1 - attackCoeff);
		limitedSignal[i] = delayLine[delayIndex] * gain;
	}

	return limitedSignal;
}
