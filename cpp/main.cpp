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
#include <future>
#include <vector>
#include "SDRReader.h"
#include "FMDemodulator.h"
#include "ComplexVector.h"
#include "AudioPlayer.h"

double stationFreq = 103.3e6; // in Hz
double samplingFreq = 2.304e6; // in Hz
const unsigned int sdrDeviceIndex = 0;
PaDeviceIndex audioDeviceIndex = 1;
const unsigned int timeoutFactor = 10;


void main() 
{
	using namespace std;

	try
	{
		ComplexVector<float> iqSamples;
		vector<float> audioSamples;
		future<ComplexVector<float>> iqSamplesFuture;
		future<vector<float>> audioSamplesFuture;

		FMDemodulatorSettings demodulatorSettings;
		FMDemodulator<float> fmDemodulator(stationFreq, samplingFreq, demodulatorSettings);

		CAudioPlayer<float> audioPlayer(audioDeviceIndex, fmDemodulator.GetAudioSamplingFreq());
		AudioSP::Audio::CSDRReader<float> sdrDevice(sdrDeviceIndex, fmDemodulator.GetCenterFreq(), samplingFreq, fmDemodulator.GetNumSignalSamplesForSDR(), 
			fmDemodulator.GetMaxNumSignalSamplesForProcessing());

		bool isFirst = true;
		while (true)
		{
			ThreadExceptionStorage::getInstance().throwStoredException();

			if (!isFirst)
			{
				auto status = iqSamplesFuture.wait_for(timeoutFactor * fmDemodulator.GetSampleTimePeriod());
				if (status == future_status::ready)
				{
					iqSamples = iqSamplesFuture.get();
				} else {
					iqSamples = ComplexVector<float>();
				}
			}

			iqSamplesFuture = sdrDevice.GetSignalData();

			if (!isFirst)
			{
				auto status = audioSamplesFuture.wait_for(timeoutFactor * fmDemodulator.GetSampleTimePeriod());
				if (status == future_status::ready)
				{
					audioSamples = audioSamplesFuture.get();
					audioPlayer.AddAudioData(move(audioSamples));
				}
			}

			audioSamplesFuture = fmDemodulator.Perform(move(iqSamples));

			isFirst = false;
		}
	}
	catch (const exception& e)
	{
		string errorMessage;
		try
		{
			ThreadExceptionStorage::getInstance().throwStoredException();
			errorMessage = e.what();
		}
		catch (const exception& e)
		{
			errorMessage = e.what();
		}

		cout << "ERROR: " << errorMessage << endl;
		exit(-1);
	}
}