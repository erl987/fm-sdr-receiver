#  FM-radio - software defined radio using RTL-SDR
#  Copyright (C) 2019-2020 Ralf Rettig
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.0)

project(FM-SDR-Receiver)

# global build system settings
if(EXISTS ${CMAKE_BINARY_DIR}/conan_paths.cmake)
	include(${CMAKE_BINARY_DIR}/conan_paths.cmake)
else()
	message(WARNING "The file conan_paths.cmake doesn't exist, you have to run conan install first")
endif()

set( CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeModules" ${CMAKE_MODULE_PATH} )

find_package( librtlsdr REQUIRED )
find_package( alglib REQUIRED )
find_package( Portaudio REQUIRED )
find_package( Boost COMPONENTS 
	system
REQUIRED )

set( SOURCE
	main.cpp
	AudioDevice.cpp
	SearchTransferFunc.cpp
)

set( HEADERS
	Polynomial.h
	Limiter.h
	AudioPlayer.h
	AudioDevice.h
	SDRReader.h
	FMDemodulator.h
	FIRfilter.h
	FilterUtils.h
	IIRfilter.h
	ComplexFIRfilter.h
	ComplexVector.h
	DataProcessing.h
	FFT.h
	SearchTransferFunc.h
	ThreadExceptionStorage.h
)

add_executable( ${PROJECT_NAME} ${SOURCE} ${HEADERS} )
target_link_libraries( ${PROJECT_NAME} 
	PRIVATE
		boost::boost
		portaudio::portaudio
		alglib::alglib
		librtlsdr::librtlsdr
)