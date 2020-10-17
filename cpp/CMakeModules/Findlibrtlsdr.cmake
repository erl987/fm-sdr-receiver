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

find_path(
	librtlsdr_DIR
	NAMES rtl-sdr.h
	HINTS
		/usr/include
		/usr/local/include
		ENV RTL_SDR_ROOT
		${PROJECT_SOURCE_DIR}/libraries/librtlsdr/include
)

set( librtlsdr_INCLUDE_DIR ${librtlsdr_DIR} )

# find the libraries
find_library(
	librtlsdr_LIBRARY
	NAMES rtlsdr
	PATHS 
		/usr/lib/x86_64-linux-gnu
		/usr/lib
		/usr/local/lib
		${librtlsdr_DIR}/../build/src/Release
)

# check that all files and directories were properly found
include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( 
	librtlsdr
	DEFAULT_MSG
	librtlsdr_INCLUDE_DIR
	librtlsdr_LIBRARY
)

if (librtlsdr_FOUND AND NOT TARGET librtlsdr::librtlsdr)
	add_library( librtlsdr::librtlsdr UNKNOWN IMPORTED )
	set_target_properties( librtlsdr::librtlsdr PROPERTIES
			IMPORTED_LOCATION "${librtlsdr_LIBRARY}"
			INTERFACE_INCLUDE_DIRECTORIES "${librtlsdr_INCLUDE_DIR}" )
endif()