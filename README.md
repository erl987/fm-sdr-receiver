# Realtime Software Defined FM-Radio implemented in different languages

This repository contains implementations for software define FM-radios in different languages. The programs are using a
simple DVB-T USB-stick for receiving the radio signals via RTL-SDR.

# Implementations

* Python
* C++


# Prerequisites

* DVB-T USB-stick compatible to RTL-SDR (see https://www.rtl-sdr.com/ for details)
* RTL-SDR installed on the machine (this requires custom driver installation, see https://www.rtl-sdr.com/ for details)


# How to get started

## Python

The realtime implementation is started by running the file `realtime.py`. All settings have to be
changed within the file. **Usually it is sufficient to adjust the frequency of the station.**

A non-realtime implementation is also available in the file `non_realtime.py`.


## C++

This is a realtime implementation that runs on low processor load. **The station frequency is hard-coded
into the file `main.cpp`.**

The program uses `conan`-dependency management, but not all dependencies are available on public `conan`-repositories.
It is required to download and build them manually and possibly to provide the via `conan`.


### Dependencies

* `librtlsdr`
* `boost`
* `alglib`
* `portaudio`


### Building

The program can be build with `cmake`.


# Performance


Both, the Python and the C++ realtime implementations have been performance-optimized quite intensively.
Even the Python-implementation runs at a processor load of ca. 7% on a typical desktop PC. The C++ implementation
runs with ca. 2% processor load.


# License

FM-radio - software defined radio using RTL-SDR
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
