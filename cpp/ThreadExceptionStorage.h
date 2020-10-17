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
#include <exception>
#include <iostream>

class ThreadExceptionStorage {
public:
	static ThreadExceptionStorage& getInstance();
	void storeException(std::exception_ptr exception);
	void throwStoredException();
private:
	ThreadExceptionStorage() = default;
	~ThreadExceptionStorage() = default;
	ThreadExceptionStorage(const ThreadExceptionStorage&) = delete;
	ThreadExceptionStorage& operator=(const ThreadExceptionStorage&) = delete;

	std::exception_ptr storedException;
};


ThreadExceptionStorage& ThreadExceptionStorage::getInstance()
{
	static ThreadExceptionStorage instance;
	return instance;
}

void ThreadExceptionStorage::storeException(std::exception_ptr exception)
{
	// no overriding of already stored exceptions
	if (storedException == nullptr) {
		storedException = exception;
	}
};

void ThreadExceptionStorage::throwStoredException()
{
	if (storedException)
	{
		std::rethrow_exception(storedException);
	}
};