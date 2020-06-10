#ifndef __FILEEXCEPTION_H
#define __FILEEXCEPTION_H

#include <stdexcept>
using std::runtime_error;
#include <string>
using std::string;

class FileException : public runtime_error {
	public:
		FileException (const string& msg) : runtime_error(msg) {}
};

#endif
