#ifndef __ARGUMENTEXCEPTION_H
#define __ARGUMENTEXCEPTION_H

#include <string>
using std::string;

class ArgumentException {
	public:
		ArgumentException (const string& _message) : message(_message) {}

		const string& getMessage() const { return message; }

	private:

		const string message;
};

#endif
