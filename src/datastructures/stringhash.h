#ifndef __STRINGHASH_H
#define __STRINGHASH_H

//#include "../headers.h"

/**
 * hasher for std::string
 */

namespace __gnu_cxx {

struct stringhash
{
	size_t operator()(const string& s) const
	{ return __stl_hash_string(s.c_str()); }
};

}

#endif
