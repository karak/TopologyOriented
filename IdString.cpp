#include "StdAfx.h"
#include "IdString.h"
#include <cstdio>
#include <ostream>
#include <algorithm>


IdString::IdString(char c, int idValue)
{
	::sprintf_s(buff,"%c%0*d", c,NUM_COLS, idValue);
}

IdString::IdString(const IdString& rhs)
{
	std::copy(rhs.buff, rhs.buff + BUFFER_SIZE, buff);
}

std::ostream& operator<<(std::ostream& os, const IdString& id)
{
	return os << id.c_str();
}