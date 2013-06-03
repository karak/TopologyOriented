#pragma once

#include <cstddef>
#include <iosfwd>


class IdString
{
public:
	IdString(char c, int idValue);

	IdString(const IdString& rhs);

	const char* c_str() const { return buff; }

private:
	static const std::size_t NUM_COLS = 5;
	static const std::size_t BUFFER_SIZE = NUM_COLS + 2;
	char buff[BUFFER_SIZE];
};

std::ostream& operator<<(std::ostream& os, const IdString& id);