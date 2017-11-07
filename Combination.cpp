// Combination.cpp

#include "Combination.h"
#include <cassert>

Combination::Combination(int src, int src2)
{
	assert(src >= src2);
	number = src;
	sNum = src2;
	// initialize
	state.resize(src2);
	for (int i(0); i < src2; i++) {
		state[i] = i;
	}
}

bool Combination::Next(int src)
{
	int i, check;
	check = sNum - 1 - src;
	for (;;) {
		state[check]++;
		if (state[check] < number + check - sNum + 1) {
			break;
		}
		check--;
		if (check < 0) {
			return false; // not exist!
		}
	}
	for (i = check + 1; i < sNum; i++) {
		state[i] = state[check] + i - check;
	}
	return true;
}

const int& Combination::operator[](int src) const
{
	assert(0 <= src && src<sNum);
	return state[src];
}

int Combination::size() const
{
	return sNum;
}

std::ostream& operator<<(std::ostream &dest, const Combination &src)
{
	int i;
	dest << src.state[0];
	for (i = 1; i<src.sNum; i++) dest << " " << src.state[i];
	return dest;
}
