// Combination.h

/* combination
* ex) two elements from (0,1,2,3)
* (0,1) -> (0,2) -> (0,3) -> (1,2) -> ...
*/

#pragma once
#include <vector>
#include <fstream>

class Combination
{
private:
	std::vector<int> state; // current state
	int number; /* number of elements in the set (4 in the above example) */
	int sNum; /* number of elements in the combination (2 in the above example) */

public:
	/* constructor
	 * param src  number of elements in the set
	 * param src2 number of elements in the combination
	 */
	Combination(int src, int src2);
	/* go to next state
 	 * return  true  success to go
	 *         false fail to go (next state does not exist)
	 */
	bool Next(int src = 0);
	const int& operator[](int src) const; /* access */
	int size() const; /* number of element */
	friend std::ostream& operator<<(std::ostream &os, const Combination &src);
};
