#include "include/Matrix.h"
#include <iostream>

#define SIZEI 6
#define SIZEJ 6
#define T double

void fillMat(doo::math::Matrix<T>& m)
{
	unsigned int count = 1;
	for(int i = 0; i<m.RowCount(); i++)
	{
		for(int j = 0; j<m.ColumnCount(); j++)
		{
			m(i,j) = count;
			//std::cout << count;
			count++;
		}
	}
}

int main()
{

	doo::math::Matrix<T> a(SIZEI,SIZEJ);
	fillMat(a);
    a.Print();

	doo::math::Matrix<T> l, u, p;

	a.LUPDecomposed(l, u, p);

	(p * l * u).Print();

	return 0;

}
