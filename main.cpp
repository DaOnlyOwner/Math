#include "include/Matrix.h"
#include <cstdio>
#include "Windows.h"

#define SIZEI 10
#define SIZEJ 10
#define T double


void printMat(const doo::math::Matrix<T>& m)
{
	for(int i = 0; i<m.RowCount(); i++)
	{
		for(int j =0 ; j<m.ColumnCount(); j++)
		{
			printf("%.2f\t", m(i,j));
		}

		printf("\n");
	};

	printf("\n");
	
}

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

	doo::math::Matrix<T> m(SIZEI,SIZEJ);
	fillMat(m);

	doo::math::Matrix<T> l, u, p;

	m.DooLittleDecomposed(l, u, p, 0.0001);

	printMat(l);
	printMat(u);


	return 0;

}
