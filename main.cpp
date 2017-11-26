#include "include/Matrix.h"
#include "include/Vector.h"
#include <iostream>

#define SIZE 15
#define T double


void printMat(doo::math::Matrix<SIZE,T>& m)
{
	for(int i = 0; i<SIZE; i++)
	{
		for(int j =0 ; j<SIZE; j++)
		{
			std::cout << m[doo::math::Matrix<SIZE,T>::Get(i, j)] << "\t";
		}

		std::cout << std::endl;
	};
	
}

int main()
{

	double in[SIZE*SIZE];

	unsigned int count = 0;
	for(int i = 0; i<SIZE; i++)
	{
		for(int j = 0; j<SIZE; j++)
		{
			in[i* SIZE + j] = count;
			count++;
		}
	}

	doo::math::Matrix<SIZE,T> m(in);

	m = m*m;

	printMat(m);

}
