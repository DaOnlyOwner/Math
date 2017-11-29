#include "include/Matrix.h"
#include "Windows.h"

#define SIZEI 7
#define SIZEJ 7
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

	//l.Print();
	//u.Print();
	//p.Print();

    (p * l * u).Print();

	//printf("Great finaleeeeee \n\n\n");
    //(p * l * u).Print();


    //getchar();

    /*const T in[] = {0,1,0,0,
                    0,0,0,1,
                    0,0,1,0,
                    1,0,0,0};

    const T in2[] = {1,2,3,4,
                     5,6,7,8,
                     9,10,11,12,
                    13,14,15,16};


    doo::math::Matrix<T> permut(4,4,in);
    doo::math::Matrix<T> test(4,4, in2);

    (permut * test).Print();*/


	return 0;

}
