#pragma once

#include <vector>
#include "cstdio"
#include "cassert"
#include "Common.h"

#define MAKE_MATRIX_OPERATION_SCALAR( op )\
Matrix<VType> operator op (VType scalar)\
{\
	Matrix<VType> out;\
	for (int i = 0; i<m_numRow*m_numCol; i++)\
	{\
		out.m_data[i] = m_data[i] op scalar;\
	}\
	return out;\
}

#define MAKE_MATRIX_OPERATION( op )\
Matrix<VType> operator op (const Matrix<VType>& other)\
{\
	assert(other.ColumnCount() == ColumnCount() && RowCount() == other.RowCount());\
	Matrix<VType> out;\
	for (int i = 0; i<m_numRow*m_numCol; i++)\
	{\
		out[i] = m_data[i] op other.m_data[i];\
	}\
	return out;\
}


#define MAKE_MATRIX_ASSIGNMENT_OPERATION(assign_op, op)\
Matrix<VType>& operator assign_op (const Matrix<VType>& other)\
{\
	*this = (*this op other);\
	return *this;\
}



namespace doo
{
	namespace math
	{
		template <typename VType = float>
		class Matrix
		{

		public:
			Matrix() = default;
			Matrix(const Matrix<VType>&) = default;
			Matrix(Matrix<VType>&&) = default;
			
			Matrix(u32 numCol, u32 numRow, const VType* data) : m_numCol(numCol), m_numRow(numRow), m_data(data, data + numCol * numRow)
			{}
			
			Matrix(u32 numCol, u32 numRow) : m_numCol(numCol), m_numRow(numRow) 
			{
				m_data.resize(numCol * numRow);
			}


			Matrix<VType>& operator=(const Matrix<VType>&) = default;
			Matrix<VType>& operator=(Matrix<VType>&&) = default;
			MAKE_MATRIX_OPERATION_SCALAR(+);
			MAKE_MATRIX_OPERATION_SCALAR(*);
			MAKE_MATRIX_OPERATION_SCALAR(/ );
			MAKE_MATRIX_OPERATION_SCALAR(-);

			MAKE_MATRIX_OPERATION(+);
			MAKE_MATRIX_OPERATION(-);

			MAKE_MATRIX_ASSIGNMENT_OPERATION(+=, +);
			MAKE_MATRIX_ASSIGNMENT_OPERATION(*=, *);

			VType& operator() (u32 i, u32 j)
			{
				return m_data[offset(i,j)];
			}

			
			const VType& operator() (u32 i, u32 j) const
			{
				return m_data[offset(i,j)];
			}

			u32 RowCount() const
			{
				return m_numRow;
			}

			u32 ColumnCount() const
			{
				return m_numCol;
			}

			Matrix<VType> operator* (const Matrix<VType>& other)
			{
				assert(ColumnCount() == other.RowCount());
				Matrix<VType> out(RowCount(), other.ColumnCount());
                #pragma omp parallel for collapse(2)
				for (int i = 0; i<RowCount(); i++)
				{
					for (int j = 0; j<other.ColumnCount(); j++)
					{
						VType result = 0;

						for (int k = 0; k<ColumnCount(); k++)
						{
							result += ro_at(i, k) * other(k,j);
						}
						out(i,j) = result;
					}
				}

				return out;
			}
			
			Matrix<VType> Transposed() const
			{
				Matrix<VType> out;
				for (int i = 0; i<m_numRow; i++)
				{
					for (int j = 0; j<m_numCol; j++)
					{
						out(j,i) = ro_at(i, j);
					}
				}

				return out;
			}

			/*VType Determinant()
			{
				if constexpr(VSize == 2) return a1() * b2() - b1() * a2();
				if constexpr(VSize == 3)
				{
					VType a = a1() * (b2() * c3() - c2()*b3());
					VType b = b1() * (a2() * c3() - c2() * a3());
					VType c = c1() * (a2() * b3() - b2() * a3());
					return a - b + c;
				}

				if constexpr(VSize > 3)
				{
					return calcDet(m_data, VSize);
				}


			}*/

            void Resize(u32 numRow, u32 numCol)
            {
                m_data.resize(numRow * numCol);
                m_numRow = numRow;
                m_numCol = numCol;
            }

			void SwapRow(u32 i, u32 i2)
			{
				for(u32 colIndex = 0; colIndex < ColumnCount(); colIndex++)
				{
					VType temp = ro_at(i, colIndex);
					(*this)(i, colIndex) = ro_at(i2, colIndex);
					(*this)(i2, colIndex) = temp;
				}
			}

			u32 LUPDecomposed(Matrix<VType>& L, Matrix<VType>& U, Matrix<VType>& P) const
			{
				std::vector<u32> pivotVec;
				Matrix<VType> LU;
				u32 out = LUPDecomposed(LU, pivotVec);
				L.Resize(LU.RowCount(),LU.ColumnCount());
				L.SetToIdentity();
				U.Resize(LU.RowCount(),LU.ColumnCount());
				P.Resize(LU.RowCount(), LU.ColumnCount());
				for(u32 i = 0; i<RowCount(); i++)
				{
					for(u32 j = 0; j<ColumnCount(); j++)
					{
						if(j>=i)
						{
							U(i, j) = LU(i, j);
						}

						else L(i, j) = LU(i, j);

					}
				}

				for(u32 i = 0; i<pivotVec.size(); i++)
				{
					P(pivotVec[i],i) = 1;
				}

				return out;

			}


			u32 LUPDecomposed(Matrix<VType>& LU, std::vector<u32>& pivotVec) const
			{
				assert(RowCount() == ColumnCount());
                pivotVec.reserve(ColumnCount());
                for(u32 i = 0; i<ColumnCount(); i++) pivotVec.push_back(i);

				LU = (*this);
				u32 numPivots = 0;

				for(u32 j = 0; j<ColumnCount(); j++)
				{
					pivot(LU, pivotVec, numPivots, j);

                    for(u32 i = j+1; i<RowCount(); i++)
					{
						// Compute the multiplicator k:
						VType k = LU(i, j) / LU(j, j);
                        //printf("%f ... %f\n",LU(i,j), LU(j,j));
                        if ( LU(i,j) == 0 && LU(j,j) == 0) k = 1; // Hack to ensure stability
                        LU(i, j) = k; // Fill  with the multiplicators
						for(u32 j_ = j+1; j_ < ColumnCount(); j_++)
						{
							LU(i, j_) = LU(i,j_) - k * LU(j,j_);
						}
					}
				}

				return numPivots;

			}

			Matrix<VType>& SetToIdentity()
			{
				assert(m_numCol == m_numRow);
				for(u32 i = 0; i<m_numRow; i++)
				{
					for(u32 j= 0; j<m_numCol; j++)
					{
						if(j == i)
							(*this)(i,i) = 1;
						else (*this)(i,j) = 0;
					}
				}
				return *this;
			}

			Matrix<VType>& SetToZero()
			{
				m_data.assign(m_data.size(), 0);
				return *this;
			}

            void Print() const
            {
                for(int i = 0; i<RowCount(); i++)
                {
                    for(int j =0 ; j<ColumnCount(); j++)
                    {
                        printf("%.2f\t",ro_at(i,j));
                    }
					printf("\n");
                };
				
				printf("\n");
			}


		private:
			void pivot(Matrix<VType>& LU, std::vector<u32>& pivotVec, u32& numPivots, u32 j) const
			{
				// Get the max element in the jth column starting from j+1 (watch the diagonal):
				VType max_elem = abs(LU(j, j));
				u32 max_row = j;
				for(u32 i = j + 1; i < RowCount(); i++)
				{
					VType elem_to_inspect = abs(LU(i, j));
					if (elem_to_inspect > max_elem) {
						max_elem = elem_to_inspect; max_row = i;
					}
				}
				// Okay now we found the max element and the corresponding row -> Switch it with the j th row.
				if(max_row != j) // The first row had the max elem -> Don't need to swap rows
				{
					LU.SwapRow(j, max_row);
					// j and max_row swap now --> The permutation of j and the permutation of max_row have to swap accordingly.
					u32 permutRow = pivotVec[j];
					pivotVec[j] = pivotVec[max_row];
					pivotVec[max_row] = permutRow;

					numPivots++;
				}
			}

			VType ro_at(u32 i, u32 j) const
			{
				return m_data[offset(i,j)];
			}

			u32 offset(u32 i, u32 j) const
			{
                if (i >= m_numRow || j>=m_numCol)
                {
                    __debugbreak();
                }
				assert(i < m_numRow && j< m_numCol);
				return i * m_numCol + j;
			}

			u32 m_numCol = 0;
			u32 m_numRow = 0;
			std::vector<VType> m_data;
		};
		
	}
}