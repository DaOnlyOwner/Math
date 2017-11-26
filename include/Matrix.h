#pragma once

#include <cstdint>
#include <vector>
#include "cassert"

#define MAKE_MATRIX_OPERATION_SCALAR( op )\
Matrix<VType> operator op (VType scalar)\
{\
	Matrix<VType> out;\
	for (int i = 0; i<m_numRow; i++)\
	{\
		for (int j = 0; j<m_numCol	; j++)\
		{\
			out(i,j) = (*this)(i,j) op scalar;\
		}\
	}\
	return out;\
}

#define MAKE_MATRIX_OPERATION( op )\
Matrix<VType> operator op (const Matrix<VType>& other)\
{\
	assert(other.ColumnCount() == ColumnCount() && RowCount() == other.RowCount());\
	Matrix<VType> out;\
	for (int i = 0; i<m_numRow; i++)\
	{\
		for (int j = 0; j<m_numCol	; j++)\
		{\
			out(i,j) = other(i,j) op (*this)(i,j);\
		}\
	}\
	return out;\
}\


#define MAKE_MATRIX_ASSIGNMENT_OPERATION(assign_op, op)\
Matrix<VType>& operator assign_op (const Matrix<VType>& other)\
{\
	*this = (*this op other);\
	return *this;\
}



namespace doo
{


		typedef int_fast8_t s8;
		typedef int_fast16_t s16;
		typedef int_fast32_t s32;

		typedef uint_fast8_t u8;
		typedef uint_fast16_t u16;
		typedef uint_fast32_t u32;
		
		typedef float f32;
		typedef double f64;
	



	namespace math
	{
		template <typename VType = float>
		class Matrix
		{

		public:
			Matrix() = default;
			Matrix(const Matrix<VType>&) = default;
			Matrix(Matrix<VType>&&) = default;
			
			Matrix(u32 numCol, u32 numRow, const VType* data) : m_numCol(numCol), m_numRow(numRow) 
			{
				m_data.resize(numCol * numRow);
				memcpy(m_data, data, m_numCol*m_numRow	 * sizeof(VType));
			}
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
				Matrix<VType> out(RowCount(), ColumnCount());
				for (int i = 0; i<RowCount(); i++)
				{
					for (int j = 0; j<other.ColumnCount(); j++)
					{
						VType result = 0;

						for (int k = 0; k<ColumnCount(); k++)
						{
							result += ro_at(k, j) * other(i,k);
						}
						out(i,j) = result;
					}
				}

				return out;
			}
			
			Matrix<VType> Transposed()
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

/*			void CalculateLUP(Matrix<VType>& L, Matrix<VType>& U, Matrix<VType>& P)
			{
				// Eliminations in each column 
				for(u32 n = VSize-1; n>0; n--)
				{
					// Start where the first zero should be at
					for(u32 i = (VSize - n); i<VSize; i++)
					{
						// Calculate k, the multiplicator and store it inside L(


					}
				}
			}*/

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

		private:

			VType ro_at(u32 i, u32 j) const
			{
				return m_data[offset(i,j)];
			}

			u32 offset(u32 i, u32 j) const
			{
				assert(i < m_numCol && j< m_numRow);
				return i * m_numRow + j;
			}

			std::vector<VType> m_data;
			u32 m_numCol = 0;
			u32 m_numRow = 0;
		};
		
	}
}