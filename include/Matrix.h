#pragma once
#include "cstring"
#include "cassert"
#include <cstdint>

#define MAKE_MATRIX_OPERATION_SCALAR( op )\
Matrix<VSize,VType> operator op (VType p_scalar)\
{\
	Matrix<VSize,VType> out;\
	for (int i = 0; i<VSize; i++)\
	{\
		for (int j = 0; j<VSize; j++)\
		{\
			out[i*VSize + j] = m_data[i*VSize + j] op p_scalar;\
		}\
	}\
	return out;\
}

#define MAKE_MATRIX_OPERATION( op )\
Matrix<VSize,VType> operator op (const Matrix<VSize,VType>& p_other)\
{\
	Matrix<VSize,VType> out;\
	for (int i = 0; i<VSize; i++)\
	{\
		for (int j = 0; j<VSize; j++)\
		{\
			out[i*VSize + j] = p_other.m_data[i*VSize + j] op m_data[i*VSize + j];\
		}\
	}\
	return out;\
}\


#define MAKE_MATRIX_ASSIGNMENT_OPERATION(assign_op, op)\
Matrix<VSize,VType>& operator assign_op (const Matrix<VSize,VType>& p_other)\
{\
	m_data = (*this op p_other).m_data;\
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
			Matrix(const VType data[VSize]) 
			{
				memcpy(m_data, data, VSize*VSize * sizeof(VType));
			}

			MAKE_MATRIX_OPERATION_SCALAR(+);
			MAKE_MATRIX_OPERATION_SCALAR(*);
			MAKE_MATRIX_OPERATION_SCALAR(/ );
			MAKE_MATRIX_OPERATION_SCALAR(-);

			MAKE_MATRIX_OPERATION(+);
			MAKE_MATRIX_OPERATION(-);

			MAKE_MATRIX_ASSIGNMENT_OPERATION(+=, +);
			MAKE_MATRIX_ASSIGNMENT_OPERATION(*=, *);

			void Row(uint p_i, const Vec<VSize,VType>& p_vec)
			{
				for(int j = 0; j<VSize; j++)
				{
					m_data[offset(p_i, j)] = p_vec[j];
				}
			}
		
			const VType* RowView(uint p_i) const
			{
				return m_data + offset(p_i,0);
			}
			
			void Column(uint p_j, const Vec<VSize,VType>& p_vec)
			{
				for (int i = 0; i<VSize; i++)
				{
					m_data[offset(i, p_j)] = p_vec[i];
				}
			}

			const VType* ColumnView(uint p_j) const
			{
				return m_data + offset(0, p_j);
			}

			VType& operator() (uint p_i, uint p_j)
			{
				return m_data[offset(p_i,p_j)];
			}

			Matrix<VSize,VType> operator* (const Matrix<VSize,VType>& p_other)
			{
				Matrix<VSize,VType> out;
				for (int i = 0; i<VSize; i++)
				{
					for (int j = 0; j<VSize; j++)
					{
						VType result = 0;

						for (int k = 0; k<VSize; k++)
						{
							result += m_data[k*VSize + j] * p_other.m_data[i * VSize + k];
						}
						out.m_data[i * VSize + j] = result;
					}
				}

				return out;
			}
			
			Matrix<VSize,VType> Transposed()
			{
				Matrix<VSize,VType> out;
				for (int i = 0; i<VSize; i++)
				{
					for (int j = 0; j<VSize; j++)
					{
						out[j*VSize + i] = m_data[i*VSize + j];
					}
				}

				return out;
			}
			Matrix<VSize,VType>& Transpose()
			{
				*this = Transpose();
				return *this;
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

			void CalculateLUP(Matrix<VSize,VType>& p_L, Matrix<VSize,VType>& p_U, Matrix<VSize,VType>& p_P)
			{
				// Eliminations in each column 
				for(uint n = VSize-1; n>0; n--)
				{
					// Start where the first zero should be at
					for(uint i = (VSize - n); i<VSize; i++)
					{
						// Calculate k, the multiplicator and store it inside L(


					}
				}
			}


		private:
			static uint offset(uint i, uint j)
			{
				assert(i < VSize && j< VSize);
				return i * VSize + j;
			}

			VType m_data[VSize * VSize];
		};
		
	}
}