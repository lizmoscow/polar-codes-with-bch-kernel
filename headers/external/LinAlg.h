#ifndef LINALG_H
#define LINALG_H

#define ALIGNMENT 32

#include "Codec.h"
#include <gsl/gsl_rng.h>
#include <iostream>
#include <istream>
#include <ostream>

//using namespace std;

//XOR two vectors.
void XOR(tBit* pDest, const tBit* pSrc, unsigned Size);
//XOR two vectors
void XOR(tBit* pDest, const tBit* pSrc1, const tBit* pSrc2, unsigned Size);
inline void XOR(unsigned* pDest, const unsigned* pSrc, unsigned Size)
{
	XOR((tBit*)pDest, (const tBit*)pSrc, sizeof(unsigned)*Size);
}
inline void XOR(unsigned* pDest, const unsigned* pSrc1, const unsigned* pSrc2, unsigned Size)
{
	XOR((tBit*)pDest, (const tBit*)pSrc1, (const tBit*)pSrc2, sizeof(unsigned)*Size);
}


//apply Gaussian elimination to the matrix
void Gauss(tBit* pMatrix //the matrix. 
		  ,unsigned& NumOfRows //number of rows in the matrix. On return contains its rank
		  ,unsigned NumOfColumns //number of columns in the matrix
		  ,unsigned NumOfAllocatedColumns //number of allocated columns in the matrix
		  ,bool ReversePass //true if reverse pass should be performed
		  ,unsigned* pPermutation //a permutation used to select the leading columns
		  );
//apply Gaussian elimination to a matrix over GF(2^m) (from right to left), permuting rows if necessary
void GaussRowPerm(
	unsigned m //field extension
	,int* pLogTable //logarithm table of the field
	,GFValue* pGF //field elements
	,GFValue* pMatrix //matrix to be Gaussed
	,unsigned& NumOfRows //input rows count; will equal row rank after calling 
	, unsigned NumOfColumns //columns count
	, unsigned NumOfAllocatedColumns  //number of allocated elements for each row
	, unsigned* pLeadingRows //first NumOfRows elements will be indices of non-zero rows (must be able to accomodate NumOfRows elements)
	, unsigned* pLeadingColumns = nullptr //index of the rightmost non-zero element in each row from pLeadingRows
	);

//apply Gaussian elimination to a matrix over GF(2) (from right to left), permuting rows if necessary
void GaussRowPerm(
	tBit* pMatrix
	, unsigned& NumOfRows
	, unsigned NumOfColumns
	, unsigned NumOfAllocatedColumns
	, bool ReversePass
	, unsigned* pLeadingRows //first NumOfRows elements will be indices of non-zero rows (must be able to accomodate NumOfRows elements)
	, unsigned* pLeadingColumns = nullptr //index of the rightmost 1 in each row from pLeadingRows
	);

//convert the matrix as obtained by Gauss into the minimum span form
void MinSpan(tBit* pMatrix //the matrix processed by Gauss
			,unsigned NumOfRows //number of rows in the matrix
			,unsigned NumOfColumns //number of columns in the matrix
			,unsigned NumOfAllocatedColumns //number of allocated columns in the matrix
			);

//allocate memory so that it is compatible with the linear algebra engine
tBit* AlignedAlloc(unsigned Size ///number of tBit entries in the array
				  );
//deallocate memory
void AlignedFree(tBit*& ptr);
//get the number of allocated columns for a matrix, so that each row is properly aligned
unsigned AlignColumnNumber(unsigned TrueNumOfColumns);

//construct a nullspace basis for a given matrix. 
tBit* Nullspace(tBit* pMatrix ///the matrix to be considered
	, unsigned Length ///number of columns in the matrix
	, unsigned NumOfAllocatedColumns ///number of bytes per row
	, unsigned& Dimension ///number of rows in the input matrix. On return, the dimension of the nullspace
	, unsigned* pPermutation = NULL //preferred column permutation
	);

//compute the number  of codewords of different weights for a code given by a packed generator matrix
void GetWeightSpectrum(unsigned Dimension,unsigned Length, const unsigned long long* pGenMatrix, unsigned long long* pSpectrum);
//make dynamic freezing matrix from a given check matrix
void MakeDynFrozen(unsigned m//log_2(code length)
	, const unsigned* pBasis //identifies a permutation of the check matrix
	, tBit* pFreezingMatrix // output: matrix with dynamic freezing constraints
	, unsigned* pPermutation // output: first NumOfChecks elements show positions of dynamic frozen symbols
	, const unsigned* pBitReversal ///bit reversal permutation
	, const tBit* pCheckMatrix ///input check matrix 
	, unsigned NumOfChecks ///number of check symbols
	);
unsigned* GetIndexA2ToindexA8(unsigned m);

void MakeNonSingMatrix(unsigned m, gsl_rng* pRNG, unsigned* pMatrix);

// apply Gaussian elimination only to the last row. 
// return true if the last row is linearly independent from other rows
bool GaussRowPermLastRow(
	tBit* pMatrix
	, unsigned& NumOfRows //number of rows (counting the one that has been added)
	, unsigned NumOfColumns
	, unsigned NumOfAllocatedColumns
	, unsigned* pLeadingColumns //pLeadingColumns[i] is the index of the rightmost 1 in the i-th row
	);

// print dynamic freezing constraints given matrix of them
void PrintFreezings(tBit* pMatrix,
	unsigned NumOfColumns,
	unsigned NumOfAllocatedColumns,
	unsigned NumOfRows,
	std::ostream &OutStream,
	unsigned* pColumnsOrder = nullptr
	);

///compute u(F^{\otimes m})^T
void MatrixKroneckerProductMultiplyTransposed(unsigned n, const tBit* pMatrix, unsigned KroneckerDegree, tBit* pVector);

//compute a dot product of two vectors
tBit DotProduct(unsigned N, const tBit* pX, const tBit* pY);

///compute a vectorized product of a vector and a matrix
///the elements of the vector are assumed to be blocks of size Stride
void MatrixMultiply(unsigned NumOfRows ///number of rows in a matrix
	, unsigned NumOfColumns ///number of columns in a matrix
	, unsigned Stride ///size of an input block in a matrix
	, const tBit* pInputVector///the input block vector 
	, const tBit* pMatrix ///the matrix
	, tBit* pOutputVector ///the output vector
	, bool ZeroOutput=true ///true if pOutputVector must be zeroed. Otherwise, the result of the operation will be XORed with it
);
///compute a vectorized product of a vector and a transposed matrix
///the elements of the vector are assumed to be blocks of size Stride
void MatrixMultiplyTransposed(unsigned NumOfRows ///number of rows in a matrix
	, unsigned NumOfColumns ///number of columns in a matrix
	, unsigned Stride ///size of an input block in a matrix
	, const tBit* pInputVector///the input block vector 
	, const tBit* pMatrix ///the matrix
	, tBit* pOutputVector ///the output vector
	, bool ZeroOutput=true///true if pOutputVector must be zeroed. Otherwise, the result of the operation will be XORed with it
);

// read dynamic freezing constraints from given input stream
void ReadFreezings(std::istream &SpecFile, 
	unsigned N,
	unsigned k,
	unsigned* pDecisionPoints,
	unsigned** ppFreezingConstraints,
	unsigned &LineNum
	);

void RemoveFrozen(tBit* pH, unsigned k, unsigned n, unsigned N=0);

template<typename T> void PrintMatrix(T* pA, unsigned rows, unsigned cols, std::ostream& out)
{
	for (unsigned i = 0; i < rows; i++)
	{
		for (unsigned j = 0; j < cols; j++)
		{
			out << int(pA[i*cols + j]) << ' ';
		}
		out << endl;
	}
};

// inverse a matrix
tBit* Inverse(tBit* pMatrix, unsigned NumOfRows);

#endif
