#include <emmintrin.h>
#include <xmmintrin.h>
#include <algorithm>
#include <ostream>
#include "LinAlg.h"
#include "ArikanEncoder.h"

using namespace std;

#ifdef __AVX__ 
typedef __m256d MMType;
#define MMXOR(Dest,Src1,Src2) *(MMType*)(Dest) = _mm256_xor_pd(*(const MMType*)(Src1), *(const MMType*)(Src2));
#else 
typedef unsigned  MMType;
#define MMXOR(Dest,Src1,Src2) *(MMType*)(Dest) = *(const MMType*)(Src1)^*(const MMType*)(Src2);
#endif

/// allocate memory so that it is compatible with the linear algebra engine
/// @param number of tBit entries in the array
tBit* AlignedAlloc(unsigned Size)
{
    return (tBit*)_aligned_malloc(Size*sizeof(tBit), ALIGNMENT);
}

void AlignedFree(tBit*& ptr)
{
    _aligned_free(ptr);
    ptr = nullptr;
}

void XOR(tBit* pDest, const tBit* pSrc, unsigned Size)
{
	unsigned i = 0;
	unsigned SIMDSize = (Size / sizeof(MMType)) * sizeof(MMType);
	//assert(_CrtCheckMemory());
	for (i = 0; i < SIMDSize; i += sizeof(MMType))
		MMXOR(pDest + i, pDest + i, pSrc + i);
	for (; i < (Size&~3); i += 4)
	{
		*((unsigned*)(pDest + i)) ^= *((unsigned*)(pSrc + i));
	};
	for (; i < Size; i++)
		pDest[i] ^= pSrc[i];
};
void XOR(tBit* pDest, const tBit* pSrc1, const tBit* pSrc2, unsigned Size)
{
	unsigned i = 0;
	unsigned SIMDSize = (Size / sizeof(MMType)) * sizeof(MMType);
	for (i = 0; i < SIMDSize; i += sizeof(MMType))
		MMXOR(pDest + i, pSrc1 + i, pSrc2 + i);
	for (; i < Size; i++)
		pDest[i] = pSrc1[i] ^ pSrc2[i];
/*#ifndef NDEBUG
    for (uint32_t i = 0; i < Size; ++i) {
        assert(pSrc1[i] == BIT_0 || pSrc1[i] == BIT_1);
        assert(pSrc2[i] == BIT_0 || pSrc2[i] == BIT_1);
    }
#endif
*/
}

//get the number of allocated columns for a matrix, so that each row is properly aligned
unsigned AlignColumnNumber(unsigned TrueNumOfColumns)
{
    if (TrueNumOfColumns%ALIGNMENT)
        return TrueNumOfColumns + ALIGNMENT - TrueNumOfColumns%ALIGNMENT;
    else
        return TrueNumOfColumns;
}

//apply Gaussian elimination to a matrix, permuting columns if necessary
void Gauss(tBit* pMatrix, unsigned& NumOfRows, unsigned NumOfColumns, unsigned NumOfAllocatedColumns, bool ReversePass, unsigned* pPermutation)
{
    for (unsigned i = 0; i < NumOfRows; i++) {
        //identify the leading column
        unsigned c = i;
        bool Success = false;
        for (c; c < NumOfColumns; c++) {
            unsigned C = (pPermutation) ? pPermutation[c] : c;
            for (unsigned j = i; j<NumOfRows; j++) {
                if (pMatrix[j*NumOfAllocatedColumns + C]) {
                    Success = true;
                    if (j>i)
                        XOR(pMatrix + i*NumOfAllocatedColumns, pMatrix + j*NumOfAllocatedColumns, NumOfColumns);
                    break;
                }
            }
            if (Success) {
                if ((c != i) && pPermutation)
                    swap(pPermutation[c], pPermutation[i]);
                break;
            }
        }
        if (!Success) {
            NumOfRows = i;
            break;
        }
        unsigned LoopStart = (ReversePass) ? 0 : (i + 1);
        unsigned C = (pPermutation) ? pPermutation[i] : c;
        for (unsigned j = LoopStart; j < NumOfRows; j++) {
            if (j == i) continue;
            if (pMatrix[j*NumOfAllocatedColumns + C])
                XOR(pMatrix + j*NumOfAllocatedColumns, pMatrix + i*NumOfAllocatedColumns, NumOfColumns);
        }
    }
}

//apply Gaussian elimination to a matrix over GF(2^m) (from right to left), permuting rows if necessary
void GaussRowPerm(
    unsigned m //field extension
    , int* pLogTable //logarithm table of the field
    , GFValue* pGF //field elements
    , GFValue* pMatrix //matrix to be Gaussed
    , unsigned& NumOfRows //input rows count; will equal row rank after calling 
    , unsigned NumOfColumns //columns count
    , unsigned NumOfAllocatedColumns  //number of allocated elements for each row
    , unsigned* pLeadingRows //first NumOfRows elements will be indices of non-zero rows (must be able to accomodate NumOfRows elements)
    , unsigned* pLeadingColumns //index of the rightmost non-zero element in each row from pLeadingRows
    )
{
    unsigned r = 0;
    for (r = 0; r < NumOfRows; r++) {
        pLeadingRows[r] = r;
    }
    int c = NumOfColumns, col;
    GFValue *pLeadRow, *pCurRow;
    int activeRowsCount = 0;
    unsigned n = 1 << m;
    //assert(_CrtCheckMemory());
    while (--c >= 0) {
        //find row with non-zero element in column c
        for (r = activeRowsCount; r < NumOfRows; r++) {
            if (pMatrix[pLeadingRows[r] * NumOfAllocatedColumns + c]) break;
        }
        //if no such row, go to next column index
        if (r == NumOfRows) {
            continue;
        }
        //mark current leading column
        if (pLeadingColumns) {
            pLeadingColumns[activeRowsCount] = c;
        }
        //mark current leading row as active
        pLeadRow = pMatrix + (pLeadingRows[r] * NumOfAllocatedColumns);
        swap(pLeadingRows[r], pLeadingRows[activeRowsCount]);
        activeRowsCount++;
        unsigned L = pLogTable[pLeadRow[c]];
        pLeadRow[c] = 1;
        //normalize leading row
        for (col = c - 1; col >= 0; col--) {
            if (pLeadRow[col]) {
                int L1 = pLogTable[pLeadRow[col]] - L;
                L1 = (L1 < 0) ? (L1 + n - 1) : L1;
                pLeadRow[col] = pGF[1 + L1];
            }
        }
        //null column c in other rows
        for (r++; r < NumOfRows; r++) {
            pCurRow = pMatrix + (pLeadingRows[r] * NumOfAllocatedColumns);
            if (pCurRow[c]) {
                int L1 = pLogTable[pCurRow[c]];
                pCurRow[c] = 0;
                //normalize current row
                for (col = c - 1; col >= 0; col--) {
                    if (pCurRow[col]) {
                        int L2 = pLogTable[pCurRow[col]] - L1;
                        L2 = (L2 < 0) ? (L2 + n - 1) : L2;
                        pCurRow[col] = pGF[1 + L2];
                    }
                }
                //subtract leading row from current row
                for (col = c - 1; col >= 0; col--) {
                    pCurRow[col] ^= pLeadRow[col];
                }
            }
        }
    }
    NumOfRows = activeRowsCount;
}

//apply Gaussian elimination to a matrix over GF(2) (from right to left), permuting rows if necessary
void GaussRowPerm(
    tBit* pMatrix
    , unsigned& NumOfRows
    , unsigned NumOfColumns
    , unsigned NumOfAllocatedColumns
    , bool ReversePass
    , unsigned* pLeadingRows //first NumOfRows elements will be indices of non-zero rows (must be able to accomodate NumOfRows elements)
    , unsigned* pLeadingColumns //index of the rightmost 1 in each row from pLeadingRows
    )
{
    unsigned r;
    int c = NumOfColumns;
    tBit *pLeadRow, *pCurRow;
    int activeRowsCount = 0;
    while (--c >= 0) {
        //find row with non-zero element in column c
        for (r = activeRowsCount; r < NumOfRows; r++) {
            if (pMatrix[pLeadingRows[r] * NumOfAllocatedColumns + c]) break;
        }
        //if no such row, go to next column index
        if (r == NumOfRows) {
            continue;
        }
        //mark current leading column
        if (pLeadingColumns) {
            pLeadingColumns[activeRowsCount] = c;
        }
        //mark current leading row as active
        swap(pLeadingRows[r], pLeadingRows[activeRowsCount]);
        pLeadRow = pMatrix + (pLeadingRows[activeRowsCount] * NumOfAllocatedColumns);
        //null column c in other rows
        if (ReversePass) {
            for (r = 0; r < NumOfRows; r++) {
                if (r == activeRowsCount) continue;
                pCurRow = pMatrix + (pLeadingRows[r] * NumOfAllocatedColumns);
                if (pCurRow[c]) {
                    XOR(pCurRow, pLeadRow, c + 1);
                }
            }
        } else {
            for (r++; r < NumOfRows; r++) {
                pCurRow = pMatrix + (pLeadingRows[r] * NumOfAllocatedColumns);
                if (pCurRow[c]) {
                    XOR(pCurRow, pLeadRow, c + 1);
                }
            }
        }
        activeRowsCount++;
    }
    NumOfRows = activeRowsCount;
}

//convert the check matrix into minimum span form
void MinSpan(tBit* pMatrix, unsigned NumOfRows, unsigned NumOfColumns, unsigned NumOfAllocatedColumns)
{
	for (int i = NumOfRows - 1; i > 0; i--)
	{
		//identify the leftmost 1
		unsigned c = 0;
		while (!pMatrix[i*NumOfAllocatedColumns + c]) c++;
		for (unsigned j = 0; j < i; j++)
			if (pMatrix[j*NumOfAllocatedColumns + c])
				XOR(pMatrix + j*NumOfAllocatedColumns, pMatrix + i*NumOfAllocatedColumns, NumOfColumns);
	}
};


//construct a nullspace basis for a given matrix. 
tBit* Nullspace(tBit* pMatrix ///the matrix to be considered
    , unsigned Length ///number of columns in the matrix
    , unsigned NumOfAllocatedColumns ///number of tBits per row
    , unsigned& Dimension ///number of rows in the input matrix. On return, the dimension of the nullspace
    , unsigned* pPermutation  //preferred column permutation
    )
{
    unsigned* pPerm = nullptr;
    if (pPermutation)
        pPerm = pPermutation;
    else {
        pPerm = new unsigned[Length];
        for (unsigned i = 0; i < Length; i++)
            pPerm[i] = i;
    }
    Gauss(pMatrix, Dimension, Length, NumOfAllocatedColumns, true, pPerm);
    tBit* pNullspace = AlignedAlloc(NumOfAllocatedColumns*(Length - Dimension));
    memset(pNullspace, 0, sizeof(tBit)*NumOfAllocatedColumns*((Length - Dimension)));
    for (unsigned i = 0; i < Length - Dimension; i++) {
        pNullspace[i*NumOfAllocatedColumns + pPerm[Dimension + i]] = 1;
        for (unsigned j = 0; j < Dimension; j++) {
            pNullspace[i*NumOfAllocatedColumns + pPerm[j]] = pMatrix[j*NumOfAllocatedColumns + pPerm[Dimension + i]];
        }
    }
    Dimension = Length - Dimension;
    if (pPermutation == nullptr)
        delete[] pPerm;
    return pNullspace;
}

#ifdef __AVX__
//this should actually be SSE3
#if !defined(_M_X64)
#define Weight(c) _mm_popcnt_u32(c&0xFFFFFFFFul)+_mm_popcnt_u32(c>>32)
#else
#define Weight(c) _mm_popcnt_u64(c)
#endif
#else
//compute weight by bit masking
unsigned Weight(unsigned long long c)
{
    //split even and odd bits
    unsigned long long c1 = c & 0x5555555555555555ull;
    unsigned long long c2 = (c >> 1) & 0x5555555555555555ull;
    c1 += c2;//now we have weight of bit tuples
    c2 = (c1 >> 2) & 0x3333333333333333ull;
    c1 &= 0x3333333333333333ull;
    c1 += c2;//now we have weight of bit quadruples
    c2 = (c1 >> 4) & 0x0F0F0F0F0F0F0F0Full;
    c1 &= 0x0F0F0F0F0F0F0F0Full;
    c1 += c2;//now we have weight of bytes
    c1 += c1 >> 32;//sum byte weights from low and upper dwords
    c1 += c1 >> 16;
    c1 += c1 >> 8;
    return c1 & 0xFF;
}
#endif

//compute weight spectrum of a code given by a packed generator matrix
void GetWeightSpectrumPrim(unsigned Dimension, const unsigned long long * pGenMatrix, unsigned long long* pSpectrum)
{
    unsigned long long c = 0;
    pSpectrum[0] = 1;
    for (unsigned long long i = 1; i < 1ull << Dimension; i++) {
        unsigned q = 0;
        while (!(i&(1ull << q)))
            q++;
        c ^= pGenMatrix[q];
        unsigned w = (unsigned)Weight(c);
        pSpectrum[w]++;
    }
}

/*Evaluate Kravchuk polynomial P_k(x,n)*/
long long Kravchuk(unsigned k, unsigned x, unsigned n)
{
    long long R = 0;
    unsigned j;
    for (j = 0; j <= k; j++) {
        unsigned __int64 B1, B2, P;
        B1 = Binomial(x, j);
        B2 = Binomial(n - x, k - j);
        P = B1*B2;
        if (j & 1)
            R -= P;
        else
            R += P;
    }
    return R;
}


//compute weight spectrum of a code given by a packed generator matrix
void GetWeightSpectrum(unsigned Dimension, unsigned Length, const unsigned long long * pGenMatrix, unsigned long long* pSpectrum)
{

    if (Dimension < Length / 2) {
        GetWeightSpectrumPrim(Dimension, pGenMatrix, pSpectrum);
        return;
    }
    //for hogh-rate codes it is easier to compute the spectrum of a dual code
    unsigned N = AlignColumnNumber(Length);
    tBit* pGeneratorMatrix = AlignedAlloc(Dimension*N);
    //unpack the matrix
    memset(pGeneratorMatrix, 0, sizeof(tBit)*Dimension*N);
    for (unsigned i = 0; i < Dimension; i++) {
        unsigned long long R = pGenMatrix[i];
        tBit* pDest = pGeneratorMatrix + i*N;
        for (unsigned j = 0; j < Length; j++) {
            pDest[j] = (R >> j) & 1;
        }
    }
    unsigned k = Dimension;
    tBit* pNullspace = Nullspace(pGeneratorMatrix, Length, N, k, NULL);
    AlignedFree(pGeneratorMatrix);
    //pack the check matrix
    unsigned long long* pCheckMatrix = new unsigned long long[k];
    for (unsigned i = 0; i < k; i++) {
        unsigned long long R = 0;
        const tBit* pSrc = pNullspace + i*N;
        for (unsigned j = 0; j < Length; j++)
        if (pSrc[j])
            R |= 1ull << j;
        pCheckMatrix[i] = R;
    }
    AlignedFree(pNullspace);
    unsigned long long *pDualSpectrum = new unsigned long long[Length + 1];
    memset(pDualSpectrum, 0, sizeof(pDualSpectrum[0])*(Length + 1));
    GetWeightSpectrumPrim(k, pCheckMatrix, pDualSpectrum);
    delete[] pCheckMatrix;
    //apply McWilliams identities
    for (unsigned j = 0; j <= Length; j++) {
        __int64 A = 0;
        for (unsigned i = 0; i <= Length; i++)
            A += pDualSpectrum[i] * Kravchuk(j, i, Length);
        pSpectrum[j] = (unsigned)(A >> (Length - Dimension));
    }
    delete[] pDualSpectrum;
}

//permute basis of given check matrix
void MakeDynFrozen(unsigned m//log_2(code length)
    , const unsigned* pBasis //identifies a permutation of the check matrix
    , tBit* pFreezingMatrix // output: matrix with dynamic freezing constraints
    , unsigned* pPermutation // output: first NumOfChecks elements show positions of dynamic frozen symbols
    , const unsigned* pBitReversal ///bit reversal permutation
    , const tBit* pCheckMatrix ///input check matrix 
    , unsigned NumOfChecks ///number of check symbols
    )
{
    unsigned LengthAligned = AlignColumnNumber(1u << m);
    for (unsigned i = 0; i < 1u << m; i++) {
        unsigned x = 0;
        for (unsigned j = 0; j < m; j++)
        if ((i >> j) & 1)
            x ^= pBasis[j];
        for (unsigned j = 0; j < NumOfChecks; j++)
            pFreezingMatrix[(j << m) + pBitReversal[i]] = pCheckMatrix[j*LengthAligned/*(j << m)*/ + x];
    }
    for (unsigned j = 0; j < NumOfChecks; j++)
        ArikanTransposed(m, pFreezingMatrix + (j << m));
    //Arikan(pPermutedMatrix+(j<<m),1u<<m);
    //identify frozen channels
    for (unsigned i = 0; i < 1u << m; i++)
        pPermutation[i] = (1u << m) - 1 - i;
    unsigned NumOfChecks2 = NumOfChecks;
    Gauss(pFreezingMatrix, NumOfChecks2, 1u << m, 1u << m, true, pPermutation);
}

unsigned* GetIndexA2ToindexA8(unsigned m)
{
    unsigned* pPermutation = new unsigned[1ull << m];
    for (unsigned i = 0; i < (1u << m); i++) {
        // Index=x,x,x,x,x, where x=0..7, the last x corresponds to least significant bits
        pPermutation[i] = i;
        for (int w = m - 3; w >= 0; w -= 3) {
            unsigned x = (i >> w) & 7;
            if ((x == 3) || (x == 4)) pPermutation[i] ^= 7 << w;
        }
    }
    return pPermutation;
}

void MakeNonSingMatrix(unsigned m, gsl_rng* pRNG, unsigned* pMatrix)
{
    unsigned* pNormalizedMatrix = (unsigned*)alloca(sizeof(unsigned)*m);
    unsigned* pPermutation = (unsigned*)alloca(sizeof(unsigned)*m);
    pNormalizedMatrix[0] = pMatrix[0] = 1;
    pPermutation[0] = 0;
    for (unsigned i = 0; i < m; i++) {
        unsigned Y, X;
        do {
            X = Y = 1 + gsl_rng_uniform_int(pRNG, (1u << m) - 1);
            //check if it is linearly dependent
            for (unsigned j = 0; j < i; j++)
            if ((X >> pPermutation[j]) & 1)
                X ^= pNormalizedMatrix[j];
        } while (!X);
        pNormalizedMatrix[i] = pMatrix[i] = X;
        for (unsigned j = 0; j < m; j++)
        if ((X >> j) & 1) {
            pPermutation[i] = j;
            for (unsigned s = 0; s < i; s++)
            if (pNormalizedMatrix[s] & (1u << j))
                pNormalizedMatrix[s] ^= X;
            break;
        }
    }
}

// apply Gaussian elimination only to the last row. 
// return true if the last row is linearly independent from other rows
bool GaussRowPermLastRow(
	tBit* pMatrix
	, unsigned& NumOfRows //number of rows (counting the one that has been added)
	, unsigned NumOfColumns
	, unsigned NumOfAllocatedColumns
	, unsigned* pLeadingColumns //pLeadingColumns[i] is the index of the rightmost 1 in the i-th row
);

/*void MakeDynFrozenSortedArikan8(unsigned m//log_2(code length)
	, const unsigned* pBasis //identifies a permutation of the check matrix
	, tBit* pFreezingMatrix // output: matrix with dynamic freezing constraints
	, unsigned* pPermutation // output: first NumOfChecks elements show positions of dynamic frozen symbols
	, const unsigned* pBitReversalSortedArikan8 ///bit reversal permutation pBitReversal_m3 = GetBitReversal(m - 3)
	, const tBit* pCheckMatrix ///input check matrix
	, unsigned NumOfChecks ///number of check symbols
	, tBit* pT // temporary array of length (1u<<m)
	)
>>>>>>>
{
    tBit *pLastRow = pMatrix + (NumOfRows - 1) * NumOfAllocatedColumns;
    int c;
    for (c = NumOfColumns - 1; c >= 0; c--) {
        if (!pLastRow[c]) continue;
        unsigned rowToSubtract = (unsigned)(find(pLeadingColumns, pLeadingColumns + NumOfRows - 1, (unsigned)c) - pLeadingColumns);
        if (rowToSubtract != NumOfRows - 1) {
            XOR(pLastRow, pMatrix + rowToSubtract * NumOfAllocatedColumns, c + 1);
        }
    }
    c = NumOfColumns - 1;
    while ((c > 0) && (!pLastRow[c])) {
        c--;
    }
    if (!pLastRow[c]) {
        NumOfRows--;
        return false;
    }
    pLeadingColumns[NumOfRows - 1] = c;
    return true;
}

// print dynamic freezing constraints given matrix of them
void PrintFreezings(tBit* pMatrix,
    unsigned NumOfColumns,
    unsigned NumOfAllocatedColumns,
    unsigned NumOfRows,
    ostream &OutStream,
	unsigned* pColumnsOrder
    )
{
<<<<<<<
    tBit* pRow = pMatrix;
    for (unsigned i = 0; i < NumOfRows; i++) {
        unsigned w = 0;
        for (unsigned j = 0; j < NumOfColumns; j++)	w += ((pRow[j] != 0) ? 1 : 0);
        OutStream << std::endl << w;
		if (pColumnsOrder == nullptr)
=======
	unsigned* pNormalizedMatrix = (unsigned*)alloca(sizeof(unsigned)*m);
	unsigned* pPermutation = (unsigned*)alloca(sizeof(unsigned)*m);
	pNormalizedMatrix[0] = pMatrix[0] = 1;
	pPermutation[0] = 0;
	for (unsigned i = 0; i < m; i++)
	{
		unsigned Y, X;
		do
>>>>>>>
		{
<<<<<<<
			for (unsigned j = 0; j < NumOfColumns; j++)
			if (pRow[j]) OutStream << ' ' << j;
		}
		else
		{
			for (unsigned j = 0; j < NumOfColumns; j++)
			if (pRow[pColumnsOrder[j]]) OutStream << ' ' << pColumnsOrder[j];
		}
=======
			X = Y = 1 + gsl_rng_uniform_int(pRNG, (1u << m) - 1);
			//check if it is linearly dependent
			for (unsigned j = 0; j < i; j++)
				if ((X >> pPermutation[j]) & 1)
					X ^= pNormalizedMatrix[j];
		} while (!X);
		pNormalizedMatrix[i] = pMatrix[i] = X;
		for (unsigned j = 0; j < m; j++)
			if ((X >> j) & 1)
			{
				pPermutation[i] = j;
				for (unsigned s = 0; s < i; s++)
					if (pNormalizedMatrix[s] & (1u << j))
						pNormalizedMatrix[s] ^= X;
				break;
			};
>>>>>>>

        pRow += NumOfAllocatedColumns;
    }
};
*/
// print dynamic freezing constraints given matrix of them
void ReadFreezings(istream &SpecFile,
	unsigned N,
	unsigned k,
	unsigned* pDecisionPoints,
	unsigned** ppFreezingConstraints,
	unsigned &LineNum
)
{
	memset(pDecisionPoints, ~0, N*sizeof(*pDecisionPoints));
	// read code constraints
	for (unsigned i = 0; i < N-k; i++)
	{
		LineNum++;
		unsigned w;
		SpecFile >> w;
		if (!SpecFile)
			throw Exception("error: CDynFrozenEncoderXHong: File read error");
		if (/*(w > m_Dimension) || */!w)
			throw Exception("error: CDynFrozenEncoderXHong: Invalid freezing constraint in line %d", (LineNum));
		if (w == 1)
			ppFreezingConstraints[i] = nullptr;
		else
			ppFreezingConstraints[i] = new unsigned[w];
		unsigned F;
		for (unsigned j = 0; j < w - 1; j++)
			SpecFile >> ppFreezingConstraints[i][j];
		SpecFile >> F;
		if (!SpecFile)
			throw Exception("error: CDynFrozenEncoderXHong: File read error");
		if (F >= N)
		{
			delete[] ppFreezingConstraints[i];
			i--;
			continue;
		}

		if (w > 1)
			ppFreezingConstraints[i][w - 1] = ~0;
		if (pDecisionPoints[F] != ~0)
			throw Exception("error: CDynFrozenEncoderXHong: Duplicated freezing constraint on symbol %d", F);
		pDecisionPoints[F] = i;
	}
};

void RemoveFrozen(tBit* pH, unsigned k, unsigned n, unsigned N)
{
	if (N == 0)
	{
		N = n;
	}
	bool staticFrozen = true;
	while (staticFrozen)
	{
		staticFrozen = false;
		for (unsigned r = 0; r < k; r++)
		{
			unsigned w = 0;
			unsigned f = 0;
			for (unsigned c = 0; c < n; c++)
			{
				if (pH[r*N + c] != 0)
				{
					w++;
					f = c;
				}
			}
			if (w == 1)
			{
				for (unsigned r1 = 0; r1 < k; r1++)
				{
					if (r1 == r) continue;
					if (pH[r1*N + f] != BIT_0)
					{
						staticFrozen = true;	
						pH[r1*N + f] = BIT_0;
					}
				}
			}
		}
	}
};

// print dynamic freezing constraints given matrix of them
void PrintFreezings(tBit* pMatrix,
	unsigned NumOfColumns,
	unsigned NumOfAllocatedColumns,
	unsigned NumOfRows,
	ostream &OutStream
)
{
	tBit* pRow = pMatrix;
	for (unsigned i = 0; i < NumOfRows; i++)
	{
		unsigned w = 0;
		for (unsigned j = 0; j < NumOfColumns; j++)	w += (pRow[j] != 0) ? 1 : 0;
		OutStream << std::endl << w;
		for (unsigned j = 0; j < NumOfColumns; j++)
			if (pRow[j]) OutStream << ' ' << j;

		pRow += NumOfAllocatedColumns;
	}
};


tBit DotProduct(unsigned N, const tBit* pX, const tBit* pY)
{
	unsigned __int64 P = 0;
	unsigned N0 = N / sizeof(__int64);
	const unsigned __int64 *x = (const unsigned __int64 *)pX;
	const unsigned __int64 *y = (const unsigned __int64 *)pY;
	for (unsigned i = 0; i < N0; i++, x++, y++)
	{
		P ^= *x&*y;
	};
	for (unsigned i = N0 * sizeof(__int64); i < N; i++)
		P ^= pX[i] & pY[i];
	return _mm_popcnt_u64(P) & 1;

}


///compute a vectorized product of a vector and a matrix
///the elements of the vector are assumed to be blocks of size Stride
void MatrixMultiply(unsigned NumOfRows ///number of rows in a matrix
	, unsigned NumOfColumns ///number of columns in a matrix
	, unsigned Stride ///size of an input block in a matrix
	, const tBit* pInputVector///the input block vector 
	, const tBit* pMatrix ///the matrix
	, tBit* pOutputVector ///the output vector
	, bool ZeroOutput///true if pOutputVector must be zeroed. Otherwise, the result of the operation will be XORed with it
)
{
	if(ZeroOutput)
		memset(pOutputVector, 0, sizeof(tBit)*Stride*NumOfColumns);
	for (unsigned j = 0; j < NumOfRows; j++)
	{
		const tBit* pS = pInputVector + j*Stride;
		for (unsigned i = 0; i < NumOfColumns; i++)
		{
			if (pMatrix[j*NumOfColumns + i])
			{
				tBit* pD = pOutputVector + i*Stride;
				XOR(pD, pS, Stride);
			};
		}
	}

}

///compute a vectorized product of a vector and a transposed matrix
///the elements of the vector are assumed to be blocks of size Stride
void MatrixMultiplyTransposed(unsigned NumOfRows ///number of rows in a matrix
	, unsigned NumOfColumns ///number of columns in a matrix
	, unsigned Stride ///size of an input block in a matrix
	, const tBit* pInputVector///the input block vector 
	, const tBit* pMatrix ///the matrix
	, tBit* pOutputVector ///the output vector
	, bool ZeroOutput///true if pOutputVector must be zeroed. Otherwise, the result of the operation will be XORed with it
)
{
	if (ZeroOutput)
		memset(pOutputVector, 0, sizeof(tBit)*Stride*NumOfRows);
	for (unsigned j = 0; j < NumOfColumns; j++)
	{
		const tBit* pS = pInputVector + j*Stride;
		for (unsigned i = 0; i < NumOfRows; i++)
		{
			if (pMatrix[i*NumOfColumns + j])
			{
				tBit* pD = pOutputVector + i*Stride;
				XOR(pD, pS, Stride);
			};
		}
	}

}

//compute Fv inplace for a N*N matrix a
void MatrixVectorMultiplicationInplace(unsigned N, const tBit* pMatrix, unsigned Stride, tBit* pVector)
{
	tBit* pTemp = (tBit*)alloca(sizeof(tBit)*N);
	for (unsigned i = 0; i < N; i++, pMatrix += N)
	{
		tBit P = 0;
		for (unsigned j = 0; j < N; j++)
			P ^= pMatrix[j] & pVector[j*Stride];
		pTemp[i] = P;
	};
	for (unsigned i = 0; i < N; i++)
		pVector[i*Stride] = pTemp[i];
}

void MatrixKroneckerProductMultiplyTransposed(unsigned n, const tBit* pMatrix, unsigned KroneckerDegree, tBit* pVector)
{
	unsigned Stride = 1;
	unsigned N = pow(n, KroneckerDegree);
	while (N > Stride)
	{
		for (unsigned j = 0; j<Stride; j++)
		{
			for (unsigned i = 0; i < N; i += Stride*n)
			{
				MatrixVectorMultiplicationInplace(n, pMatrix, Stride, pVector + j + i);
			};
		};
		Stride *= n;
	};
};

// inverse a matrix
tBit* Inverse(tBit* pMatrix, unsigned NumOfRows)
{
	unsigned n = NumOfRows;
	tBit* pInverse = new tBit[n*n];
	tBit* P = new tBit[n*n * 2];
	memset(P, 0, n*n * 2 * sizeof(*P));
	for (unsigned i = 0; i < n; i++)
	{
		memcpy(P, pMatrix, n * sizeof(*P));
		P[n + i] = 1;
		P += 2 * n;
		pMatrix += n;
	}
	unsigned* pPerm = new unsigned[2 * n];
	for (unsigned i = 0; i < 2 * n; i++)
	{
		pPerm[i] = i;
	}
	P -= 2 * n*n;
	Gauss(P, n, 2 * n, 2 * n, true, pPerm);
	if (n != NumOfRows)
	{
		throw Exception("Error: inverse a singular matrix\n");
	}
	for (size_t i = 0; i < n; i++)
	{
		memcpy(pInverse, P + n, n * sizeof(*P));
		pInverse += n;
		P += 2 * n;
	}
	delete[] P;
};

