#include "StructuredLinAlg.h"
#include <iostream>
#include "SoftSIMD.h"
#include "SoftProcessing.h"

Word OneArray[65];
unsigned ShortOneArray[33];
bool OneBoolArray[65][8];
bool isOneArrayInitialized = false;
#ifdef AVX
__m256 SignMask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF));
__m256 Mask2 = _mm256_castsi256_ps(_mm256_set1_epi32(0xFFFFFFFF));
#endif

const __m128i Shuffle0 = _mm_set_epi8((char)0xFF, (char)0xFF, (char)0xFF, (char)0xFF, (char)0xFF, (char)0xFF, (char)0xFF, (char)0xFF,
	(char)0xFF, (char)0xFF, (char)0xFF, (char)0xFF, (char)15, (char)11, (char)7, (char)3);
const __m128i Shuffle1 = _mm_set_epi8((char)0xFF, (char)0xFF, (char)0xFF, (char)0xFF, (char)0xFF, (char)0xFF, (char)0xFF, (char)0xFF,
	(char)15, (char)11, (char)7, (char)3, (char)0xFF, (char)0xFF, (char)0xFF, (char)0xFF);


unsigned CalculateSyndrome(const unsigned* pCheckMatrix, const MType* LLRs, unsigned N, tBit* pC)
{
#if defined(AVX) && defined(USE_FLOAT)
	const MType* pLLR = LLRs;
	tBit* pHardDecision = pC;
	unsigned syndrome = 0;
	unsigned W = 0;
	unsigned i = 0;
	for (i = 0; i < (N&~7); i += 8)
	{
		__m256 L = _mm256_load_ps(pLLR + i);
		__m256 Sign = _mm256_andnot_ps(SignMask, L);//extract signs
		__m128i Sign0 = _mm_castps_si128(_mm256_extractf128_ps(Sign, 0));//unpack them
		__m128i Sign1 = _mm_castps_si128(_mm256_extractf128_ps(Sign, 1));//unpack them
		Sign0 = _mm_shuffle_epi8(Sign0, Shuffle0);
		Sign1 = _mm_shuffle_epi8(Sign1, Shuffle1);
		Sign1 = _mm_or_si128(Sign0, Sign1);
		*(unsigned __int64*)(pHardDecision + i) = _mm_cvtsi128_si64(Sign1);
		//calculate hard decisions
		__m256 cmp = _mm256_maskload_ps((const float*)(pCheckMatrix + i), _mm256_castps_si256(Sign));
		cmp = _mm256_xor_ps(cmp, _mm256_shuffle_ps(cmp, cmp, 0x0e));
		cmp = _mm256_xor_ps(cmp, _mm256_shuffle_ps(cmp, cmp, 0x01));
		__m256i res = _mm256_castps_si256(cmp);
		syndrome ^= res.m256i_u32[0] ^ res.m256i_u32[4];
	}

	return syndrome;
#else
	unsigned Syndrome = 0;
	bool isParity = !GetHardDecisionsAndParity(LLRs, N, pC);
	if (isParity){
		for (unsigned i = 0; i < N; ++i){
			if (pC[i])
				Syndrome ^= pCheckMatrix[i];
		}
		return Syndrome;
	}
	else{
		return 1;
	}
#endif
}

KType ComputeSumAbsLLRs(const MType *pLLR, unsigned n)
{
	KType result = 0;
    uint32_t i = 0;
	SUM_BLOCK_COUNT(n - 1);
#if defined(AVX)
#if defined(USE_FLOAT)
    __m256 llr;
	while (i < (n & (~0x7)))
    {
        llr = _mm256_load_ps(pLLR + i);
		llr = _mm256_and_ps(llr, SignMask);
        llr = _mm256_hadd_ps(llr, llr);
        llr = _mm256_hadd_ps(llr, llr);
        result += llr.m256_f32[0] + llr.m256_f32[4];
        i += 8;
	}
#elif defined(USE_INT8)
    __m128i LLR_128;
    while (i < (n & (~0x7)))  // process 8 llrs
    {
        LLR_128 = _mm_load_si128(reinterpret_cast<const __m128i*>(pLLR + i));
        LLR_128 = _mm_cvtepi8_epi16(LLR_128); // first 8 llrs
        LLR_128 = _mm_abs_epi16(LLR_128);
        LLR_128 = _mm_hadd_epi16(LLR_128, LLR_128);           // 8 -> 4
        LLR_128 = _mm_hadd_epi16(LLR_128, LLR_128);           // 4 -> 2
        result += LLR_128.m128i_i16[0] + LLR_128.m128i_i16[1];
        i += 8;
    }
#endif
#endif
	while (i < n)
	{
		result += abs_value(pLLR[i]);
        ++i;
	}
	return result;
}
//convert a matrix into the systematic form on the positions given by pColumnPermutation[i].second
void DiagonalizeMatrix(Word* pGeneratorMatrix, unsigned Length, unsigned Dimension, tPairOfMTypeUnsigned* pColumnPermutation){

	if (Length > 64){
		throw Exception("Diagonalize Matrix :Length Must be <= 64");
	}

	for (unsigned i = 0; i < Dimension; i++)
	{
		Word& pCurRow = pGeneratorMatrix[i];
		//find a column with 1 below the staircase
		unsigned Column;
		for (unsigned c = i; c < Length; c++)
		{
			bool Finished = false;
			Column = pColumnPermutation[c].second;
			for (unsigned j = i; j<Dimension; j++)
			{
				const Word pRow = pGeneratorMatrix[j];
				if (pRow & OneArray[Column])
				{
					std::swap(pColumnPermutation[i], pColumnPermutation[c]);
					if (j>i) pCurRow ^= pRow;

					Finished = true;
					break;
				};
			};
			if (Finished) break;
		};
		//eliminate everything below the diagonal
		for (unsigned j = i + 1; j<Dimension; j++)
		{
			Word& pRow = pGeneratorMatrix[j];
			if (pRow & OneArray[Column]) pRow ^= pCurRow;
		};
		//eliminate everything above the diagonal
		for (unsigned j = 0; j<i; j++)
		{
			Word& pRow = pGeneratorMatrix[j];
			if (pRow & OneArray[Column]) pRow ^= pCurRow;
		};
	};
};

//construct a generator matrix in packed form
void GetGeneratorMatrix(Word* pGenMatrix, unsigned Length, unsigned Dimension, tBit* pGeneratorMatrix)
{
	if (!isOneArrayInitialized){
		//????? why i=64
		for (int i = 0; i <= 64; ++i){
			OneArray[i] = ONE << i;
		}
		isOneArrayInitialized = true;
	}
	unsigned index = 0;
	Word* pDest = pGenMatrix;
	for (unsigned i = 0; i < Dimension; ++i){
		Word tmp = 0;
		Word two = 1;
		for (unsigned j = 0; j < Length; ++j){
			tmp += two * ((pGeneratorMatrix[i * Length + j] == BIT_1)? 1 : 0);
			pDest[index] = tmp;
			two <<= 1;
		}
		index++;
	}
};
void ComputeWeights(Word *pDecisionsBuffer, unsigned BufferSize, unsigned Length, KType *pScores, KType *pAbsLLRs, unsigned *pShortDecisions)
{
	unsigned j, i, k;
	unsigned bound = 1;
#ifdef AVX
	bound = AVXBound;
#else
	bound = SSEBound;
#endif
	for (j = 0; j < BufferSize; j += bound)
	{
#ifdef USE_FLOAT
		if (pShortDecisions)
		{

#ifndef AVX
			__m128i mask;
			__m128 scores = _mm_setzero_ps();
			__m128i zeroReg = _mm_setzero_si128();
			__m128 tmpScore;
			__m128i var = _mm_load_si128((__m128i*)(pShortDecisions + j));

			CMP_BLOCK_COUNT(Length);	
			SUM_BLOCK_COUNT(Length);
			for (i = 0; i < Length; ++i)
			{
				mask = _mm_set1_epi32((int)OneArray[i]);
				mask = _mm_and_si128(mask, var);
				mask = _mm_cmpeq_epi32(mask, zeroReg);
				tmpScore = _mm_set_ps1(pAbsLLRs[i]);
				tmpScore = _mm_blendv_ps(tmpScore, _mm_castsi128_ps(zeroReg), _mm_castsi128_ps(mask));
				scores = _mm_sub_ps(scores, tmpScore);
			}
			_mm_store_ps(pScores + j, scores);
#else
			__m256 oneReg = _mm256_set1_ps(1.);
			__m256 mask;
			__m256 scores = _mm256_setzero_ps();
			__m256 tmpScore;
			__m256i var = _mm256_load_si256((__m256i*)(pShortDecisions + j));

			CMP_BLOCK_COUNT(Length);	
			SUM_BLOCK_COUNT(Length);
			for (i = 0; i < Length; ++i)
			{
				mask = _mm256_castsi256_ps(_mm256_set1_epi32((int)OneArray[i]));
				mask = _mm256_and_ps(mask, _mm256_castsi256_ps(var));
				mask = _mm256_xor_ps(mask, oneReg);
				mask = _mm256_cmp_ps(mask, oneReg, _CMP_EQ_UQ);
				tmpScore = _mm256_set1_ps(pAbsLLRs[i]);
				tmpScore = _mm256_blendv_ps(tmpScore, _mm256_setzero_ps(), mask);
				scores = _mm256_sub_ps(scores, tmpScore);
			}
			_mm256_store_ps(pScores + j, scores);
#endif
		}
		else
		{
			for (k = 0; k < SSEBound; ++k)
			{
				pScores[j + k] = 0;
				CMP_BLOCK_COUNT(Length);	
				SUM_BLOCK_COUNT(Length);
				for (i = 0; i < Length; ++i)
				{
					if (pDecisionsBuffer[j + k] & (int)OneArray[i])
					{
						pScores[j + k] -= pAbsLLRs[i];
					}
				}
			}
		}
#else
		//Stuff for integer score values
		//Intermediate solution:
		for (k = 0; k < SSEBound; ++k)
		{
			CMP_BLOCK_COUNT(Length);
			SUM_BLOCK_COUNT(Length);
			pScores[j + k] = 0;
			for (i = 0; i < Length; ++i)
			{
				if (pDecisionsBuffer[j + k] & OneArray[i])
				{
					pScores[j + k] -= pAbsLLRs[i];
				}
			}
		}
#endif
	}
}

void findLRIndexes(int row, tBit* matrix, int n, unsigned& l, unsigned& r){
	for (int i = 0; i < n; ++i){
		if (matrix[row*n + i] == BIT_1){
			l = i;
			break;
		}
	}

	for (int i = n - 1; i >= 0; --i){
		if (matrix[row*n + i] == BIT_1){
			r = i;
			break;
		}
	}
}

void xorMatrixRow(int rowA, int rowB, tBit* matrix, int n){
	for (int i = 0; i < n; ++i){
		matrix[n*rowA + i] ^= matrix[n*rowB + i];
	}
}

using std::cout;
using std::endl;

void getMSGM(unsigned n, unsigned k, tBit* matrix)
{
	//position of leftmost nonzero entry
	unsigned* L = new unsigned[k];
	//position of rightmost nonzero entry
	unsigned* R = new unsigned[k];

	for (unsigned i = 0; i < k; ++i){
		findLRIndexes(i, matrix, n, L[i], R[i]);
	}

	bool foundFlag = false;

	while (!foundFlag){
		foundFlag = true;
		for (unsigned i = 0; i < k; ++i){
			for (unsigned j = 0; j < k; ++j){
				if (i != j){
					if (((L[i] == L[j]) && (R[i] <= R[j])) || ((L[i] >= L[j]) && (R[i] == R[j]))){
						foundFlag = false;
						//XOR(matrix + n*j, matrix + n*i, n);
						xorMatrixRow(j, i, matrix, n);
						findLRIndexes(j, matrix, n, L[j], R[j]);
					}
				}
			}
		}
	}

	delete[] L;
	delete[] R;
}
