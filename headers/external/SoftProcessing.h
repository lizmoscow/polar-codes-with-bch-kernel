#ifndef SOFTPROCESSING_H
#define SOFTPROCESSING_H

#include "SeqConfig.h"
#include "misc.h"

/// Q(a, b) = sign(a)*sign(b)*min(|a|,|b|)
void SoftXOR(MType *M_l, const MType *M_l_1, const MType *M_l_2, const uint32_t &N);
inline void SoftXOR(MType *M_l, const MType *M_l_1, const uint32_t &N)
{
	SoftXOR(M_l, M_l_1, M_l_1 + N, N);
};
/// Calculation of Q function:
/// Q(a, b, c, d) = sign(a)*sign(b)*sign(c)*sign(d)*min(|a|,|b|,|c|,|d|), as well as sign(a)*sign(c)*min(|a|,|c|), sign(b)*sign(d)*min(|b|,|d|)
/// a,b,c,d are assumed to be located in the same array in consecutive blocks of size n
void SoftXOR(float *ABCD, float* AC, const float *A, const uint32_t &N);


/// Calculation of Q function:
/// Q(a, b) = sign(a)*sign(b)*min(a,b)
/// return \sum_i |Q(a_i,b_i)|
KType SoftXORWithSum(MType *M_l, const MType *M_l_1, const MType *M_l_2, const uint32_t &N);

/// Calculation of G function: 
/// g(lambda, beta) = (-1)^C_l[beta][0]*M_{l-1}[2*beta] + M_{l-1}[2*beta+1]
void SoftCombine(MType *M_l, const MType *M_l_1, const MType *M_l_2, const tBit *C_l, const uint32_t &N);
inline void SoftCombine(MType *M_l, const MType *M_l_1, const tBit *C_l, const uint32_t &N)
{
	SoftCombine(M_l, M_l_1, M_l_1 + N, C_l, N);
};

void SoftCombine2Elements(MType *M_l, const MType &M_l_1, const MType &M_l_2, const tBit &C_l);
void SoftXOR2Elements(MType *M_l, const MType &M_l_1, const MType &M_l_2);

/// Calculation of G function: 
/// g(lambda, beta) = (-1)^C_l[beta][0]*M_{l-1}[2*beta] + M_{l-1}[2*beta+1]
/// return \sum_beta	 |g(lambda,beta)|
MType SoftCombineWithSum(MType *M_l, const MType *M_l_1, const MType*M_l_2, const tBit *C_l, const uint32_t &N);

///compute \sum_{i=0} |L_i|
KType SumAbs(const MType* pSrc, unsigned N);

//change LLR signs using a given mask
//if C[i]=1, L[i]=-L[i]
void ChangeLLRSigns(MType* pDestLLRs, const MType* pSrcLLRs, const tBit* pMask, unsigned  N);


void SoftXOREvenOdd(MType *M_l, const MType *M_l_1, const uint32_t &N);
void SoftXOR_0312(MType *M_l, const MType *M_l_1, const uint32_t &N);
void SoftCombineEvenOdd(MType *M_l, const MType *M_l_1, const tBit *C_l, const uint32_t &N);

#endif