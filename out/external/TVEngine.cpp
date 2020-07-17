#include <iostream>
#include <algorithm>

#include "misc.h"
#include "TVEngine.h"
#include "TalVardyDecoder.h"
#include "LinAlg.h"
#include "SoftProcessing.h"

using namespace std;

CTVEngine::CTVEngine(unsigned MaxNumOfPaths, unsigned BitsPerBlock) : m_BitsPerBlock(BitsPerBlock), CTVMemoryEngine(MaxNumOfPaths)
{
    m_ppArrayPointer_C = nullptr;
}

CTVEngine_Prob::CTVEngine_Prob(unsigned MaxNumOfPaths) : CTVEngine(MaxNumOfPaths)
{
    m_ppArrayPointer_P = nullptr;
}

CTVEngine_S::CTVEngine_S(unsigned MaxNumOfPaths, unsigned BitsPerBlock) : CTVEngine(MaxNumOfPaths, BitsPerBlock)
{
    m_ppArrayPointer_S = nullptr;
#ifdef SSE2
    int64_t mask = 1ull << 63;
    SIGN_MASK = _mm_set_pd(
        *reinterpret_cast<double*>(&mask),
        *reinterpret_cast<double*>(&mask));
#endif
}

void CTVEngine::Init(unsigned LogLength, unsigned UpperLayer, unsigned MemoryLimit)
{
    CTVMemoryEngine::Init(LogLength, UpperLayer, MemoryLimit);

    m_ppArrayPointer_C = new tBit *[(m_UpperLayer + 1) * m_MaxNumOfPaths];

    for (unsigned lambda = 0; lambda <= UpperLayer; lambda++) {
        for (unsigned s = 0; s < m_MaxNumOfPaths; s++) {
#ifdef SMART_MEMORY
            m_ppArrayPointer_C[s * (m_UpperLayer + 1) + lambda] = nullptr;
#else
            m_ppArrayPointer_C[s*(m_UpperLayer+1)+ lambda] = m_Memory.Allocate<TV_TYPE_C>(BitsPerBlock << (m_UpperLayer - lambda));
#endif
        }
    }
}



CTVEngine::~CTVEngine()
{
#ifndef SPECIALIZED_DECODER_MODE
    if (m_ppArrayPointer_C) {
        for (unsigned lambda = 0; lambda <= m_UpperLayer; lambda++) {
            for (unsigned s = 0; s < m_MaxNumOfPaths; s++)
                m_Memory.Free(m_ppArrayPointer_C[s * (m_UpperLayer + 1) + lambda]);
        }
    }
    delete[] m_ppArrayPointer_C;
#endif
}

CTVEngine_Prob::~CTVEngine_Prob()
{
#ifndef SPECIALIZED_DECODER_MODE
    if (m_ppArrayPointer_P) {
        for (unsigned lambda = 0; lambda <= m_UpperLayer; lambda++) {
            for (unsigned s = 0; s < m_MaxNumOfPaths; s++)
                m_Memory.Free(m_ppArrayPointer_P[lambda * m_MaxNumOfPaths + s]);
        }
    }
    delete[]m_ppArrayPointer_P;
#endif
}

CTVEngine_S::~CTVEngine_S()
{
#ifndef SPECIALIZED_DECODER_MODE
    if (m_ppArrayPointer_S) {
        for (unsigned lambda = 0; lambda <= m_UpperLayer; lambda++) {
            for (unsigned s = 0; s < m_MaxNumOfPaths; s++) {
                m_Memory.Free(m_ppArrayPointer_S[s*(m_UpperLayer + 1) + lambda]);
            }
        }
    }
    delete[] m_ppArrayPointer_S;
#endif
}

tPairOfDouble* CTVEngine_Prob::GetArrayPointer_P(unsigned lambda,//layer id
    unsigned l //path index
    )
{

    unsigned* pIndexArray = GetIndexArrayPointer(l);
    unsigned p = UpdateIndex(lambda, pIndexArray);
    return m_ppArrayPointer_P[p];
}

const tPairOfDouble* CTVEngine_Prob::CGetArrayPointer_P(unsigned lambda,//layer id
    unsigned l //path index
    ) const
{
    unsigned* pIndexArray = GetIndexArrayPointer(l);
    unsigned p = pIndexArray[lambda];
    return m_ppArrayPointer_P[p];
}

#ifdef SMART_MEMORY
void CTVEngine::AllocateEntry(unsigned lambda, unsigned s1)
{
    // For this AllocateEntry we do not need to check if pointer is already allocated? Why?
    //#ifndef USE_MEMORY_MANAGER
    //    if (m_ppArrayPointer_C[s1] != nullptr)
    //        return;
    //#endif
    m_ppArrayPointer_C[s1] = m_Memory.Allocate<TV_TYPE_C>(m_BitsPerBlock << (m_LogLength - lambda));
}

void CTVEngine_Prob::AllocateEntry(unsigned lambda, unsigned s1)
{
#ifndef USE_MEMORY_MANAGER
    if (m_ppArrayPointer_P[s1] != nullptr)
        return;
#endif
    m_ppArrayPointer_P[s1] = m_Memory.Allocate<TV_TYPE_P>(m_BitsPerBlock << (m_LogLength - lambda));
    CTVEngine::AllocateEntry(lambda, s1);
}

void CTVEngine_S::AllocateEntry(unsigned lambda, unsigned s1)
{
#ifndef USE_MEMORY_MANAGER
    if (m_ppArrayPointer_S[s1] != nullptr)
        return;
#endif
#ifdef MEMORY_COUNTING
    m_pBlocks[lambda]++;
#endif

    m_ppArrayPointer_S[s1] = m_Memory.Allocate<TV_TYPE_S>(m_BitsPerBlock << (m_LogLength - lambda));
    CTVEngine::AllocateEntry(lambda, s1);
}
#endif

//allocate one temporary array at a given layer
unsigned CTVEngine_S::GetTempIndex(unsigned lambda)
{
    unsigned s1 = Pop(m_pInactiveArrayIndices + lambda*(m_MaxNumOfPaths + 1));
#ifdef SMART_MEMORY
    AllocateEntry(lambda, s1);
#endif
    return s1;
}

//allocate one temporary array at a given layer
void CTVEngine_S::ReleaseTemp(unsigned lambda, unsigned s)
{
    Push(s, m_pInactiveArrayIndices + lambda*(m_MaxNumOfPaths + 1));
}

//recover a dynamic frozen symbol from the interim variables
tBit CTVEngine::EvaluateDynFrozen(unsigned l ///Path ID
    , unsigned NumOfTerms ///number of terms in the expression
    , const Triple* pExpression ///list of terms to be used for evaluation
    )
{
    tBit CurValue = 0;
    const tBit* pC0 = 0;
    unsigned lambda0 = ~0u;
    const unsigned* pArrayPtr = GetIndexArrayPointer(l);
    for (unsigned i = 0; i<NumOfTerms; i++) {
        const Triple&T = pExpression[i];
        if (T.lambda != lambda0) {
            pC0 = CGetArrayPointer_C(T.lambda, pArrayPtr, 0);
            lambda0 = T.lambda;
        };
        //const tBit* pC1=pC0+(T.phi << (m_LogLength - T.lambda));
        _ASSERT((pC0[T.beta] & 0x7f)==0);
        CurValue ^= pC0[T.beta];
    }
    return CurValue;
}


#if defined(PACKED_DF_EVALUATION)
//recover a dynamic frozen symbol using packed expression
tBit CTVEngine::EvaluateDynFrozenPacked(
    unsigned l ///path ID
    , const ReverseDFExpression& Expression ///list of terms to be used for evaluation
    )
{
    if (!Expression.pExpression)
        return 0;
    const tBit* pC0 = 0;
    unsigned lambda0 = ~0u;
    const unsigned* pArrayPtr = GetIndexArrayPointer(l);
    tDFMask Z = m_pDFMask[l]; 
    const Triple* pSrc = Expression.pExpression;
    const Triple* pLast= Expression.pExpression+Expression.NumOfTerms;
    for (pSrc; pSrc<pLast; pSrc++)
    {
        if (pSrc->lambda != lambda0)
        {
            pC0 = CGetArrayPointer_C(pSrc->lambda, pArrayPtr, 0);
            lambda0 = pSrc->lambda;
        };
        if (pC0[pSrc->beta])
            Z ^= pSrc->phi;
    };
    tBit Res = Z& 1;
    m_pDFMask[l] = Z>>1;
    return Res;
}
#endif
//#define JACOBI_ZERO

//#define RATIONAL_APPROX1
#define RATIONAL_APPROX2
//compute log(1+exp(-x)), x>=0;
inline double JacobiLog(double x)
{
#if defined(JACOBI_ZERO)
    return 0;
#elif defined(RATIONAL_APPROX1)
    if (x<4.0)
        return 1.343 / (1.405 + x) - 0.25;
    else
        return 0;
#elif defined (RATIONAL_APPROX2)
    if (x<5.5)
        return .131478060787483 + (3.43245248951091 - 1.74351875423031*x) / (6.11620503468962 + (2.21612896897222 + x)*x);
    else
        return 0;
#else
    return log(1 + exp(-x));
#endif
}



void CTVEngine_Prob::RecursivelyCalcP(unsigned l,//path ID
    unsigned lambda,//layer ID
    unsigned phi //phase
    )
{
    if (lambda == 0)
        return;
    unsigned psi = phi / 2;
    //recurse first, if needed
    if (!(phi & 1)) RecursivelyCalcP(l, lambda - 1, psi);
    //do the calculation
    //    double sigma=0;
    tPairOfDouble* P_l = GetArrayPointer_P(lambda, l);
    const tPairOfDouble* P_l_1 = CGetArrayPointer_P(lambda - 1, l);
    const tBit* C_l = CGetArrayPointer_C(lambda, l, 0);
    unsigned N = (m_BitsPerBlock << (m_LogLength - lambda));
    for (unsigned beta = 0; beta<N; beta++) {
        if ((phi & 1)) {
            unsigned u1 = C_l[beta]>0;
            P_l[beta].P[0] = P_l_1[beta].P[u1] + P_l_1[beta + N].P[0];
            P_l[beta].P[1] = P_l_1[beta].P[u1 ^ 1] + P_l_1[beta + N].P[1];

        } else {
            if (P_l_1[beta].P[1] < -INF_BOUND) {
                P_l[beta].P[0] = P_l_1[beta].P[0] + P_l_1[beta + N].P[0];
                P_l[beta].P[1] = P_l_1[beta].P[0] + P_l_1[beta + N].P[1];
            } else if (P_l_1[beta].P[0] < -INF_BOUND) {
                P_l[beta].P[0] = P_l_1[beta].P[1] + P_l_1[beta + N].P[1];
                P_l[beta].P[1] = P_l_1[beta].P[1] + P_l_1[beta + N].P[0];
            } else if (P_l_1[beta + N].P[1] < -INF_BOUND) {
                P_l[beta].P[0] = P_l_1[beta].P[0] + P_l_1[beta + N].P[0];
                P_l[beta].P[1] = P_l_1[beta].P[1] + P_l_1[beta + N].P[0];
            } else if (P_l_1[beta + N].P[0] < -INF_BOUND) {
                P_l[beta].P[0] = P_l_1[beta].P[1] + P_l_1[beta + N].P[1];
                P_l[beta].P[1] = P_l_1[beta].P[0] + P_l_1[beta + N].P[1];
            } else {
                for (unsigned u1 = 0; u1 < 2; u1++) {
                    //max-log-map implementation
                    double X[2];
                    for (unsigned u2 = 0; u2<2; u2++)
                        X[u2] = P_l_1[beta].P[u1^u2] + P_l_1[beta + N].P[u2];
                    if (X[0]>X[1])
                        P_l[beta].P[u1] = X[0] + JacobiLog(X[0] - X[1]);
                    else
                        P_l[beta].P[u1] = X[1] + JacobiLog(X[1] - X[0]);
                };
            }
            /*double p = P_l_1[2*beta].P[0]+P_l_1[2*beta+1].P[0];
            double q = P_l_1[2*beta].P[1]+P_l_1[2*beta+1].P[1];
            if (p>q) P_l[beta].P[0] = p; else P_l[beta].P[0] = q;
            p = P_l_1[2*beta].P[0]+P_l_1[2*beta+1].P[1];
            q = P_l_1[2*beta].P[1]+P_l_1[2*beta+1].P[0];
            if (p>q) P_l[beta].P[1] = p; else P_l[beta].P[1] = q;*/
        }
    }
}

/*const unsigned int pMasks[64] = {
    0x1u <<  0, 0x1u <<  1, 0x1u <<  2, 0x1u <<  3, 0x1u <<  4,
    0x1u <<  5, 0x1u <<  6, 0x1u <<  7, 0x1u <<  8, 0x1u <<  9,
    0x1u << 10, 0x1u << 11, 0x1u << 12, 0x1u << 13, 0x1u << 14,
    0x1u << 15, 0x1u << 16, 0x1u << 17, 0x1u << 18, 0x1u << 19,
    0x1u << 20, 0x1u << 21, 0x1u << 22, 0x1u << 23, 0x1u << 24,
    0x1u << 25, 0x1u << 26, 0x1u << 27, 0x1u << 28, 0x1u << 29,
    0x1u << 30, 0x1u << 31
    };
    */
MType* CTVEngine_S::IterativelyCalcS(unsigned* pIndexArray//path index array
    , unsigned lambda//layer ID
    , int Depth //number of  iterations to do
    , bool ShouldCombine //true if SoftCombine should be called
    )
{
    unsigned lambda0 = lambda - Depth;
    if (lambda0 == 0)
        return GetArrayPointer_S(lambda0, pIndexArray);

    const MType *M_l_1 = CGetArrayPointer_S(lambda0 - 1, pIndexArray);
    unsigned N = (m_BitsPerBlock << (m_LogLength - lambda0));
    MType *M_l = nullptr;
    if (ShouldCombine) {
        const tBit* C_l = CGetArrayPointer_C(lambda0, pIndexArray, 0);
        M_l = GetArrayPointer_S(lambda0, pIndexArray);
        SoftCombine(M_l, M_l_1, M_l_1 + N, C_l, N);
        lambda0++;
        M_l_1 = M_l;
        N >>= 1;
    }

    for (lambda0; lambda0<= lambda; lambda0++) {
        //do the calculation
        M_l = GetArrayPointer_S(lambda0, pIndexArray);
        SoftXOR(M_l, M_l_1, M_l_1 + N, N);
        M_l_1 = M_l;
        N >>= 1;
    }
    return M_l;
}

void CTVEngine::IterativelyUpdateC_Prim(unsigned*pIndexArray, unsigned lambda, unsigned Depth)
{
    if (Depth == 0)
        return;
    unsigned lambda0 = lambda - Depth;
    const tBit* pC0 = CGetArrayPointer_C(lambda, pIndexArray, 0);
    //write everything to its ultimate destination
    tBit* pC_Out = GetArrayPointer_C(lambda0, pIndexArray, 0);
    unsigned N = (m_BitsPerBlock << (m_LogLength - lambda));
    unsigned FinalLength = N << Depth;

    const tBit* pC1 = pC_Out + FinalLength - N;
    pC_Out += FinalLength - 2 * N;
    XOR(pC_Out, pC0, pC1, N);
    lambda--;
    for (lambda; lambda > lambda0; lambda--) {
        N <<= 1;
        pC1 = pC_Out;
        pC_Out -= N;
        const tBit* pC0 = CGetArrayPointer_C(lambda, pIndexArray, 0);
        _ASSERT((pC0[0] & 0x7f) == 0);
        _ASSERT((pC1[0] & 0x7f) == 0);
        XOR(pC_Out, pC0, pC1, N);
    }
}

void CTVEngine_Prob::Init(unsigned LogLength, unsigned UpperLayer, unsigned MemoryLimit)
{
    CTVEngine::Init(LogLength, UpperLayer, MemoryLimit);
    m_ppArrayPointer_P = new tPairOfDouble*[(m_UpperLayer + 1)*m_MaxNumOfPaths];
    for (unsigned lambda = 0; lambda <= UpperLayer; lambda++) {
        for (unsigned s = 0; s < m_MaxNumOfPaths; s++) {
#ifdef SMART_MEMORY
            m_ppArrayPointer_P[lambda*m_MaxNumOfPaths + s] = nullptr;
#else
            m_ppArrayPointer_P[lambda*m_MaxNumOfPaths + s] = m_Memory.Allocate<TV_TYPE_P>(1u << (m_LogLength - lambda));
#endif
        }
    }
}


void CTVEngine_S::Init(unsigned LogLength, unsigned UpperLayer, unsigned MemoryLimit)
{
    CTVEngine::Init(LogLength, UpperLayer, MemoryLimit);
    m_ppArrayPointer_S = new MType *[(m_UpperLayer + 1) * m_MaxNumOfPaths];

    for (unsigned lambda = 0; lambda <= UpperLayer; lambda++) {
        for (unsigned s = 0; s < m_MaxNumOfPaths; s++) {
#ifdef SMART_MEMORY
            m_ppArrayPointer_S[s * (m_UpperLayer + 1) + lambda] = nullptr;
#else
            m_ppArrayPointer_S[s * (m_UpperLayer + 1) + lambda] = m_Memory.Allocate<TV_TYPE_S>(BitsPerBlock << (m_LogLength - lambda));
#endif
        }
    }
}
