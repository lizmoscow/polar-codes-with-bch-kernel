#include <algorithm>
#include "TrellisKernelProcessor.h"
#include "LinAlg.h"

using namespace std;
typedef unsigned __int64 Word;
void MinimumSpan(unsigned K, unsigned N, Word* pMatrix, unsigned* pRowStart, unsigned* pRowEnd)
{
	memset(pRowStart, ~0, sizeof(unsigned)*N);
	memset(pRowEnd, ~0, sizeof(unsigned)*N);
	unsigned C = 0;

	for (unsigned i = 0; i < K; i++)
	{
		bool Found = false;
		for (; C < N; C++)
		{
			if (!((pMatrix[i]>>C)&1))
			{
				for(unsigned j=i+1;j<K;j++)
					if ((pMatrix[j]>>C) & 1)
					{
						pMatrix[i] ^= pMatrix[j];
						Found = true;
						break;
					}
				if (Found)
				{
					pRowStart[C] = i;
					break;
				}
			}
			else
			{
				pRowStart[C] = i;
				Found = true;
				break;
			};
		};
		if (!Found)
			throw Exception("Matrix is not full rank");
		for (unsigned j = i + 1; j < K; j++)
		{
			if ((pMatrix[j]>>C) & 1ull)
				pMatrix[j] ^= pMatrix[i];
		}

	};
	for (int i = K - 1; i >= 0; i--)
	{
		//find the last non-zero entry
		for (int j = N - 1; j >= 0; j--)
		{
			if ((pMatrix[i] >> j) & 1)
			{
				pRowEnd[j] = i;
				for (int s = 0; s < i; s++)
				{
					if ((pMatrix[s] >> j) & 1)
						pMatrix[s] ^= pMatrix[i];

				}
				break;
			}
		}
	}
}

CTrellisKernelProcessor::CTrellisKernelProcessor(const CMatrixBinaryKernel& K):CKernProcLLR(K)
{
	if (K.Size() >= 64)
		throw Exception("Kernel is too big for the trellis processor");
	m_pppTrellis = new TrellisNode**[K.Size()];
	Word* pExtGM = new Word[K.Size()];
	unsigned N = K.Size() + 1;
	const tBit* pK = K.GetKernel();
	unsigned* pRowStart = new unsigned[N];
	unsigned* pRowEnd = new unsigned[N];
	unsigned* pActiveBits= new unsigned[K.Size()];///the list of active input bits
	unsigned __int64* pCW0 = new unsigned __int64[1ull<<(K.Size() / 2 + 1)];
	unsigned __int64* pCW1 = new unsigned __int64[1ull << (K.Size() / 2 + 1)];
	m_ppNumOfActiveBits = new unsigned*[K.Size()];
	m_MaxNumOfAB = 0;
	_ASSERT(_CrtCheckMemory());
	for (unsigned i = 0; i < K.Size(); i++)
	{
		m_ppNumOfActiveBits[i] = new unsigned[N];
		//construct an extended generator matrix
		for (unsigned j = i; j < K.Size(); j++)
		{
			const tBit* pSrc = pK + j*K.Size();
			Word& Row = pExtGM[j - i];
			Row = 0;
			for (unsigned s = 0; s < K.Size(); s++)
				if (pSrc[s])
					Row |= 1ull << s;
		};
		pExtGM[0] |= 1ull << K.Size();
		_ASSERT(_CrtCheckMemory());
		MinimumSpan(K.Size() - i, N, pExtGM, pRowStart, pRowEnd);
		_ASSERT(_CrtCheckMemory());
		m_pppTrellis[i] = new TrellisNode*[N];
		unsigned NumOfActiveBits = 0;
		pCW0[0] = 0;
		m_ppNumOfActiveBits[i][0] = 0;
		_ASSERT(_CrtCheckMemory());
		for (unsigned j = 0; j <= K.Size(); j++)
		{
			size_t B = find(pActiveBits, pActiveBits + NumOfActiveBits, pRowEnd[j]) - pActiveBits;
			unsigned long long ExtractionMask = (pRowEnd[j]==~0)?~0ull:( (1ull << B) - 1);
			m_pppTrellis[i][j] = new TrellisNode[1ull << NumOfActiveBits];
			if (pRowStart[j] == ~0)
			{
					//no new rows are starting here
					for (unsigned long long State = 0; State < 1ull << NumOfActiveBits; State++)
					{
						Word NextBit = (pCW0[State] >> j) & 1;
						unsigned long long NextState = (State&ExtractionMask) | ((State >> 1)&~ExtractionMask);
						pCW1[NextState] = pCW0[State];
						m_pppTrellis[i][j][State].NextState[NextBit] = NextState;
						m_pppTrellis[i][j][State].NextState[1 - NextBit] = ~0;
					};
					_ASSERT(_CrtCheckMemory());
			}
			else
			{

					for (unsigned long long State = 0; State < 1ull << NumOfActiveBits; State++)
					{
						unsigned long long NextState0 = State;
						unsigned long long NextState1 = State^(1ull<<NumOfActiveBits);
						//compress the state variables
						NextState0 = (NextState0 &ExtractionMask) | ((NextState0 >> 1)&~ExtractionMask);
						NextState1 = (NextState1 &ExtractionMask) | ((NextState1 >> 1)&~ExtractionMask);
						Word C1= pCW0[State] ^ pExtGM[pRowStart[j]];
						pCW1[NextState0] = pCW0[State];
						pCW1[NextState1] = C1;
						Word NextBit0 = (pCW0[State] >> j) & 1;
						Word NextBit1 = (C1 >> j) & 1;
						_ASSERT((NextBit0^NextBit1) == 1);
						m_pppTrellis[i][j][State].NextState[NextBit0] = NextState0;
						m_pppTrellis[i][j][State].NextState[NextBit1] = NextState1;
						_ASSERT(_CrtCheckMemory());
					}
					pActiveBits[NumOfActiveBits++] = pRowStart[j];
					_ASSERT(_CrtCheckMemory());
			};
			swap(pCW0, pCW1);
			if (pRowEnd[j] != ~0)
			{
				memmove(pActiveBits + B, pActiveBits + B + 1, sizeof(unsigned)*(NumOfActiveBits - B));
				NumOfActiveBits--;
			};
			m_ppNumOfActiveBits[i][j+1] = NumOfActiveBits;
			if (NumOfActiveBits > m_MaxNumOfAB)
				m_MaxNumOfAB = NumOfActiveBits;
			_ASSERT(_CrtCheckMemory());
		};
	/*	cout << "digraph G{\n";
		for (unsigned j = 0; j < N; j++)
		{
			for (unsigned s = 0; s < 1u << m_ppNumOfActiveBits[i][j]; s++)
			{
				if (m_pppTrellis[i][j][s].NextState[0] != ~0)
					cout << "S_" << j << '_' << s << "->" << "S_" << (j + 1) << '_' << m_pppTrellis[i][j][s].NextState[0] << "[style=dashed];\n";
				if (m_pppTrellis[i][j][s].NextState[1] != ~0)
					cout << "S_" << j << '_' << s << "->" << "S_" << (j + 1) << '_' << m_pppTrellis[i][j][s].NextState[1] << "[style=solid];\n";
			}
		}
		cout << "}\n===========\n";*/
		_ASSERT(_CrtCheckMemory());
	};
	delete[]pCW0;
	delete[]pCW1;
	delete[]pActiveBits;
	delete[]pRowEnd;
	delete[]pRowStart;
	delete[]pExtGM;
}

CTrellisKernelProcessor::~CTrellisKernelProcessor()
{
	for (unsigned i = 0; i < m_pKernel->Size(); i++)
	{
		for (unsigned j = 0; j <= m_pKernel->Size(); j++)
			delete[]m_pppTrellis[i][j];
		delete[]m_pppTrellis[i];
		delete[]m_ppNumOfActiveBits[i];
	};
	delete[]m_pppTrellis;
	delete[]m_ppNumOfActiveBits;
}	



///allocate state variables if necessary
void* CTrellisKernelProcessor::GetState(unsigned Stride  ///enable processing of data blocks of size Stride
					)const
{
	return AlignedAlloc(m_pKernel->Size() * Stride);
};



void CTrellisKernelProcessor::FreeState(void * State)const
{
	tBit* S = (tBit*)State;
	AlignedFree(S);
}
///allocate temporary variables if necessary
void* CTrellisKernelProcessor::GetTemp(unsigned Stride ///enable processing of data blocks of size Stride
)const
{
	return _aligned_malloc(sizeof(MType) * 2  << m_MaxNumOfAB,ALIGN);
}

///deallocate temporary variables
void CTrellisKernelProcessor::FreeTemp(void * Temp)const
{
	_aligned_free(Temp);
}

///copy the state variable
void CTrellisKernelProcessor::CopyState(unsigned phase ///the number of known kernel input symbols
	, unsigned Stride ///the state variable is assumed to be configured for processing of Stride LLRs
	, void* pDest ///destination state variable 
	, const void* pSrc ///source state variable
)const
{
	memcpy(pDest, pSrc, sizeof(tBit)*m_pKernel->Size() * Stride);
};

///compute LLRs for a given kernel phase
void CTrellisKernelProcessor::GetLLRs(unsigned Stride ///number of LLRs to be computed 
	, unsigned phase ///local kernel phase
	, const tBit* pKnownInputSymbols ///known input symbols of the kernel. This is a vector of size Stride*m_pKernel->size
	, const MType* pChannelLLRs ///the LLRs for kernel output symbols. This is a vector of size Stride*m_pKernel->size
	, MType* pLLRs ///output: the LLRs for kernel input symbol phase. This is a vector of size Stride
	, void* pState  ///state variable
	, void* pTemp
)const
{
	KType* pStateMetric0 = (KType*)pTemp;
	KType* pStateMetric1 = pStateMetric0 + (1ull << m_MaxNumOfAB);
	tBit* pOffset = (tBit*)pState;
	if (!phase)
		memset(pOffset, 0, sizeof(tBit)*m_pKernel->Size()*Stride);
	else
	{
			const tBit* pRow = m_pKernel->GetKernel() + (phase-1)*m_pKernel->Size();
			for (unsigned i = 0; i < phase; i++)
				for (unsigned j = 0; j < Stride; j++)
					_ASSERT((pKnownInputSymbols[i*Stride + j] & ~BIT_1) == 0);
			for (unsigned i = 0; i < m_pKernel->Size(); i++)
				for (unsigned j = 0; j < Stride; j++)
					_ASSERT((pOffset[i*Stride + j] & ~BIT_1) == 0);

			MatrixMultiply(1, m_pKernel->Size(), Stride, pKnownInputSymbols + (phase-1)*Stride, pRow, pOffset,false);
	}
	for (unsigned i = 0; i < Stride; i++)
	{
		pStateMetric0[0] = 0;
		for (unsigned j = 0; j < m_pKernel->Size(); j++)
		{
			KType Y = pChannelLLRs[j*Stride + i];
			if (pOffset[j*Stride + i])
				Y = -Y;
			for (unsigned S = 0; S < 1ull << m_ppNumOfActiveBits[phase][j + 1]; S++)
				pStateMetric1[S] = HUGE_VAL;
			unsigned HD = Y < 0;
			for (unsigned long long S = 0; S < 1ull << m_ppNumOfActiveBits[phase][j]; S++)
			{
				for (unsigned z = 0; z < 2; z++)
				{
					unsigned long long S1 = m_pppTrellis[phase][j][S].NextState[z];
					if (S1!= ~0)
					{
						KType NextScore;
						if (z^HD)
						{
							NextScore = pStateMetric0[S] + fabs(Y);
							SUM_COUNT;
						}else  NextScore=pStateMetric0[S];
						CMP_COUNT;
						if (NextScore < pStateMetric1[S1])
							pStateMetric1[S1] = NextScore;
					};
				};
			}
			swap(pStateMetric0, pStateMetric1);
		};
		pLLRs[i] = pStateMetric0[1] - pStateMetric0[0];
	};

};