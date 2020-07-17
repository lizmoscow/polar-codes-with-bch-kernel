#include <set>
#include "LinAlg.h"
#include "KernelListEngine.h"

using namespace std;
CListKernelEngine::CListKernelEngine(std::istream& SpecFile, unsigned MaxNumOfPaths):CMixedKernelEncoder(SpecFile),CTVMemoryEngine(MaxNumOfPaths)
{
	m_pDFCorrection = new tDFMask[m_UnshortenedLength];
	set<unsigned> AvailableBits;
	for (unsigned i = 0; i < 8 * sizeof(tDFMask); i++)
	{
		AvailableBits.insert(i);
	};
	m_pDFValueBits = new unsigned[m_UnshortenedLength];
	memset(m_pDFCorrection, 0, sizeof(tDFMask)*m_UnshortenedLength);
	memset(m_pDFValueBits, ~0, sizeof(unsigned)*m_UnshortenedLength);
	for (unsigned i = 0; i<m_UnshortenedLength; i++)
	{
		if (m_pDFValueBits[i] < 8 * sizeof(tDFMask))
			//this bit is no longer needed
			AvailableBits.insert(m_pDFValueBits[i]);
		//check if this symbol is a starting one in any of the following dynamic freezing constraints
		for (unsigned j = i + 1; j < m_UnshortenedLength; j++)
		{
			if ((m_pDecisionPoints[j] != ~0) && m_ppFreezingConstraints[m_pDecisionPoints[j]] && (*(m_ppFreezingConstraints[m_pDecisionPoints[j]]) == i))
			{
				if (AvailableBits.size() == 0)
					throw Exception("Too many dynamic freezing constraints are simultaneously active");
				unsigned B = *AvailableBits.begin();
				m_pDFValueBits[j] = B;
				AvailableBits.erase(B);
				unsigned* pFC = m_ppFreezingConstraints[m_pDecisionPoints[j]];
				do
				{
					m_pDFCorrection[*pFC] ^= 1ull << B;
				} while (*(++pFC) != j);
				//cleanup the bit in order to allow its reuse
				m_pDFCorrection[j] ^= 1ull << B;
			};
		};
	};


	m_pppLLREngines = new const  CKernProcLLR**[m_NumOfLayers];
	m_pOuterLengths = new unsigned[m_NumOfLayers+1];
	m_pMaxEngineID = new unsigned[m_NumOfLayers];
	memset(m_pMaxEngineID, 0, sizeof(unsigned)*m_NumOfLayers);
#ifdef ENABLE_SIMPLIFIED_ENGINES
	eFreezingState* pFS = new eFreezingState[m_UnshortenedLength];
	for (unsigned i = 0; i < m_UnshortenedLength; i++)
	{
		if ((m_pDecisionPoints[i] == ~0) && !m_pDFCorrection[i])
		{
			pFS[i] = efsFullyUnfrozen;
			continue;
		};
		if ((m_pDecisionPoints[i] != ~0) && (m_ppFreezingConstraints[m_pDecisionPoints[i]][0] == i) && !m_pDFCorrection[i])
			//only static frozen symbols can be used for simplification
			pFS[i] = efsFullyFrozen;
		else
		{
			//the whole block should be excluded from consideration
			unsigned i0 = i - (i%m_ppKernels[m_NumOfLayers - 1]->Size());
			for (unsigned j = 0; j<m_ppKernels[m_NumOfLayers - 1]->Size(); j++)
				pFS[i0 + j] = efsConstrained;
			i = i0 + m_ppKernels[m_NumOfLayers - 1]->Size() - 1;
		};
	};
#endif
	_ASSERT(_CrtCheckMemory());

	unsigned N = m_UnshortenedLength;
	unsigned Stride = 1;
	m_ppProcessorID = new unsigned*[m_NumOfLayers];
	m_pppEngineTemp = new void**[m_NumOfLayers];
	for (int i = m_NumOfLayers - 1; i >= 0; i--)
	{
		N /= m_ppKernels[i]->Size();
		m_ppProcessorID[i] = new unsigned[N];
		_ASSERT(_CrtCheckMemory());

#ifdef ENABLE_SIMPLIFIED_ENGINES
		eFreezingState* pCurFS = pFS;
		for (unsigned j = 0; j < N; j++)
		{
			m_ppProcessorID[i][j] = m_ppKernels[i]->GetProcessorID(pCurFS);
			if (m_ppProcessorID[i][j])
			{
					cout << "Enabled simplified processing ("<< m_ppProcessorID[i][j]<<") for block " << j << " at layer " << i << endl;
					m_pMaxEngineID[i] = max(m_pMaxEngineID[i], m_ppProcessorID[i][j]);
			};
			unsigned F = 0;
			unsigned U = 0;
			bool Constrained = false;
			for (unsigned t = 0; t < m_ppKernels[i]->Size(); t++)
			{
				switch (pCurFS[t])
				{
				case efsFullyFrozen:F++; break;
				case efsFullyUnfrozen:U++; break;
				case efsConstrained:Constrained = true; break;
				};
			};
			if (Constrained)
				pFS[j] = efsConstrained;
			else
				if (F == m_ppKernels[i]->Size())
					pFS[j] = efsFullyFrozen;
				else
					if (U == m_ppKernels[i]->Size())
						pFS[j] = efsFullyUnfrozen;
					else pFS[j] = efsPartiallyFrozen;
					pCurFS += m_ppKernels[i]->Size();
		};
#else
		memset(m_ppProcessorID[i], 0, sizeof(unsigned)*N );

#endif
		_ASSERT(_CrtCheckMemory());

		m_pppLLREngines[i] = new const CKernProcLLR*[m_pMaxEngineID[i] + 1];
		m_pppEngineTemp[i] = new void*[m_pMaxEngineID[i] + 1];
		for (unsigned e = 0; e <= m_pMaxEngineID[i]; e++)
		{
			m_pppLLREngines[i][e] = KernelFactory.GetKernelProcessor(m_ppKernels[i], e);
			m_pppEngineTemp[i][e] = m_pppLLREngines[i][e]->GetTemp(Stride);
		};
		Stride *= m_ppKernels[i]->Size();
	};
#ifdef ENABLE_SIMPLIFIED_ENGINES
	delete[]pFS;
#endif
	_ASSERT(_CrtCheckMemory());

	unsigned MaxNumOfEngines = 1;
	m_pOuterLengths[0] = m_UnshortenedLength;
	for (unsigned i = 0; i < m_NumOfLayers; i++)
	{
		MaxNumOfEngines = max(MaxNumOfEngines, m_pMaxEngineID[i] + 1);
		m_pOuterLengths[i + 1] = m_pOuterLengths[i]/ m_ppKernels[i]->Size();
	};

	m_pppEngineState = new void**[MaxNumOfEngines];
	for (unsigned e = 0; e < MaxNumOfEngines; e++)
	{
		m_pppEngineState[e] = new void*[(m_NumOfLayers + 1)*m_MaxNumOfPaths];
		memset(m_pppEngineState[e], 0, sizeof(void*)*(m_NumOfLayers + 1)*m_MaxNumOfPaths);
	};
	Init(m_NumOfLayers, m_NumOfLayers, 0);
	m_pR=new KType[m_MaxNumOfPaths];
	m_ppArrayPointer_C = new tBit*[(m_NumOfLayers + 1)*m_MaxNumOfPaths];
	m_ppArrayPointer_S=new MType*[(m_NumOfLayers + 1)*m_MaxNumOfPaths];
	memset(m_ppArrayPointer_C, 0, sizeof(tBit*)*(m_NumOfLayers + 1)*m_MaxNumOfPaths);
	memset(m_ppArrayPointer_S, 0, sizeof(MType*)*(m_NumOfLayers + 1)*m_MaxNumOfPaths);
	_ASSERT(_CrtCheckMemory());
#ifdef  CALL_ACCOUNTING
	for(unsigned i=0;i<m_NumOfLayers;i++)
		for (unsigned e = 0; e <= m_pMaxEngineID[i]; e++)
		{
			if (m_CallAccountTable.find(m_pppLLREngines[i][e]) == m_CallAccountTable.end())
			{
				unsigned* pC = m_CallAccountTable[m_pppLLREngines[i][e]] = new unsigned[m_ppKernels[i]->Size()];
				memset(pC, 0, sizeof(unsigned)*m_ppKernels[i]->Size());
			}
		}

	m_NumOfDecoderCalls = 0;
#endif

}

CListKernelEngine::~CListKernelEngine()
{
	unsigned E = 0;
	for (unsigned i = 0; i <= m_NumOfLayers; i++)
	{
		for (unsigned s = 0; s < m_MaxNumOfPaths; s++)
		{
			m_Memory.Free(m_ppArrayPointer_C[s*(m_NumOfLayers + 1) + i]);
			m_Memory.Free(m_ppArrayPointer_S[s*(m_NumOfLayers + 1) + i]);
			if (i < m_NumOfLayers)
			{
				for(unsigned e=0;e<=m_pMaxEngineID[i];e++)
					if (m_pppLLREngines[i][e])
						m_pppLLREngines[i][e]->FreeState(m_pppEngineState[e][s*(m_NumOfLayers + 1) + i]);
			};
		};
		if (i < m_NumOfLayers)
		{
			E = max(E, m_pMaxEngineID[i]);
			for (unsigned e = 0; e <= m_pMaxEngineID[i]; e++)
				if (m_pppLLREngines[i][e])
					m_pppLLREngines[i][e]->FreeTemp(m_pppEngineTemp[i][e]);
			delete[]m_ppProcessorID[i];
			delete[]m_pppEngineTemp[i];
			delete[]m_pppLLREngines[i];
		}
	};
	E++;
	for (unsigned e = 0; e < E; e++)
	{
		delete[]m_pppEngineState[e];
	};
	delete[]m_ppArrayPointer_C;
	delete[]m_ppArrayPointer_S;
	delete[]m_pppLLREngines;
	delete[]m_pOuterLengths;
	delete[]m_pppEngineTemp;
	delete[]m_pppEngineState;
	delete[]m_pDFCorrection;
	delete[]m_pDFValueBits;
	delete[]m_pR;
/*#ifdef ENABLE_SIMPLIFIED_ENGINES
	delete[]m_ppSimplifiedState;
	delete[]m_ppSimplificationPossible;
	delete[]m_ppSimplifiedLLREngines;
#endif
*/
#ifdef  CALL_ACCOUNTING
	tAccountMap::const_iterator it = m_CallAccountTable.begin();
	while (it != m_CallAccountTable.end())
	{
		cout << it->first->GetName() << ':';
		unsigned l = it->first->Size();
		for (unsigned i = 0; i < l; i++)
		{
			cout << double(it->second[i])/m_NumOfDecoderCalls<<' ';
		};
		cout << endl;
		it++;
	}
#endif
}

///allocate a new array at layer lambda
void CListKernelEngine::AllocateEntry(unsigned lambda, unsigned s1)
{
#ifndef USE_MEMORY_MANAGER
	if (m_ppArrayPointer_S[s1] != nullptr)
		return;
#endif

	//unsigned Offset = lambda*(m_NumOfLayers + 1) + s1;
	unsigned CSize = m_pOuterLengths[lambda];
	if (lambda > 0)
		CSize *= m_ppKernels[lambda - 1]->Size();
	m_ppArrayPointer_C[s1] = m_Memory.Allocate<tBit>(CSize);

	m_ppArrayPointer_S[s1] = m_Memory.Allocate<MType>(m_pOuterLengths[lambda]);

	if (lambda >0)
	{
		for (unsigned e = 0; e <= m_pMaxEngineID[lambda - 1]; e++)
		{
			m_pppEngineState[e][s1] = m_pppLLREngines[lambda - 1][e]->GetState(m_pOuterLengths[lambda]);
/*#ifdef ENABLE_SIMPLIFIED_ENGINES
			if (m_ppSimplifiedLLREngines[lambda - 1])
				m_ppSimplifiedState[s1] = m_ppSimplifiedLLREngines[lambda - 1]->GetState(m_pOuterLengths[lambda]);
#endif*/
		};
	};

}


void CListKernelEngine::IterativelyUpdateC(
	unsigned l      //path ID
	,unsigned phi     //phase
	,unsigned Layers2Skip 
)
{
	unsigned lambda = m_NumOfLayers;
	unsigned Stride = 1;
	while ((lambda>0)&&!((phi + 1) % m_ppKernels[lambda - 1]->Size()))
	{
		unsigned Psi = phi / m_ppKernels[lambda - 1]->Size();
		unsigned NextStride = Stride*m_ppKernels[lambda - 1]->Size();
		if (!Layers2Skip)
		{
			unsigned Phi0 = (lambda > 1) ? (Psi%m_ppKernels[lambda - 2]->Size())*NextStride : 0;
			const tBit* pSrc = CGetArrayPointer_C(lambda, l);
			tBit* pDest = GetArrayPointer_C(lambda - 1, l, phi);
			//		_ASSERT(_CrtCheckMemory());
			const CKernProcLLR* pEngine=m_pppLLREngines[lambda-1][m_ppProcessorID[lambda-1][Psi]];
/*#ifdef ENABLE_SIMPLIFIED_ENGINES
			pEngine = (m_ppSimplificationPossible[lambda - 1]&&m_ppSimplificationPossible[lambda-1][Psi]) ? m_ppSimplifiedLLREngines[lambda-1] : m_ppLLREngines[lambda-1];
#else 
			pEngine = m_ppLLREngines[lambda-1];
#endif*/
			pEngine->GetKernelOutput(Stride, pSrc, pDest + Phi0);
#ifdef ENABLE_LOGGING
			cout << "C" << lambda << ':';
			for (unsigned i = 0; i < NextStride; i++)
				cout << pDest[Phi0 + i];
			cout << endl;
#endif
		}
		else
		{
			Layers2Skip--;
			if (!Layers2Skip)
			{
				//it is time to zero the output
				tBit* pDest = GetArrayPointer_C(lambda - 1, l, phi);
				unsigned Phi0 = (lambda > 1) ? (Psi%m_ppKernels[lambda - 2]->Size())*NextStride : 0;
				memset(pDest + Phi0, 0, sizeof(tBit)*NextStride);
			}
		};
		Stride = NextStride;
		phi = Psi;
		lambda--;
	//	_ASSERT(_CrtCheckMemory());

	};
}

///copy the data if the array is cloned
unsigned CListKernelEngine::UpdateIndex(unsigned lambda ///layer ID
	, unsigned *pIndexArray ///index array for the path
	, unsigned CurPathPhase ///current path phase
)
{
	unsigned OldIndex = pIndexArray[lambda];

	unsigned NewIndex = CTVMemoryEngine::UpdateIndex(lambda, pIndexArray);
	if ((OldIndex!=m_EmptyIndex)&&(OldIndex != NewIndex))
	{
		//copy C variables
		const tBit* pOldC = m_ppArrayPointer_C[OldIndex];
		tBit* pNewC= m_ppArrayPointer_C[NewIndex];
		unsigned N = (lambda > 0) ? (m_pOuterLengths[lambda - 1]) : m_pOuterLengths[lambda];
		memcpy(pNewC, pOldC, sizeof(tBit)*N);///check this

		//copy engine state
		if (lambda > 0)
		{
			unsigned f = CurPathPhase;
			for (int t = m_NumOfLayers; t >= int(lambda); t--)
				f /= m_ppKernels[t - 1]->Size();
			unsigned e = m_ppProcessorID[lambda - 1][f];

			const void* pOldState = m_pppEngineState[e][OldIndex];
			void* pNewState = m_pppEngineState[e][NewIndex];
			//copy engine state
			const CKernProcLLR* pEngine=m_pppLLREngines[lambda-1][e];
/*#ifdef ENABLE_SIMPLIFIED_ENGINES
			if (m_ppSimplificationPossible[lambda - 1]&&m_ppSimplificationPossible[lambda - 1][f])
			{
				pEngine = m_ppSimplifiedLLREngines[lambda - 1];
				pNewState = m_ppSimplifiedState[NewIndex];
				pOldState = m_ppSimplifiedState[OldIndex];
			}
			else
			{
				pEngine=m_ppLLREngines[lambda - 1];
			}

#else 
			pEngine = m_ppLLREngines[lambda-1];
#endif*/
			pEngine->CopyState(CurPathPhase%m_ppKernels[lambda-1]->Size(), m_pOuterLengths[lambda], pNewState, pOldState);
		}


	}
	return NewIndex;
}

///compute the LLR for the input symbol of a polarizing transformation
const MType* CListKernelEngine::IterativelyCalcS(
	unsigned l      ///path ID
	,unsigned phi    ///phase
	,unsigned Layers2Skip ///number of last layers to be skipped
)
{
	unsigned m = m_NumOfLayers-1;
	unsigned phi0 = phi;
	unsigned* pB = (unsigned*)alloca(sizeof(unsigned)*m_NumOfLayers);

	pB[m] = phi/ m_ppKernels[m]->Size();
	while (m && !(phi%m_ppKernels[m ]->Size()))
	{
		phi /= m_ppKernels[m ]->Size();
		m--;
		pB[m] = pB[m+1]/ m_ppKernels[m]->Size();
	};
	unsigned* pIdx = GetIndexArrayPointer(l);
	MType* pDest=0;
	const MType* pSrc = CGetArrayPointer_S(m, pIdx);
#ifdef ENABLE_LOGGING
	cout << phi << '#';
#endif

	//phi = phi0;
	for (unsigned j = m; j < m_NumOfLayers-Layers2Skip; j++)
	{
		pDest = GetArrayPointer_S(j + 1, pIdx,phi);
		const tBit* pC = CGetArrayPointer_C(j + 1, pIdx);
		unsigned e = m_ppProcessorID[j][pB[j]];
		const CKernProcLLR* pEngine=m_pppLLREngines[j][e];
		void* pKernelState=m_pppEngineState[e][pIdx[j+1]];
/*
#ifdef ENABLE_SIMPLIFIED_ENGINES
		if (m_ppSimplificationPossible[j]&&m_ppSimplificationPossible[j][pB[j]])
		{
			pEngine = m_ppSimplifiedLLREngines[j];
			pKernelState = m_ppSimplifiedState[pIdx[j + 1]];
		}
		else
		{
			pEngine = m_ppLLREngines[j];
			pKernelState = m_ppEngineState[pIdx[j + 1]];
		}

#else 
		pEngine = m_ppLLREngines[j];
		pKernelState = m_ppEngineState[pIdx[j + 1]];

#endif*/
		unsigned LocalPhase = phi%m_ppKernels[j]->Size();
#ifdef ENABLE_LOGGING
		cout << "C:";
		for (unsigned i = 0; i < LocalPhase*m_pOuterLengths[j + 1]; i++)
			cout << pC[i];
		cout << "\nS:";
		for (unsigned i = 0; i < m_ppKernels[j]->Size()*m_pOuterLengths[j + 1]; i++)
			cout << pSrc[i] << ' ';
		cout << endl;
#endif
		pEngine->GetLLRs(m_pOuterLengths[j+1], LocalPhase, pC, pSrc, pDest, pKernelState, m_pppEngineTemp[j][e]);
#ifdef  CALL_ACCOUNTING
		unsigned* pAcct = m_CallAccountTable[pEngine];
		pAcct[LocalPhase]++;
#endif
#ifdef ENABLE_LOGGING
		for (unsigned t = 0; t < m_pOuterLengths[j + 1]; t++)
			cout << pDest[t] << ' ';
		cout << '\n';
#endif
		phi =0;
		pSrc = pDest;
		//_ASSERT(_CrtCheckMemory());

	}
	_ASSERT(*pDest >= 0 || *pDest <= 0);
	return pDest;
}
