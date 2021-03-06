// TrellisKernelListDecoder.cpp : Defines the entry point for the console application.
// TrellisKernelListDecoder.cpp : Defines the entry point for the console application.
//
#include <algorithm>
#include "MixedKernelListDecoder.h"
#include "misc.h"

using namespace std;

CMixedKernelListDecoder::CMixedKernelListDecoder(std::istream& Spec, unsigned ListSize):CListKernelEngine(Spec,ListSize)
{
	m_pSortingBuffer = new tPairOfMTypeUnsigned[2 * ListSize];
	m_pActivePaths = new bool[ListSize];
/*#ifdef	ZERO_SKIP
	m_pZeroAdvancePhases = new tPairOfUnsigned[m_UnshortenedLength - m_Dimension + 1];
	unsigned I = 0;
	unsigned F = 0;
	for (unsigned b = 0; b < m_UnshortenedLength ; b++)
	{
		if ((m_pDecisionPoints[b] != ~0) && (m_ppFreezingConstraints[m_pDecisionPoints[b]][0] == b))
		{
			F++;
		}
		else
		{
				unsigned Depth = 0;
				m_pZeroAdvancePhases[I].first = b - F;
				//advance to the next block
				while (F /= m_ppKernels[m_NumOfLayers - Depth - 1]->Size())
				{
					Depth++;
				};
				if (Depth > 0)
				{
					m_pZeroAdvancePhases[I].second = Depth;
					unsigned n = 1;
					for (unsigned i = 0; i < Depth; i++)
						n *= m_ppKernels[m_NumOfLayers - 1 - i]->Size();
					b = m_pZeroAdvancePhases[I].first + n - 1;;
				//	cout << "Skipping over " << m_pZeroAdvancePhases[I].second << " layers from phase " << m_pZeroAdvancePhases[I].first << endl;
					I++;
				};
		}
	}
	m_pZeroAdvancePhases[I].first = ~0;//boundary condition;
#endif
*/
	_ASSERT(_CrtCheckMemory());

}

CMixedKernelListDecoder::~CMixedKernelListDecoder()
{
	delete[]m_pSortingBuffer;
	delete[]m_pActivePaths;
#ifdef	ZERO_SKIP
	delete[]m_pZeroAdvancePhases;
#endif
}

void CMixedKernelListDecoder::ContinuePathsFrozen(unsigned phi)
{
#ifdef ENABLE_LOGGING
	cout << phi << ':';
#endif
	for (unsigned l = 0; l < m_MaxNumOfPaths; l++)
	{
		if (m_pActivePaths[l])
		{
			//_ASSERT(_CrtCheckMemory());
			const MType* pLLR = IterativelyCalcS(l, phi);
#ifdef ENABLE_LOGGING
			cout << *pLLR << ' ';
#endif
			tBit C=BIT_0;
			if (m_pDFValueBits[phi]!=~0) 
				C= ((m_pDFMask[l] >> m_pDFValueBits[phi]) & 1) ? BIT_1 : BIT_0;
			//_ASSERT(_CrtCheckMemory());

			if ((C > 0) ^ (*pLLR < 0))
			{
				m_pR[l] -= fabs(*pLLR);
				SUM_COUNT;
			};
			tBit* pC = GetArrayPointer_C(m_NumOfLayers, l, phi);
			//_ASSERT(_CrtCheckMemory());

			pC[phi%m_ppKernels[m_NumOfLayers - 1]->Size()] = C;
			if (C)
				m_pDFMask[l] ^= m_pDFCorrection[phi];
			IterativelyUpdateC(l, phi);
			//_ASSERT(_CrtCheckMemory());
		}
	};
#ifdef ENABLE_LOGGING
	cout << endl;
#endif
}

void CMixedKernelListDecoder::ContinuePathsUnfrozen(unsigned phi)
{
	unsigned J = 0;
#ifdef ENABLE_LOGGING
	cout << phi << ':';
#endif
	for (unsigned l = 0; l < m_MaxNumOfPaths; l++)
	{
		if (m_pActivePaths[l])
		{
			const MType* pLLR = IterativelyCalcS(l, phi);
#ifdef ENABLE_LOGGING
			cout << *pLLR << ' ';
#endif
			bool D = *pLLR < 0;
			m_pSortingBuffer[ J].first = m_pR[l];
			m_pSortingBuffer[J].second = 2 * l + D;
			m_pSortingBuffer[J + 1].first = m_pR[l] - fabs(*pLLR); SUM_COUNT;
			m_pSortingBuffer[J+1].second = 2 * l + (D^1);
			J+=2;
		}
	};
#ifdef ENABLE_LOGGING
	cout << endl;
#endif
	sort(m_pSortingBuffer, m_pSortingBuffer + J, CPairComparator());
	//sort(m_pSortingBuffer, m_pSortingBuffer + J);
	unsigned __int8* pContinuation = (unsigned __int8*)alloca(m_MaxNumOfPaths);
	memset(pContinuation, 0, sizeof(unsigned __int8)*m_MaxNumOfPaths);
	for (unsigned i = 0; i < min(J,m_MaxNumOfPaths); i++)
	{
		pContinuation[m_pSortingBuffer[i].second >> 1] |= 1u << (m_pSortingBuffer[i].second & 1);
	};
	for (unsigned i = 0; i < m_MaxNumOfPaths; i++)
	{
		if (m_pActivePaths[i] && !pContinuation[i])
		{
			KillPath(i);
			m_pActivePaths[i] = false;
		};
	};

	for (unsigned l = 0; l < m_MaxNumOfPaths; l++)
	{
		switch (pContinuation[l])
		{
		case 0: break;
		case 1://extend it with 0. Path score is not going to change
		{
			tBit* pC = GetArrayPointer_C(m_NumOfLayers, l, phi);
			pC[phi%m_ppKernels[m_NumOfLayers-1]->Size()] = BIT_0;
			break;
		};
		case 2://extend it with 1. Path score is not going to change
		{
			tBit* pC = GetArrayPointer_C(m_NumOfLayers, l, phi);
			pC[phi%m_ppKernels[m_NumOfLayers - 1]->Size()] = BIT_1;
			m_pDFMask[l] ^= m_pDFCorrection[phi];
			break;
		};
		case 3:
		{
			const MType* pLLR = CGetArrayPointer_S(m_NumOfLayers,l);
			tBit* pC = GetArrayPointer_C(m_NumOfLayers, l, phi);
			tBit C=pC[phi%m_ppKernels[m_NumOfLayers - 1]->Size()]=(*pLLR<0)?BIT_1:BIT_0;
			unsigned l1 = ClonePath(l);
			tBit* pC1 = GetArrayPointer_C(m_NumOfLayers, l1, phi);
			pC1[phi%m_ppKernels[m_NumOfLayers - 1]->Size()] = C^BIT_1;
			m_pActivePaths[l1] = true;
			m_pR[l1] = m_pR[l] - fabs(*pLLR); SUM_COUNT;
			m_pDFMask[l1] = m_pDFMask[l];
			if (C)
				m_pDFMask[l] ^= m_pDFCorrection[phi];
			else 
				m_pDFMask[l1] ^= m_pDFCorrection[phi];
			break;
		}
		}
	}
	for (unsigned l = 0; l < m_MaxNumOfPaths; l++)
	{
		if (m_pActivePaths[l])
			IterativelyUpdateC(l, phi);
	}

}
/*
#ifdef ZERO_SKIP
///process a block of static frozen symbols of a given depth. Return the size of the skipped block
unsigned CMixedKernelListDecoder::SkipZeroBlock(unsigned phi, unsigned Depth)
{
	unsigned BlockSize = 1;
	for (unsigned i = 0; i < Depth; i++)
		BlockSize *= m_ppKernels[m_NumOfLayers - 1-i]->Size();
	for (unsigned l = 0; l < m_MaxNumOfPaths; l++)
	{
		if (m_pActivePaths[l])
		{
			const MType* pLLRs = IterativelyCalcS(l, phi, Depth);
			for (unsigned i = 0; i < BlockSize; i++)
				m_pR[l] -= (pLLRs[i]<0)?fabs(pLLRs[i]):0;
			CMP_BLOCK_COUNT(BlockSize);
			IterativelyUpdateC(l, phi+BlockSize-1, Depth);
		}
	};
	return BlockSize;
}
#endif
*/
//soft-input decoder. Produces  a list of information vectors
//returns the number of vectors obtained, or -1 in case of error
int CMixedKernelListDecoder::Decode(const MType* pLLR,//log(P{c_i=1|y_i}/P{c_i=0|y_i})
	tBit* pInfVectorList, //output array. Must be sufficiently big to accomodate GetMaxListSize() vectors of length GetDimension()
	tBit* pCodewordList  //output array. If non-zero, must be sufficiently big to accomodate GetMaxListSize() vectors of length GetLength()
)
{
	_ASSERT(_CrtCheckMemory());

	Cleanup();
	memset(m_pActivePaths, 0, sizeof(bool)*m_MaxNumOfPaths);
	unsigned PID = AssignInitialPath();

	m_pActivePaths[PID] = true;
	m_pDFMask[PID] = 0;
	m_pR[PID] = 0;

	MType* pS = GetArrayPointer_S(0, PID, 0);

	LoadLLRs(pS, pLLR);
	unsigned Z = 0;
	for (unsigned phi = 0; phi  < m_UnshortenedLength; phi++)
	{
/*#ifdef ZERO_SKIP
		if (phi == m_pZeroAdvancePhases[Z].first)
		{
			phi += SkipZeroBlock(phi, m_pZeroAdvancePhases[Z].second)-1;
#ifdef ENABLE_LOGGING
			cout <<"\n!"<< phi << ' ' << m_pR[0] << ' ' << m_pR[1] << endl;
#endif
			Z++;
			continue;
		}
#endif*/
		if (m_pDecisionPoints[phi] != ~0)
			ContinuePathsFrozen(phi);
		else
			ContinuePathsUnfrozen(phi);
#ifdef ENABLE_LOGGING
		cout<<phi<<' '<<m_pR[0]<<' '<<m_pR[1]<<endl;
#endif
		//_ASSERT(_CrtCheckMemory());
	//	_ASSERT(m_pR[PID] == 0);
	};
	unsigned J = 0;
	for (unsigned l = 0; l < m_MaxNumOfPaths; l++)
		if (m_pActivePaths[l])
			m_pSortingBuffer[J++] = tPairOfMTypeUnsigned(m_pR[l], l);
	sort(m_pSortingBuffer, m_pSortingBuffer + J, CPairComparator());

	for (unsigned l = 0; l < J; l++)
	{
		const tBit* pCW = CExtractCodeword(m_pSortingBuffer[l].second);
		///update this for shortening and puncturing
		if (pCodewordList)
			Shorten(pCW,pCodewordList + l*m_Length);
		ExtractInformationBitsUnshortened(pCW, pInfVectorList + l*m_Dimension);
	}
	_ASSERT(_CrtCheckMemory());
	return J;
}

