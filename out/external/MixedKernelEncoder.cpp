#include <string>
#include "MixedKernelEncoder.h"
#include "LinAlg.h"

using namespace std;

CMixedKernelEncoder::CMixedKernelEncoder(std::istream& SpecFile)
{
	unsigned NumOfShortened, NumOfPunctured;
	SpecFile >> m_Length >> m_Dimension >> m_MinDist >> m_NumOfLayers>> NumOfShortened >> NumOfPunctured;
	if (!SpecFile)
		throw Exception("Error reading file header");
	if (m_Dimension > m_Length)
		throw Exception("Code dimension cannot exceed code length");

	m_ppKernels = new const CBinaryKernel*[m_NumOfLayers];
	//read the names of the kernels involved
	m_UnshortenedLength = 1;
	for (unsigned i = 0; i < m_NumOfLayers; i++)
	{
		string KernelName;
		SpecFile >> KernelName;
		m_ppKernels[i] = KernelFactory.GetKernelByName(KernelName);
		m_UnshortenedLength *= m_ppKernels[i]->Size();
	};
	if (m_Length + NumOfShortened + NumOfPunctured != m_UnshortenedLength)
		throw Exception("Code length mismatch");
	m_pSymbolType = 0;
	if (NumOfPunctured + NumOfShortened)
	{
		m_pSymbolType = new SymbolType[m_UnshortenedLength];
		memset(m_pSymbolType, 0, sizeof(SymbolType)*m_UnshortenedLength);
		for (unsigned i = 0; i < NumOfShortened; i++)
		{
			unsigned S;
			SpecFile >> S;
			if (!SpecFile)
				throw Exception("Error loading shortened symbols");
			if (S >= m_UnshortenedLength)
				throw Exception("Invalid shortened symbol %d", S);
			m_pSymbolType[S] = stShortened;
		};
		for (unsigned i = 0; i < NumOfPunctured; i++)
		{
			unsigned P;
			SpecFile >> P;
			if (!SpecFile)
				throw Exception("Error loading punctured symbols");
			if (P >= m_UnshortenedLength)
				throw Exception("Invalid punctured symbol %d", P);
			m_pSymbolType[P] = stPunctured;
		};

	}
	/*if (m_NumOfPunctured + m_NumOfShortened)
		throw Exception("Puncturing and shortening are not supported yet");*/
	//load freezing constraints
	m_ppFreezingConstraints = new unsigned*[m_UnshortenedLength - m_Dimension];
	m_pDecisionPoints = new unsigned[m_UnshortenedLength];
	memset(m_pDecisionPoints, ~0, sizeof(unsigned)*m_UnshortenedLength);
	for (unsigned i = 0; i < m_UnshortenedLength - m_Dimension; i++)
	{
		unsigned w;
		SpecFile >> w;
		if (!SpecFile)
			throw Exception("Error reading freezing constraint %d", i);

		m_ppFreezingConstraints[i] = new unsigned[w];
		unsigned J = 0;
		for (unsigned j = 0; j < w; j++)
		{
			SpecFile >> m_ppFreezingConstraints[i][J];
			if (!SpecFile)
				throw Exception("Error parsing freezing constraint %d", i);
			if (m_ppFreezingConstraints[i][J] >= m_UnshortenedLength)
				throw Exception("Invalid term %d in freezing constraint %d", m_ppFreezingConstraints[i][J], i);
			if ((J>0)&&(m_ppFreezingConstraints[i][J]<= m_ppFreezingConstraints[i][J-1]))
				throw Exception("Invalid freezing constraint %d", i);
			//check if this symbol is a non-trivial one
			unsigned DP = m_pDecisionPoints[m_ppFreezingConstraints[i][J]];
			if (DP != ~0)
			{
				/*if (m_ppFreezingConstraints[DP][0] == m_ppFreezingConstraints[i][J])
					//this symbol is in fact static frozen
					continue;*/
			};
			J++;
		};
		if (m_pDecisionPoints[m_ppFreezingConstraints[i][J - 1]] != ~0)
			throw Exception("Duplicate freezing constraint on symbol %d", m_ppFreezingConstraints[i][J - 1]);
		else
			m_pDecisionPoints[m_ppFreezingConstraints[i][J - 1]] = i;
	}
	m_pEncodingBuffer = AlignedAlloc(m_UnshortenedLength);
	m_pEncodingBuffer2 = AlignedAlloc(m_UnshortenedLength);
	m_pEncodedInformation = AlignedAlloc(m_UnshortenedLength);
}


CMixedKernelEncoder::~CMixedKernelEncoder()
{
	AlignedFree(m_pEncodingBuffer);
	AlignedFree(m_pEncodingBuffer2);
	AlignedFree(m_pEncodedInformation);
	delete[]m_pDecisionPoints;
	for (unsigned i = 0; i < m_UnshortenedLength - m_Dimension; i++)
		delete[]m_ppFreezingConstraints[i];
	delete[]m_ppFreezingConstraints;
	delete[]m_ppKernels;
	delete[]m_pSymbolType;

}

///shorten and puncture the original polar codeword
void CMixedKernelEncoder::Shorten(const tBit* pPolarCW ///the original polar codeword
	, tBit* pShortened ///shortened and punctured codeword
)
{
	if (m_pSymbolType)
	{
		unsigned I = 0;
		for (unsigned i = 0; i < m_UnshortenedLength; i++)
		{
			switch (m_pSymbolType[i])
			{
			case stNormal: pShortened[I++] = pPolarCW[i];
				break;
			case stShortened:
				if (pPolarCW[i])
					throw Exception("Invalid shortening specification at symbol %d", i);
				break;
			case stPunctured:
				break;
			default: _ASSERT(false);
			}
		}
	}
	else memcpy(pShortened, pPolarCW, sizeof(tBit)*m_UnshortenedLength);
}

///encode a block of data into a codeword
void CMixedKernelEncoder::Encode(const tBit* pSrc ///the information symbols to be encoded
	, tBit* pEncoded ///the codeword
)
{
	_ASSERT(_CrtCheckMemory());
	//evaluate frozen symbols
	for (unsigned i = 0; i < m_UnshortenedLength; i++)
	{
		if (m_pDecisionPoints[i] != ~0)
		{
			const unsigned* pFC = m_ppFreezingConstraints[m_pDecisionPoints[i]];
			tBit C = 0;
			while (*pFC != i)
				C ^= m_pEncodingBuffer[*(pFC++)];
			m_pEncodingBuffer[i] = C;
		}
		else m_pEncodingBuffer[i] = *(pSrc++);
	};
	memcpy(m_pEncodedInformation, m_pEncodingBuffer, sizeof(tBit)*m_UnshortenedLength);
	unsigned Stride = 1;
	for (int L = m_NumOfLayers - 1; L >= 0; L--)
	{
		unsigned NumOfBlocks = m_UnshortenedLength / (Stride*m_ppKernels[L]->Size());
		unsigned NextStride = Stride*m_ppKernels[L]->Size();
		for (unsigned i = 0; i < NumOfBlocks; i++)
		{
			m_ppKernels[L]->Multiply(Stride, m_pEncodingBuffer + i*NextStride, m_pEncodingBuffer2 + i*NextStride);
		}
		swap(m_pEncodingBuffer, m_pEncodingBuffer2);
		Stride = NextStride;
		
	}
	Shorten(m_pEncodingBuffer, pEncoded);
	_ASSERT(_CrtCheckMemory());
}

///unpuncture and unshorten the code
void CMixedKernelEncoder::LoadLLRs(MType* pDest ///the LLR vector as required by the decoder
	, const MType* pSrc
)
{
	if (!m_pSymbolType)
		memcpy(pDest, pSrc, sizeof(MType)*m_Length);
	else
	{
		unsigned I = 0;
		for (unsigned i = 0; i < m_UnshortenedLength; i++)
		{
			switch (m_pSymbolType[i])
			{
			case stNormal: pDest[i]=pSrc[I++];
				break;
			case stShortened:
				pDest[i] = MTYPE_UPPER_BOUND;
				break;
			case stPunctured:
				pDest[i] = 0;
				break;
			}
		};
	}
}

void CMixedKernelEncoder::ExtractInformationBits(const tBit* pCodeword ///the codeword
	, tBit* pInfBits ///information symbols
)
{
	throw Exception("Not implemented yet");
}

///extract information bits from a codeword
void CMixedKernelEncoder::ExtractInformationBitsUnshortened(const tBit* pCodeword ///the codeword
	, tBit* pInfBits ///information symbols
)
{
	memcpy(m_pEncodingBuffer, pCodeword, sizeof(tBit)*m_UnshortenedLength);
	unsigned NumOfBlocks = 1;
	unsigned PrevStride = m_UnshortenedLength;
	for (unsigned L = 0; L < m_NumOfLayers; L++)
	{
		unsigned Stride = PrevStride/ m_ppKernels[L]->Size();
		for (unsigned i = 0; i < NumOfBlocks; i++)
		{
			m_ppKernels[L]->InverseMultiply(Stride, m_pEncodingBuffer + i*PrevStride, m_pEncodingBuffer2 + i*PrevStride);
		}
		PrevStride = Stride;
		NumOfBlocks *= m_ppKernels[L]->Size();
		swap(m_pEncodingBuffer, m_pEncodingBuffer2);
	};
	for (unsigned i = 0; i < m_UnshortenedLength; i++)
	{
		if (m_pDecisionPoints[i] == ~0)
			*(pInfBits++) = m_pEncodingBuffer[i];
	}

}

