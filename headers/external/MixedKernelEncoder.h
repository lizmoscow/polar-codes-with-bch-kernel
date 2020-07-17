#pragma once
#include <fstream>
#include "Kernel.h"
#include "Codec.h"
#include "DynFrozenX.h"

class CMixedKernelEncoder :public CBinaryEncoder
{
protected:
	///minimum distance of the code, if it is available
	unsigned m_MinDist;
	///the number of layers for each polarizing transformation
	unsigned m_NumOfLayers;
	///size of the polarizing transformation, i.e. length of the unshortened and unpunctured code
	unsigned m_UnshortenedLength;
	///kernel for each layer
	const CBinaryKernel** m_ppKernels;
	///dynamic freezing constraints
	unsigned** m_ppFreezingConstraints;
	///if m_pDecisionPoints[i]!=~0, then u_i is (dynamic) frozen, and m_pDecisionPoints[i] gives the corresponding freezing constraint
	unsigned* m_pDecisionPoints;
	///input data  with frozen symbols	
	tBit* m_pEncodedInformation;
	///status of each codeword symbol
	SymbolType* m_pSymbolType;

	///temporary array of size m_UnshortenedLength
	tBit* m_pEncodingBuffer;
	///temporary array of size m_UnshortenedLength
	tBit* m_pEncodingBuffer2;
protected:
	///unpuncture and unshorten the code
	virtual void LoadLLRs(MType* pDest ///the LLR vector as required by the decoder
		, const MType* pSrc
	);
public:
	CMixedKernelEncoder(std::istream& SpecFile);
	virtual ~CMixedKernelEncoder();
	virtual unsigned GetMinimumDistance()const
	{
		return m_MinDist;
	}
	///encode a block of data into a codeword
	virtual void Encode(const tBit* pSrc ///the information symbols to be encoded
						,tBit* pEncoded ///the codeword
						);
	///shorten and puncture the original polar codeword
	void Shorten(const tBit* pPolarCW ///the original polar codeword
		, tBit* pShortened ///shortened and punctured codeword
	);
	///extract information bits from a codeword
	virtual void ExtractInformationBits(const tBit* pCodeword ///the codeword
										, tBit* pInfBits ///information symbols
	) ;
	///extract information bits from an unshortened and unpunctured codeword
	void ExtractInformationBitsUnshortened(const tBit* pCodeword ///the codeword
		, tBit* pInfBits ///information symbols
	);

	friend int main(int argc, char* argv[]);
};
