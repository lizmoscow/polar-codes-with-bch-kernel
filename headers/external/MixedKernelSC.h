#pragma once

#include "Codec.h"
#include "TrellisKernelProcessor.h"
#include "MixedKernelEncoder.h"
#include "TVMemoryEngine.h"

class CMixedKernelSCCodec :public CMixedKernelEncoder, public CBinarySoftDecoder
{
	MType** m_ppS;
	tBit** m_ppC;
	void** m_ppState;
	void** m_ppTemp;
	const CKernProcLLR** m_ppLLREngines;
	
	void RecursivelyCalcS(unsigned phi, unsigned m,unsigned Stride=1);
	void RecursivelyUpdateC(unsigned phi, unsigned m, unsigned Stride = 1);
public:
	CMixedKernelSCCodec(std::istream& CodeSpecification);
	virtual ~CMixedKernelSCCodec();
	//soft-input decoder. Produces  a list of information vectors
	//returns the number of vectors obtained, or -1 in case of error
	virtual int Decode(const MType* pLLR,//log(P{c_i=1|y_i}/P{c_i=0|y_i})
		tBit* pInfVectorList, //output array. Must be sufficiently big to accomodate GetMaxListSize() vectors of length GetDimension()
		tBit* pCodewordList = 0 //output array. If non-zero, must be sufficiently big to accomodate GetMaxListSize() vectors of length GetLength()
	);
	unsigned* m_pNumOfSubchannelErrors;
};