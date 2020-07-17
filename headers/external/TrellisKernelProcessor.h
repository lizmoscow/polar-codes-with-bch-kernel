#pragma once
#include "KernProc.h"


class CTrellisKernelProcessor :public CKernProcLLR
{
protected:
	struct TrellisNode
	{
		//the index of the next state. If NextState[i]=~0, then the corresponding transition is invalid
		unsigned long long NextState[2];
	};
	///for each kernel phase we need to construct an extended trellis
	TrellisNode*** m_pppTrellis;
	///number of active bits at each layer of the trellis
	unsigned** m_ppNumOfActiveBits;
	///maximal number of active bits
	unsigned m_MaxNumOfAB;
public:
	CTrellisKernelProcessor(const CMatrixBinaryKernel& K);
	virtual ~CTrellisKernelProcessor();

	///allocate state variables if necessary
	virtual void* GetState(unsigned Stride  ///enable processing of data blocks of size Stride
	)const;
	///allocate temporary variables if necessary
	virtual void* GetTemp(unsigned Stride ///enable processing of data blocks of size Stride
	)const;
	///deallocate state variables
	virtual void FreeState(void * State)const;
	///deallocate temporary variables
	virtual void FreeTemp(void * Temp)const;
	///copy the state variable
	virtual void CopyState(unsigned phase ///the number of known kernel input symbols
		, unsigned Stride ///the state variable is assumed to be configured for processing of Stride LLRs
		, void* pDest ///destination state variable 
		, const void* pSrc ///source state variable
	)const;


	///compute LLRs for a given kernel phase
	virtual void GetLLRs(unsigned Stride ///number of LLRs to be computed 
		, unsigned phase ///local kernel phase
		, const tBit* pKnownInputSymbols ///known input symbols of the kernel. This is a vector of size Stride*m_pKernel->size
		, const MType* pChannelLLRs ///the LLRs for kernel output symbols. This is a vector of size Stride*m_pKernel->size
		, MType* pLLRs ///output: the LLRs for kernel input symbol phase. This is a vector of size Stride
		, void* pState  ///state variable
		, void* pTemp ///temporary variable
	)const;
	virtual const char* GetName()const
	{
		return "TrellisProcessor";
	}


};