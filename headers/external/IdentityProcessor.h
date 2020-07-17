#pragma once

#include "KernProc.h"

///a processor for a few Arikan layers
class CIdentityProcessor :public CKernProcLLR
{
	unsigned m_Size;
public:
	CIdentityProcessor(unsigned Size) :CKernProcLLR(ArikanKernel), m_Size(Size)
	{
		m_pKernel = 0;///this is a dummy kernel
	};
	virtual ~CIdentityProcessor()
	{

	};
	virtual unsigned Size()const
	{
		return m_Size;
	}
	void GetLLRs(unsigned Stride ///number of LLRs to be computed 
		, unsigned phase ///local kernel phase
		, const tBit* pKnownInputSymbols ///known input symbols of the kernel. This is a vector of size Stride*m_pKernel->size
		, const MType* pChannelLLRs ///the LLRs for kernel output symbols. This is a vector of size Stride*m_pKernel->size
		, MType* pLLRs ///output: the LLRs for kernel input symbol phase. This is a vector of size Stride
		, void* pState = 0 ///state variable
		, void* pTemp = 0 ///temporary variable
	)const
	{
		memcpy(pLLRs, pChannelLLRs + phase*Stride, sizeof(MType)*Stride);
	}
	virtual void GetKernelOutput(unsigned Stride  ///number of symbol blocks to be processed
		, const tBit* pKnownInputSymbols ///known kernel input symbol values
		, tBit* pOutput ///kernel output symbol values
	)const
	{
		for (unsigned i = 0; i < Stride*m_Size; i++)
		{
			_ASSERT((pKnownInputSymbols[i] & ~BIT_1) == 0);
		};
		memcpy(pOutput, pKnownInputSymbols, sizeof(tBit)*Stride*m_Size);
	}
	const char* GetName()const
	{
		return "Identity";
	}

};

