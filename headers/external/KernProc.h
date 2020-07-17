#pragma once
#include "Kernel.h"
#include "SoftProcessing.h"
///a generic LLR-domain processor for a binary kernel
class CKernProcLLR
{
protected:
	const CBinaryKernel* m_pKernel;
public:
	CKernProcLLR(const CBinaryKernel& K) :m_pKernel(&K)
	{

	}
	CKernProcLLR() {};
	virtual ~CKernProcLLR()
	{

	};
	virtual const char* GetName()const = 0;
	virtual unsigned Size()const
	{
		return m_pKernel->Size();
	}
	///allocate state variables if necessary
	virtual void* GetState(unsigned Stride  ///enable processing of data blocks of size Stride
						  )const
	{
		return 0;
	};
	///allocate temporary variables if necessary
	virtual void* GetTemp(unsigned Stride ///enable processing of data blocks of size Stride
	)const
	{
		return 0;
	}
	///deallocate state variables
	virtual void FreeState(void * State)const
	{

	};
	///deallocate temporary variables
	virtual void FreeTemp(void * Temp)const
	{

	};
	///copy the state variable
	virtual void CopyState(unsigned phase ///the number of known kernel input symbols
		, unsigned Stride ///the state variable is assumed to be configured for processing of Stride LLRs
		, void* pDest ///destination state variable 
		, const void* pSrc ///source state variable
	)const
	{

	};

	///compute LLRs for a given kernel phase
	virtual void GetLLRs(unsigned Stride ///number of LLRs to be computed 
		, unsigned phase ///local kernel phase
		, const tBit* pKnownInputSymbols ///known input symbols of the kernel. This is a vector of size Stride*m_pKernel->size
		, const MType* pChannelLLRs ///the LLRs for kernel output symbols. This is a vector of size Stride*m_pKernel->size
		, MType* pLLRs ///output: the LLRs for kernel input symbol phase. This is a vector of size Stride
		, void* pState = 0 ///state variable
		, void* pTemp = 0 ///temporary variable
	)const =0;
	virtual void GetKernelOutput(unsigned Stride  ///number of symbol blocks to be processed
		, const tBit* pKnownInputSymbols ///known kernel input symbol values
		, tBit* pOutput ///kernel output symbol values
	)const
	{
		for (unsigned i = 0; i < Stride*m_pKernel->Size(); i++)
			_ASSERT((pKnownInputSymbols[i] & ~BIT_1) == 0);
		m_pKernel->Multiply(Stride, pKnownInputSymbols, pOutput);
	}

};


class CArikanKernProc_S :public CKernProcLLR
{
public:
	CArikanKernProc_S():CKernProcLLR(ArikanKernel)
	{

	};
	virtual ~CArikanKernProc_S()
	{

	};
	void GetLLRs(unsigned Stride ///number of LLRs to be computed 
		, unsigned phase ///local kernel phase
		, const tBit* pKnownInputSymbols ///known input symbols of the kernel. This is a vector of size Stride*m_pKernel->size
		, const MType* pChannelLLRs ///the LLRs for kernel output symbols. This is a vector of size Stride*m_pKernel->size
		, MType* pLLRs ///output: the LLRs for kernel input symbol phase. This is a vector of size Stride
		, void* pState = 0 ///state variable
		, void* pTemp = 0 ///temporary variable
	)const
	{
		if (!phase)
			SoftXOR(pLLRs, pChannelLLRs, Stride);
		else
			SoftCombine(pLLRs, pChannelLLRs, pKnownInputSymbols, Stride);
	};
	virtual const char* GetName()const
	{
		return "Arikan";
	}
};