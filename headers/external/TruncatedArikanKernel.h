#pragma once
#include "Kernel.h"
#include "KernProc.h"
/**
Modified Gabry-Bioglio-Land-Belfiore kernel

100
110
011

*/

class CGKernel :public CBinaryKernel
{
	CKernelCapacityEvaluator* m_ppCapacityEvaluators[ec_NUMOFCHANNELS];
public:
	CGKernel() ;
	///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l, 0\leq i<d
	///where l is kernel dimension (m_Size), and d is stride
	virtual void Multiply(unsigned Stride ///d
		, const tBit* pSrc ///input data (x)
		, tBit* pDest
	)const;
	///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l^T, 0\leq i<d
	///where l is kernel dimension (m_Size), and d is stride
	virtual void MultiplyTransposed(unsigned Stride ///d
		, const tBit* pSrc ///input data (x)
		, tBit* pDest
	)const;

	///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l^{-1}, 0\leq i<d
	///where l is kernel dimension (m_Size), and d is stride
	virtual void InverseMultiply(unsigned Stride ///d
		, const tBit* pSrc ///input data (x)
		, tBit* pDest
	)const;
	const tBit* GetKernel()const;
	///allocate an LLR computation engine
	virtual CKernProcLLR* GetProcessor(unsigned ProcessorID)const;
	///get kernel name
	const char*Name()const;
	///compute the capacity of a subchannel induced by the kernel for a given physical channel capacity
	virtual double GetSubchannelCapacity(eChannel  Channel  ///channel model to be used
		, unsigned SubchannelID ///ID of the subchannel to be analyzed
		, double PhysicalCapacity ///capacity of the underlying channel
	)const;
	///check if simplified LLR processing is possible
	virtual unsigned GetProcessorID(const eFreezingState* pFrozenFlags)const
	{
		bool AllTheSame = (pFrozenFlags[0] == pFrozenFlags[1]) && (pFrozenFlags[0] == pFrozenFlags[2]) ;
		return AllTheSame  && (pFrozenFlags[0] == efsFullyFrozen || pFrozenFlags[0] == efsFullyUnfrozen);
		//return 0;
	};


};

extern CGKernel GKernel;

class CGKernelProcessor :public CKernProcLLR
{
public:
	CGKernelProcessor() :CKernProcLLR(GKernel)
	{

	};
	virtual ~CGKernelProcessor()
	{

	};
	///allocate state variables if necessary
	virtual void* GetState(unsigned Stride  ///enable processing of data blocks of size Stride
	)const;
	///deallocate state variables
	virtual void FreeState(void * State)const;
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
		, void* pState = 0 ///state variable
		, void* pTemp = 0 ///temporary variable
	)const;
	const char* GetName()const
	{
		return "GProc";
	}
};

class CA5Kernel :public CBinaryKernel
{
	CKernelCapacityEvaluator* m_ppCapacityEvaluators[ec_NUMOFCHANNELS];
public:
	CA5Kernel();
	///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l, 0\leq i<d
	///where l is kernel dimension (m_Size), and d is stride
	virtual void Multiply(unsigned Stride ///d
		, const tBit* pSrc ///input data (x)
		, tBit* pDest
	)const;
	///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l^T, 0\leq i<d
	///where l is kernel dimension (m_Size), and d is stride
	virtual void MultiplyTransposed(unsigned Stride ///d
		, const tBit* pSrc ///input data (x)
		, tBit* pDest
	)const;

	///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l^{-1}, 0\leq i<d
	///where l is kernel dimension (m_Size), and d is stride
	virtual void InverseMultiply(unsigned Stride ///d
		, const tBit* pSrc ///input data (x)
		, tBit* pDest
	)const;
	const tBit* GetKernel()const;
	///allocate an LLR computation engine
	virtual CKernProcLLR* GetProcessor(unsigned ProcessorID)const;
	///get kernel name
	const char*Name()const;
	///compute the capacity of a subchannel induced by the kernel for a given physical channel capacity
	virtual double GetSubchannelCapacity(eChannel  Channel  ///channel model to be used
		, unsigned SubchannelID ///ID of the subchannel to be analyzed
		, double PhysicalCapacity ///capacity of the underlying channel
	)const;
	///check if simplified LLR processing is possible
	virtual unsigned GetProcessorID(const eFreezingState* pFrozenFlags)const
	{
		bool AllTheSame = (pFrozenFlags[0] == pFrozenFlags[1]) && (pFrozenFlags[0] == pFrozenFlags[2]) && (pFrozenFlags[0] == pFrozenFlags[3]) && (pFrozenFlags[0] == pFrozenFlags[4]);
		return AllTheSame && (pFrozenFlags[0] == efsFullyFrozen || pFrozenFlags[0] == efsFullyUnfrozen);
	};

};
extern CA5Kernel A5Kernel;
class CA5KernelProcessor :public CKernProcLLR
{
public:
	CA5KernelProcessor() :CKernProcLLR(A5Kernel)
	{

	};
	virtual ~CA5KernelProcessor()
	{

	};
	///allocate state variables if necessary
	virtual void* GetState(unsigned Stride  ///enable processing of data blocks of size Stride
	)const;
	///deallocate state variables
	virtual void FreeState(void * State)const;
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
		, void* pState = 0 ///state variable
		, void* pTemp = 0 ///temporary variable
	)const;
	const char* GetName()const
	{
		return "A5Proc";
	}

};
