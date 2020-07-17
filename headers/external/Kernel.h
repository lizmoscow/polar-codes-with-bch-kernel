#pragma once
#include <map>
#include "Codec.h"
#include "Simulator.h"
#include "CapacityEvaluator.h"

class CKernProcLLR;
enum eFreezingState {efsFullyFrozen  ///all descendants of a symbol are frozen
	,efsPartiallyFrozen ///some children are frozen,  some are not
	,efsFullyUnfrozen ///all descendants are unfrozen
	,efsConstrained ///there are some constraints on this symbol
};

class CBinaryKernel
{
protected:
	unsigned m_Size;
	CBinaryKernel() :m_Size(0)
	{

	}
public:
	CBinaryKernel(unsigned Size) :m_Size(Size)
	{

	}
	virtual ~CBinaryKernel()
	{

	};
	virtual const tBit* GetKernel()const = 0;
	///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l, 0\leq i<d
	///where l is kernel dimension (m_Size), and d is stride
	virtual void Multiply(unsigned Stride ///d
		, const tBit* pSrc ///input data (x)
		, tBit* pDest
	)const = 0;
	///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l^T, 0\leq i<d
	///where l is kernel dimension (m_Size), and d is stride
	virtual void MultiplyTransposed(unsigned Stride ///d
		, const tBit* pSrc ///input data (x)
		, tBit* pDest
	)const = 0;

	///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l^{-1}, 0\leq i<d
	///where l is kernel dimension (m_Size), and d is stride
	virtual void InverseMultiply(unsigned Stride ///d
		, const tBit* pSrc ///input data (x)
		, tBit* pDest
	)const = 0;

	///get the dimension of the kernel
	unsigned Size()const
	{
		return m_Size;
	}
	///allocate an LLR computation engine
	virtual CKernProcLLR* GetProcessor(unsigned ProcessorID)const = 0;
	///return the ID of the kernel processor to be used
	virtual unsigned GetProcessorID(const eFreezingState* pFrozenFlags)const
	{
		return 0;
	};
	///compute the capacity of a subchannel induced by the kernel for a given physical channel capacity
	virtual double GetSubchannelCapacity(eChannel  Channel  ///channel model to be used
										,unsigned SubchannelID ///ID of the subchannel to be analyzed
									    , double PhysicalCapacity ///capacity of the underlying channel
	)const
	{
		throw Exception("Capacity functions are not known for kernel %s",Name());
	}
	///get kernel name
	virtual const char*Name()const = 0;

};


class CArikanKernel :public CBinaryKernel
{
public:
	CArikanKernel() :CBinaryKernel(2)
	{

	};
	const tBit* GetKernel()const;
	///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+s},x_{i+2s},...,x_{i+(l-1)s})F_l, 0\leq i<d
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

	///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+s},x_{i+2s},...,x_{i+(l-1)s})F_l^{-1}, 0\leq i<s
	///where l is kernel dimension (m_Size), and d is stride
	virtual void InverseMultiply(unsigned Stride ///d
		, const tBit* pSrc ///input data (x)
		, tBit* pDest
	)const
	{
		return Multiply(Stride, pSrc, pDest);
	}
	///allocate an LLR computation engine
	virtual CKernProcLLR* GetProcessor(unsigned ProcessorID)const;
	///get kernel name
	const char*Name()const;
	///compute the capacity of a subchannel induced by the kernel for a given physical channel capacity
	virtual double GetSubchannelCapacity(eChannel  Channel  ///channel model to be used
		,  unsigned SubchannelID ///ID of the subchannel to be analyzed
		, double PhysicalCapacity ///capacity of the underlying channel
	)const;
	///return the ID of the kernel processor to be used
	virtual unsigned GetProcessorID(const eFreezingState* pFrozenFlags)const
	{
		return pFrozenFlags[0]==pFrozenFlags[1] &&((pFrozenFlags[0] == efsFullyUnfrozen)||(pFrozenFlags[0] == efsFullyFrozen));
	};


};

extern CArikanKernel ArikanKernel;


///an arbitrary binary kernel given by some matrix
class CMatrixBinaryKernel :public CBinaryKernel
{
	char* m_pName;
	bool m_DoNotDelete;
	CKernelCapacityEvaluator* m_ppCapacityEvaluators[ec_NUMOFCHANNELS];
protected:
	 tBit* m_pKernel;
	///inverse of the kernel matrix;
	tBit* m_pInverseKernel;
	//initialize the inverse kernel
	void Init();
public:
	///initialize the kernel to a pre-defined matrix. The kernel does not become owner of it
	CMatrixBinaryKernel(unsigned Size ///size of the kernel
					   ,const tBit* pKernel ///kernel matrix
					   ,const char* pName ///name of the kernel
					   ,const CKernelCapacityEvaluator**ppCapEval=0  ///capacity evaluators
					   );
	CMatrixBinaryKernel(const char* pFileName);
	virtual ~CMatrixBinaryKernel();

	const tBit* GetKernel()const
	{
		return m_pKernel;
	}
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

	///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+s},x_{i+2s},...,x_{i+(l-1)s})F_l^{-1}, 0\leq i<s
	///where l is kernel dimension (m_Size), and d is stride
	virtual void InverseMultiply(unsigned Stride ///d
		, const tBit* pSrc ///input data (x)
		, tBit* pDest
	)const;
	///allocate an LLR computation engine
	virtual CKernProcLLR* GetProcessor(unsigned ProcessorID)const;
	const char*Name()const
	{
		return m_pName;
	}
	///compute the capacity of a subchannel induced by the kernel for a given physical channel capacity
	virtual double GetSubchannelCapacity(
		eChannel  Channel  ///channel model to be used
		,unsigned SubchannelID ///ID of the subchannel to be analyzed
		, double PhysicalCapacity ///capacity of the underlying channel
	)const;

};


class CKernelFactory
{
	typedef std::map<std::string, CBinaryKernel*>  tKerMap;
	typedef std::pair<const CBinaryKernel*, unsigned> tProcPair;
	typedef std::map<tProcPair,CKernProcLLR*>  tProcMap;
	///the mapping from kernel name to kernel processors
	tKerMap m_KerMap;
	tProcMap m_KerProcMap;
public:
	CKernelFactory()
	{

	};
	~CKernelFactory();
	///get pointer to a kernel. The pointer is owned by the factory
	const CBinaryKernel* GetKernelByName(const std::string& Name ///name of the kernel. If prefixed by < or -, the kernel is loaded from a file
						);
	///get pointer to a kernel processor
	const CKernProcLLR* GetKernelProcessor(const CBinaryKernel* pKernel ///the kernel obtained from GetKernelByName
										   ,unsigned ProcessorID ///ID of the kernel processor
						);


};
extern CKernelFactory KernelFactory;
///lookup the kernel by its name. If pName is prefixed by <, the kernel is loaded from the corresponding file
CBinaryKernel* GetKernelByName(const char* pName); 


extern const unsigned NUMOFKERNELS;
//names of known kernels
extern const char* pKernels[];