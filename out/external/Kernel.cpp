#include <fstream>
#include "LinAlg.h"
#include "Kernel.h"
#include "TrellisKernelProcessor.h"
#include "TruncatedArikanKernel.h"
#include "IdentityProcessor.h"
using namespace std;
const tBit pArikan[4] =
{
	1,0,
	1,1
};

const tBit* CArikanKernel::GetKernel()const
{
	return pArikan;
}
CArikanKernel ArikanKernel;


///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l, 0\leq i<d
///where l is kernel dimension (m_Size), and d is stride
void CArikanKernel::Multiply(unsigned Stride ///d
	, const tBit* pSrc ///input data (x)
	, tBit* pDest
)const
{
	memcpy(pDest + Stride, pSrc + Stride, sizeof(tBit)*Stride);
	XOR(pDest, pSrc, pSrc + Stride, Stride);
}

///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l^T, 0\leq i<d
///where l is kernel dimension (m_Size), and d is stride
void CArikanKernel::MultiplyTransposed(unsigned Stride ///d
	, const tBit* pSrc ///input data (x)
	, tBit* pDest
)const
{
	memcpy(pDest, pSrc, sizeof(tBit)*Stride);
	XOR(pDest + Stride, pSrc, pSrc + Stride, Stride);
}

///get kernel name
const char*CArikanKernel::Name()const
{
	return "A";
}

double FrankDelta(double C)
{
	return -fabs(C - 0.5) / 32 + 1.0 / 64;
};

/**
See R1-1706674 – “FRANK polar construction: nested extension design of Polar codes based on mutual information ”
*/
double CArikanKernel::GetSubchannelCapacity(eChannel  Channel  ///channel model to be used
	, unsigned SubchannelID ///ID of the subchannel to be analyzed
	, double C ///capacity of the underlying channel
)const
{
	switch (Channel)
	{
	case ec_AWGN:
		if (SubchannelID == 0)
			return C*C + FrankDelta(C);
		else
			return 2 * C - C*C - FrankDelta(C);
		break;
	case ec_BEC:
		if (SubchannelID == 0)
			return C*C ;
		else
			return 2 * C - C*C;
		break;
	default:
		throw Exception("Cannot evaluate capacity for channel type %d", (unsigned)Channel);
	}
}


CKernProcLLR* CArikanKernel::GetProcessor(unsigned ProcessorID)const
{
	switch (ProcessorID)
	{
	case 0: return new CArikanKernProc_S;
	case 1:  return new CIdentityProcessor(2);
	default: throw Exception("Invalid processor ID %d", ProcessorID);
	}
}


CMatrixBinaryKernel::CMatrixBinaryKernel(const char* pFileName):m_DoNotDelete(false)
{
	ifstream KernelFile(pFileName);
	KernelFile >> m_Size;
	if (!KernelFile)
		throw Exception("Error reading kernel file %s", pFileName);
	m_pKernel = AlignedAlloc(2 * m_Size*m_Size);
	for (unsigned i = 0; i < m_Size*m_Size; i++)
	{
		unsigned c;
		KernelFile >> c;
		m_pKernel[i] = c;
	}
	if (!KernelFile)
		throw Exception("Error parsing kernel file %s", pFileName);
	size_t L = strlen(pFileName);
	m_pName = new char[L + 2];
	strcpy_s(m_pName+1, L+1, pFileName);
	m_pName[0] = '-';
	Init();
	memset(m_ppCapacityEvaluators, 0, sizeof(m_ppCapacityEvaluators));
	//try to load capacity data for various channels
	bool Failed=false;
	do
	{
		unsigned ChannelType;
		KernelFile >> ChannelType;
		if (!KernelFile||ChannelType>=ec_NUMOFCHANNELS)
			Failed = true;
		else
		{
			try
			{
				m_ppCapacityEvaluators[ChannelType] = LoadCapacityEvaluator((eChannel)ChannelType,m_Size, KernelFile);
			}
			catch (Exception &ex)
			{
				cerr << "Warning, loading capacity data for channel " << ChannelType << " from file " << pFileName << " failed:" << ex.what() << endl;;
				Failed = true;
			}
		}
	} while (!Failed);
}
///initialize the kernel to a pre-defined matrix. The kernel does not become owner of it
CMatrixBinaryKernel::CMatrixBinaryKernel(unsigned Size ///size of the kernel
	, const tBit* pKernel ///kernel matrix
	, const char* pName ///name of the kernel
	, const CKernelCapacityEvaluator**ppCapEval  ///capacity evaluators
):CBinaryKernel(Size),m_DoNotDelete(true)
{
	m_pKernel = AlignedAlloc(2 * Size*Size);
	memcpy(m_pKernel, pKernel, sizeof(tBit) * Size*Size);
	size_t L = strlen(pName);
	m_pName = new char[L + 1];
	strcpy_s(m_pName,L+1, pName);
	if (!ppCapEval)
		memset(m_ppCapacityEvaluators, 0, sizeof(m_ppCapacityEvaluators));
	else
		memcpy(m_ppCapacityEvaluators, ppCapEval, sizeof(CKernelCapacityEvaluator*)*ec_NUMOFCHANNELS);
	Init();
}

void CMatrixBinaryKernel::Init()
{
	//compute inverse of the kernel matrix
	tBit* pTemp = AlignedAlloc(2 * m_Size*m_Size);
	for (unsigned i = 0; i < m_Size; i++)
	{
		memcpy(pTemp + 2 * i*m_Size, m_pKernel + i*m_Size, m_Size * sizeof(tBit));
		memset(pTemp + (2 * i + 1)*m_Size, 0, sizeof(tBit)*m_Size);
		pTemp[(2 * i + 1)*m_Size + i] = 1;
	}
	unsigned R = m_Size;
	Gauss(pTemp, R, 2 * m_Size, 2 * m_Size, true, NULL);
	if (!pTemp[(m_Size - 1)*(2 * m_Size + 1)])
		throw Exception("Kernel is singular");
	m_pInverseKernel = AlignedAlloc(m_Size*m_Size);
	for (unsigned i = 0; i < m_Size; i++)
	{
		memcpy(m_pInverseKernel + i*m_Size, pTemp + (2 * i + 1)*m_Size, sizeof(tBit)*m_Size);
	};
	AlignedFree(pTemp);

};
CMatrixBinaryKernel::~CMatrixBinaryKernel()
{
	AlignedFree(m_pInverseKernel);
	delete[]m_pName;
	if (!m_DoNotDelete)
	{
		AlignedFree(m_pKernel);
		for (unsigned i = 0; i < ec_NUMOFCHANNELS; i++)
			delete m_ppCapacityEvaluators[i];
	};
};

double CMatrixBinaryKernel::GetSubchannelCapacity(
	eChannel  Channel  ///channel model to be used
	, unsigned SubchannelID ///ID of the subchannel to be analyzed
	, double PhysicalCapacity ///capacity of the underlying channel
)const
{
	if (m_ppCapacityEvaluators[Channel])
		return m_ppCapacityEvaluators[Channel]->GetSubchannelCapacity(PhysicalCapacity, SubchannelID);
	else throw Exception("Kernel %s does not support capacity evaluation for channel model %d", m_pName, (unsigned)Channel);
}

///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l, 0\leq i<d
///where l is kernel dimension (m_Size), and d is stride
void CMatrixBinaryKernel::Multiply(unsigned Stride ///d
	, const tBit* pSrc ///input data (x)
	, tBit* pDest
)const
{
	MatrixMultiply(m_Size, m_Size, Stride, pSrc, m_pKernel, pDest);
}

///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l^T, 0\leq i<d
///where l is kernel dimension (m_Size), and d is stride
void CMatrixBinaryKernel::MultiplyTransposed(unsigned Stride ///d
	, const tBit* pSrc ///input data (x)
	, tBit* pDest
)const
{
	MatrixMultiplyTransposed(m_Size, m_Size, Stride, pSrc, m_pKernel, pDest);
};
///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l^{-1}, 0\leq i<d
///where l is kernel dimension (m_Size), and d is stride
void CMatrixBinaryKernel::InverseMultiply(unsigned Stride ///d
	, const tBit* pSrc ///input data (x)
	, tBit* pDest
)const
{
	MatrixMultiply(m_Size, m_Size, Stride, pSrc, m_pInverseKernel, pDest);
}

CBinaryKernel* NewArikan()
{
	return new CArikanKernel;
}

const unsigned NUMOFKERNELS = 3;
const char* pKernels[NUMOFKERNELS] = { "A", "G","A5" };



CBinaryKernel* NewG()
{
	return new CGKernel();
}

CBinaryKernel* NewA5()
{
	return new CA5Kernel();
}

typedef CBinaryKernel* (*tNewKernel)();


tNewKernel KernelAllocators[NUMOFKERNELS] = { NewArikan,  NewG,NewA5};

///lookup the kernel by its name. If pName is prefixed by <, the kernel is loaded from the corresponding file
CBinaryKernel* GetKernelByName(const char* pName)
{
	if ((pName[0] == '<')|| (pName[0] == '-'))
	{
		pName++;
		CBinaryKernel* K = new CMatrixBinaryKernel(pName);
		return K;
	}
	else  
	{
		for (unsigned i = 0; i < NUMOFKERNELS; i++)
		{
			if (_stricmp(pName, pKernels[i]) == 0)
			{
				return KernelAllocators[i]();
			}
		};
		throw Exception("Unknown kernel %s",pName);
	}
}

CKernProcLLR* CMatrixBinaryKernel::GetProcessor(unsigned ProcessorID)const
{
	if (ProcessorID)
		throw Exception("Invalid processor ID %d",ProcessorID);
	return new CTrellisKernelProcessor(*this);
}

CKernelFactory::~CKernelFactory()
{
	tProcMap::iterator it2 = m_KerProcMap.begin();
	while (it2 != m_KerProcMap.end())
	{
		delete it2->second;
		it2++;
	}
	tKerMap::iterator it = m_KerMap.begin();
	while (it != m_KerMap.end())
	{
		delete it->second;
		it++;
	};
}

///get pointer to a kernel. The pointer is owned by the factory
const CBinaryKernel* CKernelFactory::GetKernelByName(const std::string& Name ///name of the kernel. If prefixed by <, the kernel is loaded from a file
)
{
	tKerMap::iterator it = m_KerMap.find(Name);
	if (it == m_KerMap.end())
	{
		CBinaryKernel* pK=::GetKernelByName(Name.c_str());
		m_KerMap[Name] = pK;
		return pK;
	}
	else
	{
		return  it->second;
	}
}
///get pointer to a kernel processor
const CKernProcLLR* CKernelFactory::GetKernelProcessor(const CBinaryKernel* pKernel ///the kernel obtained from GetKernelByName
													   ,unsigned ProcessorID ///ID of the processor
)
{
	tProcPair P(pKernel, ProcessorID);
	tProcMap::iterator it = m_KerProcMap.find(P);
	if (it == m_KerProcMap.end())
	{
		CKernProcLLR* pP = pKernel->GetProcessor(ProcessorID);
		m_KerProcMap[P] = pP;
		return pP;
	}
	else
		return it->second;
}



CKernelFactory KernelFactory;