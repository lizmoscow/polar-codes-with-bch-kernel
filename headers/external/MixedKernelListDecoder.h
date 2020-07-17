#pragma once
#include "SeqConfig.h"
#include "KernelListEngine.h"

KType RandomizedSelect(tPairOfMTypeUnsigned* A, int p, int r, int i);
KType RandomizedPartition(tPairOfMTypeUnsigned* A, int, int);
KType Partition(tPairOfMTypeUnsigned* A, int p, int r);

#define USE_ORDERED_STATISTICS
class CMixedKernelListDecoder :public CListKernelEngine, public CBinarySoftDecoder
{
protected:
	tPairOfMTypeUnsigned* m_pSortingBuffer;
	tPairOfMTypeUnsigned* m_pSortingBuffer2;
	bool* m_pActivePaths;
	void ContinuePathsFrozen(unsigned phi);
	void ContinuePathsUnfrozen(unsigned phi);
	/*
#ifdef ZERO_SKIP
	///the first element is a phase, where zero advancement is possible. The second element is the number of layers one should advance
	tPairOfUnsigned* m_pZeroAdvancePhases;
	///process a block of static frozen symbols of a given depth. Return the size of the skipped block
	unsigned SkipZeroBlock(unsigned phi, unsigned Depth);
#endif
*/
public:
	CMixedKernelListDecoder(std::istream& Spec, unsigned ListSize);
	~CMixedKernelListDecoder();
	virtual unsigned GetMaxListSize()const
	{
		return m_MaxNumOfPaths;
	};

	//soft-input decoder. Produces  a list of information vectors
	//returns the number of vectors obtained, or -1 in case of error
	virtual int Decode(const MType* pLLR,//log(P{c_i=1|y_i}/P{c_i=0|y_i})
		tBit* pInfVectorList, //output array. Must be sufficiently big to accomodate GetMaxListSize() vectors of length GetDimension()
		tBit* pCodewordList = 0 //output array. If non-zero, must be sufficiently big to accomodate GetMaxListSize() vectors of length GetLength()
	);


};

