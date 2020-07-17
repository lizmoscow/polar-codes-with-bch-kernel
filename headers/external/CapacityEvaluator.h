#pragma once
#include <iostream>
#include <gsl/gsl_interp.h>
#include "Simulator.h"

class CKernelCapacityEvaluator
{
protected:
	unsigned m_KernelSize;
public:
	CKernelCapacityEvaluator(unsigned Size) :m_KernelSize(Size)
	{

	};
	virtual ~CKernelCapacityEvaluator()
	{

	};
	///compute capacity (or some its surrogate) of a subchannel induced by the kernel, given the capacity of the underlying physical channel
	virtual double GetSubchannelCapacity(double PhysCap ///capacity of the underlying channel
		, unsigned SubchannelID ///ID of the subchannel
	) = 0;
};

///compute subchannel capacities using Gaussian approximation based on interpolation tables
class CAWGNCapacityInterpolator :public CKernelCapacityEvaluator
{
	///number of entries in the capacity table
	unsigned m_NumOfCapPoints;
	double* m_pCapacities;
	double* m_pPhysCapacities;
	bool m_DoNotDelete;
	gsl_interp** m_ppInterpolators;
	gsl_interp_accel** m_ppAcceletators;
	void Init(unsigned NumOfValidPoints=0);
public:
	CAWGNCapacityInterpolator(unsigned KernelSize, std::istream& DataFile);
	CAWGNCapacityInterpolator(unsigned KernelSize, unsigned NumOfPoints, double* pCapacities, double* pPhysCapacities);
	~CAWGNCapacityInterpolator();
	///compute capacity of a subchannel induced by the kernel, given the capacity of the underlying physical channel
	virtual double GetSubchannelCapacity(double PhysCap ///capacity of the underlying channel
		, unsigned SubchannelID ///ID of the subchannel
	);
};

///compute BEC subchannel capacities
class CBECCapacityEvaluator :public CKernelCapacityEvaluator
{
	///number of erasure patterns killing each subchannel
	double* m_pNumOfUncorrectableErasurePatterns;
	bool m_DoNotDelete;
public:
	CBECCapacityEvaluator(unsigned KernelSize, std::istream& DataFile);
	CBECCapacityEvaluator(unsigned KernelSize, double* pNumOfUncorrectableEP) :CKernelCapacityEvaluator(KernelSize), 
				m_pNumOfUncorrectableErasurePatterns(pNumOfUncorrectableEP), 
		m_DoNotDelete(true)
	{

	};
	~CBECCapacityEvaluator();
	///compute capacity of a subchannel induced by the kernel, given the capacity of the underlying physical channel
	virtual double GetSubchannelCapacity(double PhysCap ///capacity of the underlying channel
		, unsigned SubchannelID ///ID of the subchannel
	);

};

CKernelCapacityEvaluator* LoadCapacityEvaluator(eChannel Type, unsigned KernelSize,std::istream& DataFile);
