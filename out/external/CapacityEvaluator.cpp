#include "misc.h"
#include "CapacityEvaluator.h"
using namespace std;

CAWGNCapacityInterpolator::CAWGNCapacityInterpolator(unsigned KernelSize, std::istream& DataFile):CKernelCapacityEvaluator(KernelSize),m_DoNotDelete(false)
{
	//try to load the interpolation table
	DataFile >> m_NumOfCapPoints;
	if (!DataFile)
		throw Exception("Cannot load capacity data");
	m_NumOfCapPoints += 2;
	m_pPhysCapacities = new double[m_NumOfCapPoints];
	m_pCapacities = new double[m_NumOfCapPoints*m_KernelSize];
	//intialize boundaries
	m_pPhysCapacities[0] = 0;
	m_pPhysCapacities[m_NumOfCapPoints - 1] = 1;
	for (unsigned i = 0; i < m_KernelSize; i++)
	{
		m_pCapacities[i*m_NumOfCapPoints] = 0;
		m_pCapacities[(i + 1)*m_NumOfCapPoints - 1] = 1;
	};
	unsigned StopAt = m_NumOfCapPoints-1;
	for (unsigned i = 1; i < m_NumOfCapPoints - 1; i++)
	{
		DataFile >> m_pPhysCapacities[i];
		for (unsigned j = 0; j < m_KernelSize; j++)
		{
			DataFile >> m_pCapacities[j*m_NumOfCapPoints + i];
		}
		if (!DataFile)
			throw Exception("Error parsing capacity data");
		if (StopAt < m_NumOfCapPoints - 1) continue;
		if (m_pPhysCapacities[i] == m_pPhysCapacities[i - 1])
		{
			cerr << "Duplicate channel capacity entry " << m_pPhysCapacities[i] << " found, ignoring the remaining data in the file\n";
			StopAt = i;
		}
	};
	m_pPhysCapacities[StopAt- 1] = 1;
	for (unsigned i = 0; i < m_KernelSize; i++)
	{
		m_pCapacities[i*m_NumOfCapPoints +StopAt] = 1;
	};

	Init(StopAt);
}

void CAWGNCapacityInterpolator::Init(unsigned NumOfValidPoints)
{
	if (!NumOfValidPoints)
		NumOfValidPoints = m_NumOfCapPoints;
	m_ppInterpolators = new gsl_interp*[m_KernelSize];
	m_ppAcceletators = new gsl_interp_accel*[m_KernelSize];
	for (unsigned i = 0; i < m_KernelSize; i++)
	{
		m_ppInterpolators[i] = gsl_interp_alloc(gsl_interp_linear, NumOfValidPoints);

		gsl_interp_init(m_ppInterpolators[i], m_pPhysCapacities, m_pCapacities + i*m_NumOfCapPoints, NumOfValidPoints);
		m_ppAcceletators[i] = gsl_interp_accel_alloc();

	};

}

CAWGNCapacityInterpolator::CAWGNCapacityInterpolator(unsigned KernelSize, unsigned NumOfPoints, double* pCapacities, double* pPhysCapacities):CKernelCapacityEvaluator(KernelSize),
m_NumOfCapPoints(NumOfPoints),m_pCapacities(pCapacities),m_pPhysCapacities(pPhysCapacities),m_DoNotDelete(true)
{
	Init();
}

CAWGNCapacityInterpolator::~CAWGNCapacityInterpolator()
{
	if (!m_DoNotDelete)
	{
		delete[]m_pPhysCapacities;
		delete[]m_pCapacities;
	};
	for (unsigned i = 0; i < m_KernelSize; i++)
	{
		gsl_interp_free(m_ppInterpolators[i]);
		gsl_interp_accel_free(m_ppAcceletators[i]);
	};
	delete[]m_ppInterpolators;
	delete[]m_ppAcceletators;

}

double CAWGNCapacityInterpolator::GetSubchannelCapacity(double PhysicalCapacity ///capacity of the underlying channel
	, unsigned SubchannelID ///ID of the subchannel
)
{
	if (PhysicalCapacity > 1)
		PhysicalCapacity = 1;
	if (PhysicalCapacity <= 0)
		PhysicalCapacity = 0;
	return gsl_interp_eval(m_ppInterpolators[SubchannelID], m_pPhysCapacities, m_pCapacities + SubchannelID*m_NumOfCapPoints, PhysicalCapacity, m_ppAcceletators[SubchannelID]);

}

CBECCapacityEvaluator::CBECCapacityEvaluator(unsigned KernelSize, std::istream& DataFile):CKernelCapacityEvaluator(KernelSize)
{
	m_pNumOfUncorrectableErasurePatterns = new double[KernelSize*KernelSize];
	for (unsigned i = 0; i < KernelSize; i++)
	{
		double c;
		DataFile >> c;
		if (c != 0)
			throw Exception("Invalid polarization behaviour file");
		for (unsigned j = 0; j < KernelSize; j++)
		{
			DataFile >> m_pNumOfUncorrectableErasurePatterns[i*KernelSize + j];
		}
	};
	if (!DataFile)
		throw Exception("Error loading BEC capacity file");

}
CBECCapacityEvaluator::~CBECCapacityEvaluator()
{
	if (!m_DoNotDelete)
		delete[]m_pNumOfUncorrectableErasurePatterns;
}
///compute capacity of a subchannel induced by the kernel, given the capacity of the underlying physical channel
double CBECCapacityEvaluator::GetSubchannelCapacity(double PhysCap ///capacity of the underlying channel
	, unsigned SubchannelID ///ID of the subchannel
)
{
	double P = 1 - PhysCap;
	double S = P / PhysCap;
	//Horner rule for evaluation of the erasure probability
	double E = 0;
	double* pUEP = m_pNumOfUncorrectableErasurePatterns + SubchannelID*m_KernelSize;
	for (int i = m_KernelSize; i > 0; i--)
	{
		E *= S;
		E += pUEP[i - 1];
	};
	E *= S*pow(PhysCap,m_KernelSize);
	return 1-E;
}



CKernelCapacityEvaluator* LoadCapacityEvaluator(eChannel Type, unsigned KernelSize, std::istream& DataFile)
{
	switch (Type)
	{
	case ec_AWGN:return new CAWGNCapacityInterpolator(KernelSize, DataFile);
	case ec_BEC:return new CBECCapacityEvaluator(KernelSize, DataFile);
	default: throw Exception("Unsupported channel model %d", (unsigned)Type);
	}
}
