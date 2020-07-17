#include <gsl/gsl_math.h>
#include "misc.h"


 void *operator new(size_t Size)
{
	void *Ptr = _aligned_malloc(Size, ALIGN);
	if (!Ptr) throw std::bad_alloc();
#ifdef MEMORY_COUNTING
	CNT(CurrentCounter, Size);
#endif
	return Ptr;
}
 void *operator new[](size_t Size)
{
	void *Ptr = _aligned_malloc(Size, ALIGN);
	if (!Ptr) throw std::bad_alloc();
#ifdef MEMORY_COUNTING
	CNT(CurrentCounter, Size);
#endif
	return Ptr;
}

 void operator delete(void* ptr)
{
	_aligned_free(ptr);
};
 void operator delete[](void* ptr)
{
	_aligned_free(ptr);
};

std::jmp_buf JUMP_BUFFER;
char ErrorMessageBuffer[ERROR_MESSAGE_BUFFER_SIZE];

//operation counters
unsigned long long SumCount = 0;
unsigned long long CmpCount = 0;
unsigned long long MulCount = 0;
unsigned long long ExpCount = 0;
unsigned long long CpyCount = 0;
unsigned long long pCounters[CNT_END] = {0};
Counters CurrentCounter = CNT_UNKNOWN;

#define GENPOLY 0x18005 
#define CRCSIZE 16

unsigned CRC16(const tBit* pSrc, unsigned K)
{
	//calculate CRC-16
	unsigned CRCBuffer = 0;
	for (unsigned i = 0; i<K; i++)
	{
		CRCBuffer <<= 1;
		CRCBuffer |= (pSrc[i]>0);
		if (CRCBuffer >> CRCSIZE)
			CRCBuffer ^= GENPOLY;
	};
	for (unsigned i = 0; i<CRCSIZE; i++)
	{
		CRCBuffer <<= 1;
		if (CRCBuffer >> CRCSIZE)
			CRCBuffer ^= GENPOLY;
	};

	return CRCBuffer;
};


/*Вычисление биномиального коэффициента C_n^k*/
unsigned long long Binomial(unsigned n, unsigned k)
{
	unsigned __int64 B = 1;
	unsigned i;
	if (k>n)
		return 0;
	for (i = 1; i <= n - k; i++)
	{
		B *= k + i;
		B /= i;
	};
	return B;
}

void neg(MType * value) {
#if defined(USE_INT8)
	*value = (*value == '\x80') ? '\x7F' : -(*value);
#else
	*value = -(*value);
#endif
}

MType umin(const MType & value) {
#if defined(USE_INT8)
	return (value == '\x80') ? '\x7F' : -value;
#else
	return -value;
#endif
}

int ceilLog2(unsigned a)
{
	if (!a) return -1;
	unsigned m = (a & (a - 1)) ? 1 : 0;
	while (a >>= 1)
	{
		m++;
	}
	return m;
};

unsigned binomialCoeff(unsigned n, unsigned k){
	unsigned res = 1;
	if (k > n - k) k = n - k;

	for (unsigned i = 0; i < k; ++i){
		res = res * (n - i);
		res = res / (i + 1);
	}

	return res;
}




// Coefficients for Hermite integration
double Hermite_wp[NUMOFHERMITEINTEGRATIONPOINTS] = {
	7.64043285560789861e-06,
	1.34364574684723431e-03,
	3.38743944571450600e-02,
	2.40138611094110666e-01,
	6.10862633765332230e-01,
	6.10862633765332230e-01,
	2.40138611094110666e-01,
	3.38743944571450600e-02,
	1.34364574684723431e-03,
	7.64043285560789861e-06
};

//evaluation points for Hermite integration
double Hermite_xp[NUMOFHERMITEINTEGRATIONPOINTS] = {
	3.43615911883773739e+00,
	2.53273167423278966e+00,
	1.75668364929988163e+00,
	1.03661082978951358e+00,
	3.42901327223704588e-01,
	-3.42901327223704588e-01,
	-1.03661082978951358e+00,
	-1.75668364929988163e+00,
	-2.53273167423278966e+00,
	-3.43615911883773739e+00
};

double GetSubsetCapacity(double Sigma,
	unsigned SubconstellationSize, const double* pSubConstellation)
{
	//Use Hermite integration
	double Sum = 0;
	for (unsigned i = 0; i<SubconstellationSize; i++)
	{
		//evaluate the integral using substitution z=(y-a_i)/sqrt(2*sigma^2) 
		double CurIntegral = 0;
		double CurSubConstPoint = pSubConstellation[i];
		for (unsigned k = 0; k<NUMOFHERMITEINTEGRATIONPOINTS; k++)
		{
			//evaluate PDF of the signal 
			double CurEvaluationPoint = Hermite_xp[k];
			double CurWeight = Hermite_wp[k];

			double Numenator = exp(-CurEvaluationPoint*CurEvaluationPoint);
			double Denumenator = 0;
			for (unsigned j = 0; j<SubconstellationSize; j++)
			{
				double Temp = CurEvaluationPoint + (pSubConstellation[i] - pSubConstellation[j]) / (M_SQRT2*Sigma);
				Denumenator += exp(-Temp*Temp);
			};
			CurIntegral += CurWeight*log(Numenator*SubconstellationSize / Denumenator) / M_LN2;
		};
		Sum += CurIntegral;
	};
	return Sum / (sqrt(M_PI)*SubconstellationSize);
};


double GetBPSKAWGNChannelEntropy(double Sigma)
{
	double BPSK[2] = { 1,-1 };
	//Use Hermite integration	
	double Sum = 0;
	for (unsigned i = 0; i<2; i++)
	{
		//evaluate the integral using substitution z=(y-a_i)/sqrt(2*sigma^2) 
		double CurIntegral = 0;
		double CurSubConstPoint = BPSK[i];
		for (unsigned k = 0; k<NUMOFHERMITEINTEGRATIONPOINTS; k++)
		{
			//evaluate PDF of the signal 
			double CurEvaluationPoint = Hermite_xp[k];
			double CurWeight = Hermite_wp[k];

			double Numenator = exp(-CurEvaluationPoint*CurEvaluationPoint);
			double Denumenator = 0;
			for (unsigned j = 0; j<2; j++)
			{
				double Temp = CurEvaluationPoint + (BPSK[i] - BPSK[j]) / (M_SQRT2*Sigma);
				Denumenator += exp(-Temp*Temp);
			};
			CurIntegral += CurWeight*log(Numenator/ Denumenator) / M_LN2;
		};
		Sum -= CurIntegral;
	};
	return Sum / (sqrt(M_PI)*2);
};
