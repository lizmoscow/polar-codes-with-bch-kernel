#ifndef MODEM_H
#define MODEM_H

#include <complex>
#include "misc.h"
#include "SeqConfig.h"

typedef std::complex<double> Complex;


//generic modulator and demodulator
class DLL_API CModem
{
protected:
	//number of bits per symbol
	unsigned m_BitsPerSymbol;
	//noise variance
	double m_NoiseVariance;
    // discretization coefficient
    double m_DiscretizationCoeff;

	CModem(unsigned BitsPerSymbol, bool NoLLRs = false) :
        m_BitsPerSymbol(BitsPerSymbol), 
        m_DiscretizationCoeff(1.0)
	{
	};
public:
	virtual ~CModem()
	{
	}
	unsigned GetBitsPerSymbol()const
	{
		return m_BitsPerSymbol;
	}
	void SetNoiseVariance(double Sigma2)
	{
		m_NoiseVariance = Sigma2;
	}
    void SetDiscretizationCoefficient(double DC)
    {
        m_DiscretizationCoeff = DC;
    }
	///modulate a symbol
	virtual Complex Modulate( const tBit* pSrc ///source bits
		) = 0;
	///demodulate a number of symbols, and compute the log-likelihood ratios 
	///the log-likelihood ratios are defined as log P(1|y)/P(0|y)
	virtual void Demodulate(const Complex& Noisy ///noisy symbols
		, const Complex& ChannelData ///channel fading factor
		, MType * pLLRs ///LLRs for each bit
		) = 0;
};

class DLL_API CBPSKModem :public CModem
{
public:
	CBPSKModem() :CModem(1)
	{

	};
	Complex Modulate(const tBit* pSrc ///source bits
		)
	{
		return 2.0*pSrc[0] - 1;
	};
	void Demodulate(const Complex& Noisy ///noisy symbols. Only the real component will be taken into account
		, const Complex& ChannelData ///channel fading factor. Only the real component will be taken into account
		, MType * pLLRs ///LLRs for each bit
		)
	{
#if defined(NORMALIZE_MODEM)
#if defined(USE_INT8)
        pLLRs[0] = static_cast<MType>(- ChannelData.real()*Noisy.real()*m_DiscretizationCoeff);
#else
        pLLRs[0] = static_cast<MType>(- ChannelData.real()*Noisy.real());
#endif
#else
		pLLRs[0] =(MType) (-2 * ChannelData.real()*Noisy.real() / m_NoiseVariance);
#endif
	};
};

class DLL_API CPAMModem :public CModem
{
	double* m_pConstellation;
	double* m_pDist;
public:
	CPAMModem(unsigned BitPerSymbol);

	virtual ~CPAMModem();
	Complex Modulate(const tBit* pSrc ///source bits
		)
	{
		unsigned Bits = 0;
		for (unsigned i = 0; i < m_BitsPerSymbol; i++)
			Bits += pSrc[i] << i;
		return m_pConstellation[Bits];
	};

	void Demodulate(const Complex& Noisy ///noisy symbols. Only the real component will be taken into account
		, const Complex& ChannelData ///channel fading factor. Only the real component will be taken into account
		, MType * pLLRs ///LLRs for each bit
		);

};
#endif