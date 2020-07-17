#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <chrono>
#include <iostream>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "Modem.h"
#include "Codec.h"

void ControlC_Handler(int dummy);

enum eChannel{ec_AWGN ///simulate AWGN channel
		     ,ec_RayleighReal ///simulate real-valued Rayleigh fading channel
			 ,ec_BSC /// simulate binary symmetric channel. NoiseStdDev stores error probability
	         ,ec_BEC /// simulate the binary erasure channel
			 ,ec_NUMOFCHANNELS  ///dummy value used to identify the number of known channel types
};

class CSimulator;

/// This function will be called after each simulation iteration. If it returns false, simulations will be terminated.
/// @param pEngine -- pointer to the simulator engine
/// @param DecodingFailed -- true if the last decoding attempt has failed
typedef  bool(*tSimulatorHook)(CSimulator* pEngine, bool DecodingFailed);

class CSimulator
{
public:
    static bool KeepRunning;
	gsl_rng* m_pRNG;
	CBinaryEncoder* m_pEncoder;
	CBinarySoftDecoder* m_pDecoder;
	CModem& m_Modem;
	//the data being transmitted
	tBit* m_pData;
	//the encoded data
	tBit* m_pEncoded;
	//LLRs at channel output
	MType * m_pLLRs;
	///decoded information symbols
	tBit* m_pDecodedData;
	///decoded codeword
	tBit* m_pDecodedCW;
	///channel output
	Complex* m_pNoisy;
	///channel input
	Complex* m_pModulated;
	// channel fading factors
	Complex *m_pChannel;

	//noise standard deviation
	double m_NoiseStdDev;
	double m_ErasureProbability;
	//the number of symbols to be transmitted
	unsigned m_NumOfSymbols;
	// decrement for the number of symbols being used for calculation of noise standard deviation
	int m_DecrementForNumOfSymbols;
	///the type of the channel to be used
	eChannel m_ChannelType;

	unsigned long long m_NumOfBitErrors;
	unsigned long long m_NumOfBlockErrors;
	unsigned long long m_NumOfMLErrors;
	unsigned long long m_NumOfRawBitErrors; // errors in channel
	unsigned long long m_NumOfCodeBitErrors; // codeword bit errors
	unsigned long long m_NumOfUndetectedError;
	unsigned long long m_NumOfListErrors;//number of events, when correct codeword is not available in the list
	unsigned long long m_IterationsDone;
	unsigned long long m_DecodingTime;

	unsigned long long m_MemoryUsage;
	unsigned long long m_MemoryUsageSuccess;
	unsigned long long m_MemoryUsageStatic;

	///run one simulation iteration. Return true if decoding is correct
	bool Iterate(unsigned long long It ///Iteration ID
				,bool Indeed ///true if the data should be indeed coded and decoded
				,const tBit* pData=0 ///do this for some specific data
			);
	CSimulator(CBinaryCodec& Codec ///the codec to be used
	          ,CModem& Modem ///the modem to be used
			  ,eChannel ChannelType
			  );
	virtual ~CSimulator();
	// compute decrement for the number of symbols being used for calculation of noise standard deviation 
	void LengthDecrement(unsigned LengthForEbN0);
	//set signal-to-noise ratio per bit
	void SetEbN0(double SNR);
	void SetEsN0(double SNR);
	void SetPrecomputedNoiseStdDev(double NoiseStdDev);
	void ResetCounters()
	{
		m_NumOfBitErrors = 0;
		m_NumOfBlockErrors = 0;
		m_NumOfMLErrors = 0;
		m_NumOfRawBitErrors = 0;
		m_NumOfCodeBitErrors = 0;
		m_NumOfUndetectedError = 0;
		m_IterationsDone = 0;
		m_DecodingTime = 0;
	}
	/// run simulations 
	/// @param iterations limit. stop after this number of iterations
	/// @param errors limit. stop after this number of errors
	/// @param hook. the function to be called after each iteration
	template<typename tHook = tSimulatorHook>
	void Simulate(unsigned long long MaxIterations, unsigned long long MaxErrors, tHook Hook);

	void RunOneIteration(unsigned IterationID ///do real encoding and decoding only for the specified iteration
		, const tBit* pData = 0 ///do this for some specific data pattern
		);

	virtual double BER()const
	{
		return double(m_NumOfBitErrors) / (m_IterationsDone*m_pEncoder->GetDimension());
	};
	virtual double FER()const
	{
		return double(m_NumOfBlockErrors) / m_IterationsDone;
	};
	virtual double MLFER()const
	{
		return double(m_NumOfMLErrors) / m_IterationsDone;
	};
	virtual double ListFER()const
	{
		return double(m_NumOfListErrors) / m_IterationsDone;
	}
	virtual double ChannelBER()const
	{
		return double(m_NumOfRawBitErrors) / (m_IterationsDone*m_pEncoder->GetLength());
	};
	virtual double CodeBER()const
	{
		return double(m_NumOfCodeBitErrors) / (m_IterationsDone*m_pEncoder->GetLength());
	};

	virtual double UndetectedFER()const
	{
		return double(m_NumOfUndetectedError) / (m_IterationsDone);
	}
	virtual double CodeFAR ()const
	{
		return double(m_NumOfUndetectedError) / (m_NumOfBlockErrors);
	}
	double GetAvgDecodingTime()const;
	void SetRngSeed(unsigned long Seed = time(NULL));

    void printData(const tBit* pData0, const tBit* pData1, unsigned N);
    void printDecodedData(const tBit* pDecodedData, unsigned K, const tBit* pDecodedCW, unsigned N);

	///set it to true if erroneous iteration numbers should be printed
	bool m_PrintErrorIterations;
	///set it to true if zero data should be used
	bool m_ZeroData;
    ///set it to true and data will contain no zeroes
    bool m_OneData;
	///set it to true to force assertion failure on non-ML decoding error
	bool m_AssertNonMLError;
	///set it to true if BER is not needed
	bool m_NoBER;

    bool m_PrintDecodedData;
    bool m_PrintCodeword;
};

template<typename tHook>
void CSimulator::Simulate(unsigned long long MaxIterations, unsigned long long MaxErrors, tHook Hook)
{
	for (m_IterationsDone; m_IterationsDone < MaxIterations;) {
		bool Correct = Iterate(m_IterationsDone, true);
		m_IterationsDone++;
		if (!Hook(this, Correct) || (m_NumOfBlockErrors >= MaxErrors) || !KeepRunning) {
			break;
		}
	}
}

class CTimeHook {
public:
	CTimeHook() : start(std::chrono::high_resolution_clock::now()), t0(start) {}

	bool operator()(CSimulator* pEngine, bool DecodingFailed)
	{
		using namespace std::chrono;
		t1 = high_resolution_clock::now();
		seconds elapsed = duration_cast<seconds>(t1 - t0);
		if (elapsed.count() >= 5) {
			high_resolution_clock::duration dur = t1 - start;
			hours h = duration_cast<hours>(dur);
			minutes m = duration_cast<minutes>(dur);
			seconds s = duration_cast<seconds>(dur);
			std::cout << pEngine->m_IterationsDone << ' ' 
				<< std::scientific << pEngine->FER() << " : "
				<< h.count() << "h " << m.count() % 60 << "m " << s.count() % 60 << "s"
				<< std::endl;
			t0 = t1;
		}
		return true; // do nothing, continue simulations
	}

private:
	std::chrono::high_resolution_clock::time_point start;
	std::chrono::high_resolution_clock::time_point t0;
	std::chrono::high_resolution_clock::time_point t1;
};

class CIterationHook {
public:
	bool operator()(CSimulator* pEngine, bool DecodingFailed)
	{
		if ((pEngine->m_IterationsDone > 0) && (pEngine->m_IterationsDone % 1000) == 0) {
			std::cout << pEngine->m_IterationsDone << ' ' << std::scientific << pEngine->FER() << std::endl;
		}
		return true; // do nothing, continue simulations
	}
};

#endif