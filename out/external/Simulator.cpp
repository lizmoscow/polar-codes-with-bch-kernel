#include <math.h>
#include <iostream>
#include <fstream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#define NOMINMAX
#include <Windows.h>
#include <signal.h>
#include "Simulator.h"
#include "misc.h"

//#define _PRINT_ERROR_IT
//#define PRINT_FIRST_BIT_ERROR
using namespace std;
bool CSimulator::KeepRunning = true;
void ControlC_Handler(int dummy)
{
    CSimulator::KeepRunning = false;
}

//#define TRANSMIT_RANDOM_SYMBOLS
//#define TRANSMIT_PURE_NOISE

//#define CHECK_ENCODE_AND_EXTRACT // no noise, just encode data and extract it from the codeword

/// @param the codec to be used
/// @param
/// @param 
CSimulator::CSimulator(CBinaryCodec& Codec, CModem& Modem, eChannel ChannelType) :
m_Modem(Modem),
m_ChannelType(ChannelType),
m_NumOfBitErrors(0),
m_NumOfBlockErrors(0),
m_NumOfMLErrors(0),
m_NumOfRawBitErrors(0),
m_NumOfCodeBitErrors(0),
m_NumOfUndetectedError(0),
m_NumOfListErrors(0),
m_IterationsDone(0),
m_DecodingTime(0),
m_MemoryUsage(0), m_MemoryUsageSuccess(0),
m_PrintErrorIterations(
#ifdef _DEBUG
	//true
	false
#else
	false
#endif
),
m_ZeroData(0), m_OneData(false),
m_AssertNonMLError(false),
m_NoBER(false),
m_PrintDecodedData(
#ifdef _DEBUG
	//true
#else
	false
#endif
), 
m_PrintCodeword(
#ifdef _DEBUG
	//true
#else
	false
#endif
)
{
    m_pEncoder = dynamic_cast<CBinaryEncoder*>(&Codec);
    m_pDecoder = dynamic_cast<CBinarySoftDecoder*>(&Codec);
    m_NumOfSymbols = Codec.GetLength() / Modem.GetBitsPerSymbol();
    if (Codec.GetLength() % Modem.GetBitsPerSymbol())
        m_NumOfSymbols++;

    m_DecrementForNumOfSymbols = 0;
    m_pData = new tBit[Codec.GetDimension()];
    m_pDecodedData = new tBit[Codec.GetDimension()*m_pDecoder->GetMaxListSize()];
    m_pEncoded = new tBit[m_NumOfSymbols*Modem.GetBitsPerSymbol()];
    m_pDecodedCW = new tBit[Codec.GetLength()*m_pDecoder->GetMaxListSize()];
    m_pNoisy = new Complex[m_NumOfSymbols];
    m_pModulated = new Complex[m_NumOfSymbols];
	m_pChannel = new Complex[m_NumOfSymbols];
    memset(m_pEncoded + Codec.GetLength(), 0, sizeof(tBit)*(m_NumOfSymbols*Modem.GetBitsPerSymbol() - Codec.GetLength()));
    m_pLLRs = (MType*)_aligned_malloc(sizeof(MType)*m_NumOfSymbols*Modem.GetBitsPerSymbol(), 32);
//  gsl_rng_default_seed = time(nullptr);
    m_pRNG = gsl_rng_alloc(gsl_rng_mt19937);
    signal(SIGINT, ControlC_Handler);
}


CSimulator::~CSimulator()
{
    delete[]m_pData;
    delete[]m_pEncoded;
    delete[]m_pDecodedCW;
    delete[]m_pDecodedData;
    delete[]m_pNoisy;
    delete[]m_pModulated;
	delete[]m_pChannel;
    _aligned_free(m_pLLRs);
    gsl_rng_free(m_pRNG);
}


void CSimulator::SetEbN0(double SNR)
{
	if (m_ChannelType == ec_AWGN)
	{
		m_NoiseStdDev = sqrt(0.5*pow(10, -SNR / 10)*(m_NumOfSymbols - m_DecrementForNumOfSymbols) / m_pEncoder->GetDimension());
		m_Modem.SetNoiseVariance(m_NoiseStdDev * m_NoiseStdDev);
	}
	else
	{
		//m_ErasureProbability = SNR;
		m_Modem.SetNoiseVariance(1);
	}
}

void CSimulator::SetEsN0(double SNR)
{
	if (m_ChannelType == ec_AWGN)
	{
		m_NoiseStdDev = sqrt(pow(10, -SNR / 10));
		m_Modem.SetNoiseVariance(m_NoiseStdDev * m_NoiseStdDev);
	}
	else
	{
		//m_ErasureProbability = SNR;
		m_Modem.SetNoiseVariance(1);
	}
}

void CSimulator::SetPrecomputedNoiseStdDev(double NoiseStdDev)
{
    m_NoiseStdDev = NoiseStdDev;
    m_Modem.SetNoiseVariance(m_NoiseStdDev * m_NoiseStdDev);
}


bool CSimulator::Iterate(unsigned long long It, bool Indeed, const tBit* pData)
{
    unsigned K = m_pEncoder->GetDimension();
    unsigned N = m_pEncoder->GetLength();

#if !defined(TRANSMIT_RANDOM_SYMBOLS) && !defined(TRANSMIT_PURE_NOISE)
    if (m_ZeroData) {
        pData = m_pData;
        memset(m_pData, 0, sizeof(tBit)*K);
    } else if (m_OneData) {
        pData = m_pData;
        memset(m_pData, 1, sizeof(tBit)*K);
    } else {
        if (!pData) {
            pData = m_pData;
            for (unsigned i = 0; i < K; i++)
                m_pData[i] = static_cast<tBit>(gsl_rng_uniform_int(m_pRNG, 2));
        }
    }
	if (Indeed)
        m_pEncoder->Encode(pData, m_pEncoded);
#elif defined(TRANSMIT_RANDOM_SYMBOLS)
    for (unsigned i = 0; i < N; ++i)
        m_pEncoded[i] = static_cast<tBit>(gsl_rng_uniform_int(m_pRNG, 2));
#else
    memset(m_pEncoded, 0, sizeof(tBit)*N);
#endif
    Complex Channel = 1.;
    Complex* pNoisy = m_pNoisy;
    Complex* pModulated = m_pModulated;
	Complex* pChannel = m_pChannel;
    double EW = 0;
    for (unsigned i = 0; i < N; i += m_Modem.GetBitsPerSymbol()) {
        double Noise = gsl_ran_gaussian(m_pRNG, m_NoiseStdDev);
#ifndef TRANSMIT_PURE_NOISE
        Complex Modulated = (Indeed) ? m_Modem.Modulate(m_pEncoded + i) : 0;
#else
        Complex Modulated = 0;
#endif
        Complex Noisy; double h;
        switch (m_ChannelType) {
        case ec_AWGN:
            Noisy = Modulated.real() + Noise;
            break;
        case ec_RayleighReal:
            h = gsl_ran_rayleigh(m_pRNG, M_SQRT1_2);
            Noisy = h*Modulated.real() + Noise;
            Channel = h;
            break;
		case ec_BSC:
			h = gsl_ran_flat(m_pRNG, 0., 1.);
			Noisy = h > m_NoiseStdDev ? Modulated : -Modulated;
			break;
		case ec_BEC:
			//Noisy = (gsl_rng_uniform(m_pRNG) < m_ErasureProbability) ? 0 : Modulated;
			break;
        default:
            throw Exception("Unknown channel type");
        }
        *(pNoisy++) = Noisy;
        *(pModulated++) = Modulated;
		*(pChannel++) = Channel;
        if (Indeed) {
            m_Modem.Demodulate(Noisy, Channel, m_pLLRs + i);
            for (unsigned j = 0; j < m_Modem.GetBitsPerSymbol(); j++) {
                double Z = (m_pEncoded[i + j]) ? umin(m_pLLRs[i + j]) : m_pLLRs[i + j];
                if (Z < 0) {
                    m_NumOfRawBitErrors++;
                    EW += abs_value(m_pLLRs[i + j]);
                }
            }
        }
    }
    if (!Indeed)
        return true;
    LARGE_INTEGER OldClock, NewClock;

    if (m_PrintDecodedData) fill(m_pDecodedData, m_pDecodedData + m_pDecoder->GetDimension(), (~0));
    if (m_PrintErrorIterations) fill(m_pDecodedCW, m_pDecodedCW + m_pDecoder->GetLength(), (~0));

    QueryPerformanceCounter(&OldClock);
    int L = m_pDecoder->Decode(m_pLLRs, (m_NoBER) ? 0 : m_pDecodedData, m_pDecodedCW);
    QueryPerformanceCounter(&NewClock);
    m_DecodingTime += NewClock.QuadPart - OldClock.QuadPart;

#ifdef CHECK_ENCODE_AND_EXTRACT
    L = 1;
    memcpy(m_pDecodedCW, m_pEncoded, m_pDecoder->GetLength()*sizeof(*m_pDecodedCW));
    m_pDecoder->ExtractInformationBits(m_pDecodedCW, m_pDecodedData);
#endif

    bool HasError = false;
    unsigned CurBitErrors = 0;
    unsigned CurCodeBitErrors = 0;
    if (L <= 0) HasError = true;
    if ((L <= 0)&&m_NoBER) {
        HasError = true;
    } else {
        if (m_NoBER) {
            for (unsigned i = 0; i < N; ++i) {
                if ((m_pDecodedCW[i] != 0) != (m_pEncoded[i] != 0)) {
                    HasError = true;
                    break;
                }
            }
        } else {
            for (unsigned i = 0; i < K; i++) {
                if ((m_pDecodedData[i] != 0) ^ (pData[i] != 0)) {
                    CurBitErrors++;
                }
            }
            for (unsigned i = 0; i < N; i++) {
                if ((m_pDecodedCW[i] != 0) ^ (m_pEncoded[i] != 0)) {
                    CurCodeBitErrors++;
                }
            }
            HasError |= (CurBitErrors > 0);
        }
        if (HasError) {
            if (L)
                m_NumOfUndetectedError++;
        }
    }
    //cout << "MLCorr=" << MLCorr << "    CurCorr=" << CurCorr << endl;
	bool ListError = HasError;
    if (HasError) {
        double MLCorr = 0;
        double CurCorr = 0;
        if (m_Modem.GetBitsPerSymbol() == 1) {
            for (unsigned i = 0; i<N; i++) {
                MLCorr += (m_pEncoded[i]) ? -m_pLLRs[i] : m_pLLRs[i];
                CurCorr += (m_pDecodedCW[i]) ? -m_pLLRs[i] : m_pLLRs[i];
            }
        } else {
            //modulate the codeword, and compare the euclidean distance
            for (unsigned i = 0; i < m_NumOfSymbols; i++) {
                Complex ReModulated = m_Modem.Modulate(m_pDecodedCW + i*m_Modem.GetBitsPerSymbol());
                MLCorr -= norm(m_pModulated[i] * m_pChannel[i] - m_pNoisy[i]);
                CurCorr -= norm(ReModulated * m_pChannel[i] - m_pNoisy[i]);
            }
        }

        if (m_PrintErrorIterations)
            cerr << '!' << It << '\n';
        m_NumOfBitErrors += CurBitErrors;
        m_NumOfCodeBitErrors += CurCodeBitErrors;
        m_NumOfBlockErrors++;
        if ((CurCorr >= MLCorr)&&L)
            m_NumOfMLErrors++;
		for (int j = 1; j < L; j++)
		{
			//check if the correct codeword is in the list
			bool CurError = false;
			for (unsigned i = 0; i < N; ++i) 
			{
				if ((m_pDecodedCW[j*N+i] != 0) != (m_pEncoded[i] != 0)) 
				{
					CurError = true;
					break;
				}
			};
			if (!CurError)
			{
				ListError = false;
				break;
			}
		}
		if (m_AssertNonMLError)
            assert(false);

        if (m_AssertNonMLError) assert(false);

    } 
	else
	{
        if (CNT_INFO(CNT_MEMORY_CS) > m_MemoryUsageSuccess)
            m_MemoryUsageSuccess = CNT_INFO(CNT_MEMORY_CS);
    }
	if (ListError)
		m_NumOfListErrors++;
	if (m_PrintDecodedData)
	{
		cout << "Tx/Rx data:\n";
		printData(pData, m_pDecodedData, K);
	}
	if (m_PrintCodeword)
	{
		cout << "Tx/Rx codeword:\n";
		printData(m_pEncoded, m_pDecodedCW, N);
	}
    if (m_PrintDecodedData || m_PrintCodeword) cin.get();

    if (CNT_INFO(CNT_MEMORY_CS) > m_MemoryUsage)
        m_MemoryUsage = CNT_INFO(CNT_MEMORY_CS);
    m_MemoryUsageStatic = CNT_INFO(CNT_MEMORY_OUTER);
    return !HasError;
}


/// print two tBit vectors to compare them visually
void CSimulator::printData(const tBit* pData0, const tBit* pData1, unsigned N)
{
    for (unsigned i = 0; i < N; ++i)
        cout << ((pData0[i] == BIT_0 || pData0[i] == 0) ? '0' : ((pData0[i] == BIT_1 || pData0[i] == 1) ? '1' : '?'));
    cout << endl;
    for (unsigned i = 0; i < N; ++i)
        cout << ((pData1[i] == BIT_0 || pData1[i] == 0) ? '0' : ((pData1[i] == BIT_1 || pData1[i] == 1) ? '1' : '?'));
    cout << endl;
    bool match = true;
    for (unsigned i = 0; i < N; ++i) {
        if ((pData0[i] == BIT_0) != (pData1[i] == BIT_0)) {
            match = false;
            break;
        }
    }
    cout << "equals? " << ((match) ? "YES" : "NO") << endl;
}


/// Do real encoding and decoding only for the specified iteration.
/// @param iteration number
/// @param some specific data pattern
void CSimulator::RunOneIteration(unsigned IterationID, const tBit* pData)
{
    for (m_IterationsDone; m_IterationsDone<=IterationID;) {
        if (m_IterationsDone!=IterationID)
            Iterate(m_IterationsDone++, false, pData);
        else {
            Iterate(m_IterationsDone++, true, pData);
            break;
        }
    }
}


double CSimulator::GetAvgDecodingTime()const
{
    LARGE_INTEGER Freq;
    QueryPerformanceFrequency(&Freq);
    return double(m_DecodingTime) /(Freq.QuadPart*m_IterationsDone);
}


void CSimulator::SetRngSeed(unsigned long Seed /*= std::chrono::system_clock::now().time_since_epoch().count()*/)
{
    gsl_rng_set(m_pRNG, Seed);
}


void CSimulator::LengthDecrement(unsigned LengthForEbN0)
{
    int l = LengthForEbN0 / m_Modem.GetBitsPerSymbol();
    if (LengthForEbN0 % m_Modem.GetBitsPerSymbol())
        l++;
    m_DecrementForNumOfSymbols = (int)m_NumOfSymbols-l;
}
