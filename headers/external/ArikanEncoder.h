#ifndef ARIKANENCODER_H
#define ARIKANENCODER_H
//encoder for Arikan polar code
#include "Codec.h"
#include "SeqConfig.h"
#include <algorithm>
#include "LinAlg.h"

/// get position of a particular element in in an ordered vector.
int VectorFind(const unsigned * pStart,/*!<sorted vector to be searched*/
                unsigned VectorSize,/*!<Size of the vector*/
                unsigned Element/*!<Element to be found*/
                );



class CArikanEncoder:virtual public CBinaryEncoder
{
    //true if systematic encoding is needed
    bool m_bSystematic;
	//temp codeword for nonsystematic encoding
	tBit *m_pTempCodeword;
protected:
    //true if the i-th channel is frozen
    bool* m_pFrozen;
    //log_2(code length)
    unsigned m_LogCodeLength;
    //bit reversal permutation
    unsigned* m_pBitReversal;
    //sorted sequence of non-frozen channels
    const unsigned* m_pNonFrozenChannels;
public:
    //initialize polar code (2^m,k) with a given list of frozen channels
    CArikanEncoder(unsigned m,///log_2(code length)
                    unsigned Dimension,//k
                    const unsigned* pNonFrozenChannels, //nonfrozen channel IDs sorted in ascending order
                    bool Systematic //true if systematic encoding is needed
                    );
    virtual ~CArikanEncoder();
    //encode a block of data into a codeword
    void Encode(const tBit* pSrc,tBit* pEncoded);
	void ExtractInformationBits(const tBit* pCodeword, tBit* pInfBits);
    unsigned GetMinimumDistance()const
    {
        unsigned x=~0u;
        for(unsigned i=0;i<m_Dimension;i++)
        {
            unsigned w=_mm_popcnt_u32(m_pNonFrozenChannels[i]);
            if (w<x) x=w;
        };
        return 1u<<x;
    };
	bool IsFrozen(unsigned ChannelID)const {
		return m_pFrozen[ChannelID];
	};

};


//Gaussian approximation for density evolution
void GetSubchannelQuality(unsigned m,//log_2(code length)
                          double sigma,//AWGN standard deviation
                          double* pErrProb ///estimated error probability
                          );
//Gaussian approximation for density evolution for extended polar codes
void GetSubchannelQualityExtended(unsigned m //number of layers in the polarizing transformation
	, unsigned Length//actual length of the code
	, double sigma//AWGN standard deviation
	, double* pErrProb ///estimated error probability
	);
//Gaussian approximation for density evolution
void GetSubchannelQualityShortened(unsigned m //number of layers in the polarizing transformation
	, unsigned Length//actual length of the code
	, double sigma//AWGN standard deviation
	, double* pErrProb ///estimated error probability
	);
//Gaussian approximation for density evolution
void GetSubchannelQualityPunctured(unsigned m //number of layers in the polarizing transformation
	, unsigned Length//actual length of the code
	, double sigma//AWGN standard deviation
	, double* pErrProb ///estimated error probability
	);

//Gaussian approximation for density evolution. 
void GetSubchannelQualityGAPrim(unsigned m //number of layers in the polarizing transformation
	, double* pLLRMean ///on input the mean values of each codeword symbol LLR. On output the error probabilities for the input symbols of the polarizing transformation
	);

void GetSubchannelQualityWithMax(unsigned m,//log_2(code length)
                          double sigma,//AWGN standard deviation
                          double* pErrProb ///estimated error probability
                          );

//subchannel qualities (the less the better) given by Polarizing Weights method
void GetSubchannelQualityPolarizingWeight(unsigned m, double* pW);

//find a number of good channels in 2^m-dimensional Arikan transformation
unsigned* GetGoodChannels(unsigned m,//log_2(code length)
                          double sigma,//AWGN standard deviation
                          unsigned NumOfGoodChannels, //number of good channels to be constructed
                          bool BEC, //if true, sigma is the erasure probability
                          unsigned TargetMinDist=0 //exclude bit subchannels leading to minimum distance <TargetMinDist
                          );
void GetSubchannelQualityBEC(unsigned m,//log_2(code length)
                          double ErasureProb,//erasure probability
                          double* pErrProb ///estimated error probability
                          );
void GetSubchannelQualityBECShortened(unsigned m, double* pErasureProb);

void GetPolarChannelNoise(unsigned m,//log_2(code length)
                          double sigma,//AWGN standard deviation
                          double* pNoiseStdDev ///noise standard deviation for each subchannel
                          );
//Gaussian approximation for density evolution
void GetSubchannelLLRMean(unsigned m,//log_2(code length)
	double sigma,//AWGN standard deviation
	double* pLLRMean ///estimated LLR mean, sigma^2=2*LLRMean
	);
//channel selection according to the balanced criterion
unsigned* MakeGoodChannelsBalanced(unsigned LogLength,unsigned Dimension);
unsigned* GetBitReversal(unsigned m);

double GetListErasureDecodingErrProb(unsigned m,//code of length 2^m will be considered
                         unsigned Dimension, //number of non-frozen symbols
                         unsigned MaxListDimension,///maximal dimension of the list space
                         const unsigned* pNonFrozenSymbols,//frozen subchannels (must be sorted)
                         double ChannelErasureProb //channel erasure probability
                         );
unsigned* OptimizeCodeBEC(unsigned m,unsigned TargetDimension,unsigned MaxPathSpaceDim,double ChannelErasureProb);
void PolarSysEncode(tBit* pSrc,unsigned InputSize,const bool* pFrozenChannels);


//inplace calculation of xF^{\otimes m}, where F is the Arikan kernel
template <class T> void Arikan(unsigned LogLength, T* pData)
{
	unsigned N = 1u << LogLength;
	//implement Arikan's butterfly
	unsigned L = 1;
	while (N>1)
	{
		N >>= 1;
		for (unsigned k = 0; k<L; k++)
		{
			/*for (unsigned j = 2 * k*N; j<(2 * k + 1)*N; j++)
			{
				pData[j] ^= pData[j + N];
			};*/
			XOR(pData + 2 * k*N, pData + 2 * k*N+N, N);
		};
		L <<= 1;
	};
};

//inplace calculation of xF^{\otimes m}, where F is the Arikan kernel
template <class T> void ArikanWithPermutationInsideLayer(unsigned LogLength, unsigned layerID, unsigned** ppOuterPermutation, T* pData, T* pTemp)
{
	unsigned N = 1/*1u << LogLength*/;
	//implement Arikan's butterfly
	unsigned L = 1u << LogLength;
	unsigned id = 1;
	while (L>1)
	{
		L >>= 1;
		for (unsigned k = 0; k<L; k++)
		{
			for (unsigned j = 2 * k*N; j<(2 * k + 1)*N; j++)
			{
				pData[j] ^= pData[j + N];
			};
		};
		N <<= 1;

		// for puncturing
		if (id == layerID) {
			//unsigned *pRB = GetBitReversal(layerID);
			memcpy(pTemp, pData, sizeof(tBit)*(1u << LogLength));
			// perform permutation for symbols of outer codes codewords
			for (unsigned k = 0; k < L; k++) {
				for (unsigned j = 0; j < N; j++) pData[k * N + j] = pTemp[/*pRB[ppOuterPermutation[k][pRB[j]]] * N + j*//*k * N + pRB[ppOuterPermutation[k][pRB[j]]]*/k * N + ppOuterPermutation[k][j]];
			}
		}
		id++;
	};
	/*unsigned *pRB = GetBitReversal(LogLength);
	memcpy(pTemp, pData, sizeof(tBit)*(1u << LogLength));
	for (unsigned k = 0; k < (1u << LogLength); k++) pData[k] = pTemp[k];*/
};

//inplace calculation of xF^{\otimes m}, where F is the Arikan kernel
template <class T> void ArikanWithInversePermutationInsideLayer(unsigned LogLength, unsigned layerID, unsigned** ppOuterPermutation, T* pData, T* pTemp)
{
	unsigned N = 1u << LogLength;
	//implement Arikan's butterfly
	unsigned L = 1;
	unsigned id = 1;
	while (N>1)
	{
		N >>= 1;
		for (unsigned k = 0; k<L; k++)
		{
			for (unsigned j = 2 * k*N; j<(2 * k + 1)*N; j++)
			{
				pData[j] ^= pData[j + N];
			};
		};
		L <<= 1;

		// for puncturing
		if (LogLength-id == layerID) {
			memcpy(pTemp, pData, sizeof(tBit)*(1u << LogLength));
			// perform permutation for symbols of outer codes codewords
			for (unsigned k = 0; k < L; k++) {
				for (unsigned j = 0; j < N; j++) pData[ppOuterPermutation[k][j] + k*N] = pTemp[j + k*N];
			}
		}
		id++;
	};
};

//inplace calculation of x(F^{\otimes m})^T, where F is the Arikan kernel
template <class T>void ArikanTransposed(unsigned LogLength, T* pVector)
{

	unsigned N = 1u << LogLength;
	unsigned l = 1;
	while (N>0)
	{
		N >>= 1;
		for (unsigned j = 0; j<l; j++)
		{
			for (unsigned i = 0; i<N; i++)
				pVector[j * 2 * N + i + N] ^= pVector[j * 2 * N + i];
		};
		l += l;
	};

};

//compute LLR mean for x+y mod 2 given LLR mean for x,y
double XORDensity(double SrcLLRMean);
double XORDensity(double SrcLLRMean1, double SrcLLRMean2);
//compute LLR mean for x+y  given LLR mean for x,y
double SumDensity(double SrcLLRMean);
double SumDensity(double SrcLLRMean1, double SrcLLRMean2);

double ChungPhi(double x);
double ChungPhiInverse(double y);
#endif