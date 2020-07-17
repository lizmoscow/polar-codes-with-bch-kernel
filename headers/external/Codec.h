#ifndef CODEC_H
#define CODEC_H
#include <string.h>
#include "misc.h"
#include "Simulation/ChannelTemplates.h"

//disabv
#pragma warning(disable:4250)

class DLL_API CBlockCodec
{
protected:
    //length of the code
    unsigned m_Length;
    //dimension of the code
    unsigned m_Dimension;

	CBlockCodec(unsigned Length,unsigned Dimension):m_Length(Length),m_Dimension(Dimension)
    {
    };
	CBlockCodec():m_Length(0),m_Dimension(0)
    {

    };

public:
    virtual ~CBlockCodec()
    {
    };
    virtual unsigned GetLength()const
    {
        return m_Length;
    };
    virtual unsigned GetDimension()const
    {
        return m_Dimension;
    };
    //get the code minimum distance. Return lower bound
	virtual unsigned GetMinimumDistance()const { return 0; }

    //constructs a generator matrix in packed form. Each row has size  ceil(GetLength()/(8*sizeof(Word)) 
	virtual Word* GetGeneratorMatrix() { return nullptr; }
    //constructs check matrix in packed column-wise format
    //each column will have size ceil((m_CodeLength-m_CodeDimension)/(8*sizeof(Word)) 
    Word* GetCheckMatrix();
};

typedef CBlockCodec CBinaryCodec;


//codec for a binary code
class DLL_API CBinaryEncoder :public virtual CBinaryCodec
{
protected:
	CBinaryEncoder() :CBinaryCodec(0, 0)
	{
	};
public:
    virtual ~CBinaryEncoder()
    {
    };
	CBinaryEncoder(unsigned Length, unsigned Dimension) :CBinaryCodec(Length, Dimension)
	{
	};

    //construct a generator matrix in packed form. Each row will have size ceil(GetLength/(8*sizeof(Word)) if GetExtension()==1
    Word* GetGeneratorMatrix();
    //encode a block of data into a codeword
    virtual void Encode(const tBit* pSrc,tBit* pEncoded)=0;
	//extract information bits from a codeword
	virtual void ExtractInformationBits(const tBit* pCodeword, tBit* pInfBits) = 0;
};

//codec for a binary code given by a generator matrix
class DLL_API CBinaryLinearEncoder :public CBinaryEncoder
{
protected:
    //column permutation used to obtain the canonical form
    unsigned* m_pInformationSet;
    //generator matrix in the canonical form. It may not be the same as given to the constructor
    tBit* m_pGeneratorMatrix;
	//linear transformation for message extraction
	Word* m_pMessageExtractionMatrix;
    //minimum distance, if known
    unsigned m_MinDist;
public:
	CBinaryLinearEncoder(unsigned Length, unsigned Dimension, const tBit* pGenMatrix, unsigned MinDist = 0);
    ~CBinaryLinearEncoder();
    void Encode(const tBit* pSrc,tBit* pEncoded);
    void ExtractInformationBits(const tBit* pCodeword,tBit* pInfBits);
    unsigned GetMinimumDistance()const
    {
        return m_MinDist;
    };

};

#pragma vtordisp(push, 2) 
//a generic soft-input list-output decoder
class DLL_API CBinarySoftDecoder :public virtual CBinaryCodec
{
public:
    virtual unsigned GetMaxListSize()const
    {
        return 1;
    };
    //soft-input decoder. Produces  a list of information vectors
    //returns the number of vectors obtained, or -1 in case of error
    virtual int Decode(const MType* pLLR,
                       //log(P{c_i=1|y_i}/P{c_i=0|y_i})
                       tBit* pInfVectorList,
                       //output array. Must be sufficiently big to accomodate GetMaxListSize() vectors of length GetDimension()
                       tBit* pCodewordList=0 //output array. If non-zero, must be sufficiently big to accomodate GetMaxListSize() vectors of length GetLength()
    )=0;

	///set noise standard deviation for the decoder
	virtual void SetAWGNStdDev(double Sigma) {}
};
#pragma vtordisp(pop)  

/**
 * \brief Codec which represents a modulated code
 */
class DLL_API CModulatedCodec : public CBlockCodec
{
protected:
	/// number of bits per one channel symbol
	unsigned m_BitsPerSymbol;
	CModulatedCodec() : CBlockCodec(0, 0), m_BitsPerSymbol(0) {}
public:

	/**
	 * \brief Construct the modulated codec
	 * \param Length code length in channel symbols
	 * \param Dimension code dimension in bits
	 * \param BitsPerSymbol number of bits per one channel symbol
	 */
	CModulatedCodec(unsigned Length, unsigned Dimension, unsigned BitsPerSymbol)
	: CBlockCodec(Length, Dimension), m_BitsPerSymbol(BitsPerSymbol) {}

	/**
	 * \brief Get number of bits per one channel symbol
	 * \return bits per channel symbol
	 */
	unsigned GetBitsPerSymbol() const
	{
		return m_BitsPerSymbol;
	}
};

class DLL_API CModulatedEncoder : public virtual CModulatedCodec
{
protected:
	CModulatedEncoder() : CModulatedCodec(0, 0, 0) {}
public:
	CModulatedEncoder(unsigned Length, unsigned Dimension, unsigned BitsPerSymbol)
	: CModulatedCodec(Length, Dimension, BitsPerSymbol) {}

	// encode a block of data into a modulated codeword
	virtual void Encode(const tBit* pSrc, tChannelSymbol *pEncoded) = 0;
};

namespace nsSimulation {
	struct CChannelDataInterface;
};

class DLL_API CModulatedDecoder : public virtual CModulatedCodec
{
protected:
	CModulatedDecoder() : CModulatedCodec(0, 0, 0) {}
public:
	/**
	 * \brief Modulated decoder constructor
	 * \param Length Code length in channel symbols
	 * \param Dimension Code dimension in bits
	 * \param BitsPerSymbol Number of bits per one channel symbol
	 */
	CModulatedDecoder(unsigned Length, unsigned Dimension, unsigned BitsPerSymbol)
		: CModulatedCodec(Length, Dimension, BitsPerSymbol) {}

	virtual unsigned GetMaxListSize() const { return 1;}

	/**
	 * \brief Soft-input modulated decoder. Produces a list of information vectors
	 * \param pChannelData   channel output with noised symbols and additional data
	 * \param pInfVectorList output array. Must be sufficiently big to accomodate GetMaxListSize() vectors of length GetDimension()
	 * \param pCodewordList output array. If non-zero, must be sufficiently big to accomodate GetMaxListSize() vectors of length GetLength()
	 * \return the number of vectors obtained, or -1 in case of error
	 */
	virtual int Decode(const nsSimulation::CChannelDataInterface *pChannelData, tBit *pInfVectorList, tChannelSymbol *pCodewordList = nullptr) = 0;

	
	/**
	 * \brief Adjust the decoder to the current channel
	 * \param Channel current channel
	 */
	virtual void SetChannel(nsSimulation::CChannelInterface &Channel) {}
};
#endif
