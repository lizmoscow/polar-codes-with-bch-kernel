#ifndef DYNFROZENX_H
#define DYNFROZENX_H
#include <fstream>

#include "Codec.h"

enum SymbolType {stNormal, stShortened, stPunctured};

typedef std::pair<unsigned, unsigned> tPairOfUnsigned;

#define FROZENDISTANCESCHEDULER

//Arikan dynamic frozen encoder with extension according to X construction
class CDynFrozenEncoderX :public CBinaryEncoder
{
protected:
	//for each global phase we need to store the corresponding engine ID and local phase
	tPairOfUnsigned* m_pEngineMapping;

public:
	//number of polarizing transformations
	unsigned m_NumOfPolarizers;
	///number of layers in each Arikan transformation
	unsigned* m_pNumOfArikanLayers;
	///IDs of the freezing constraints, 
    ///m_ppFreezingConstraints[m_pDecisionPoints[phi]] - constraints on phase phi
	unsigned* m_pDecisionPoints;
    ///Array of freezing constraints: [p0 p1 ... phi], where phi is the constraint's phase
	unsigned** m_ppFreezingConstraints;
	///encoding buffer
	tBit* m_pEncodingBuffer;
	///minimum distance of the code
	unsigned m_MinDist;
	//number of shortened symbols
	unsigned m_NumOfShortened;
    //number of punctured symbols
    unsigned m_NumOfPunctured;
	//types of symbols
	SymbolType* m_pSymbolTypes;
    //length of unshortened and unpunctured code
    unsigned m_ArikanLength;

	CDynFrozenEncoderX() {}
	CDynFrozenEncoderX(std::istream& SpecFile);
	virtual ~CDynFrozenEncoderX();
	virtual void Encode(const tBit* pSrc, tBit* pEncoded);
	virtual void ExtractInformationBits(const tBit* pCodeword, tBit* pInfBits);
	void ExtractInformationBitsUnshortened(const tBit* pCodeword, tBit* pInfBits);
	virtual unsigned GetMinimumDistance()const
	{
		return m_MinDist;
	}

	unsigned GetNumOfPolarizers() const
	{
		return m_NumOfPolarizers;
	}

	unsigned GetNumOfArikanLayers(unsigned TransformNum) const {
		return m_pNumOfArikanLayers[TransformNum];
	}

};


#endif