#ifndef ARIKAN_CODEC_H
#define ARIKAN_CODEC_H

#include "ArikanEncoder.h"

class CArikanCodec:public CArikanEncoder
{
    //LLRs for each layer
    double** m_ppLLRs;
    //decisions for each layer
    tPairOfBool** m_ppC;
    tBit* m_pDecoderBuffer;
    void RecursivelyCalcLLR(unsigned lambda,unsigned phi);
    void RecursivelyUpdateC(unsigned lambda,unsigned phi);

public:
    CArikanCodec(unsigned m,///log_2(code length)
                 unsigned Dimension,//k
                 const unsigned* pNonFrozenChannels //nonfrozen channel IDs sorted in ascending order
                 );
    ~CArikanCodec();
    //soft-input decoder. Produces  an information vector
    //returns 1 on success
    virtual  int Decode(const double* pLLR,//log(P{c_i=1|y_i}/P{c_i=0|y_i})
                       tBit* pInfVectorList //output array. Must be sufficiently big to accomodate GetMaxListSize() vectors of length GetDimension()
                       );
    virtual unsigned GetMinimumDistance()const
    {
        return 0;//TBD
    };

};
//use approximate implementation of -log tanh(x/2)
#define GALLAGER_APPROXIMATION

///Gallager's function used by the decoder
double Gallager(const double& x);

#endif