//
// Created by Moskovskaya Elizaveta on 10/06/2018.
//

#ifndef BCH_KANEKO_H
#define BCH_KANEKO_H


#include <utility>
#include "Decoder.h"
//#include "external/KernProc.h"

class KanekoKernelProcessor /*: CKernProcLLR*/ {

private:
    long n;
    long t;
    double sd;

    Decoder decoder;

    double *alpha;
    std::pair<double, int> *alphaSorted;
    unsigned char *yH;
    unsigned char *x;
    unsigned char *err;
    long m;
    long m0;
    double state;

    unsigned long comparisonCount{0};
    unsigned long summCount{0};
    unsigned long decodingCount{0};

    void calcError(long i);
    long calcM();
    double calcL() const;
    double calcRightSide() const;
    double calcT(long j);

public:

    KanekoKernelProcessor(long pw, long n, long t, long k, unsigned long *antilogarithms, unsigned long *logarithms, double signalToNoiseRatio);
    virtual ~KanekoKernelProcessor();

    void decode(unsigned char *res);
    void decode(const double *word, unsigned char *res);
    void set(const double *word) const;

    long getN() const;
    long getT() const;
    long getK() const;
    double calcL(const unsigned char *word) const;

    unsigned long getComparisonCount() const;
    unsigned long getSummCount() const;
    unsigned long getDecodingCount() const;
    void setDecodingCount(unsigned long c = 0);
    void setComparisonCount(unsigned long comparisonCount = 0);
    void setSummCount(unsigned long summCount = 0);
    //virtual void* GetState(unsigned Stride)const;

    virtual const char* GetName()const { return "KanekoProcessor"; }



};


#endif //BCH_KANEKO_H
