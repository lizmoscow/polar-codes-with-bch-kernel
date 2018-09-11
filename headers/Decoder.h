//
// Created by Moskovskaya Elizaveta on 12.05.18.
//

#ifndef BCH_DECODER_H
#define BCH_DECODER_H


class Decoder {

private:

    struct matrix {
        unsigned long *ff;
        unsigned long *fs;
        unsigned long *sf;
        unsigned long *ss;
        int sizeFF;
        int sizeFS;
        int sizeSF;
        int sizeSS;
    };

    long power;
    long n;
    long t;
    long k;
    unsigned long *antilogarithms;
    unsigned long *logarithms;

    unsigned long *lambda;
    int sizeLambda;
    unsigned long *locators;
    unsigned char *oldPoly;


    long l;

    unsigned long *s; int sizeS;
    unsigned long *p; int sizeP;
    unsigned long *q; int sizeQ;
    unsigned long *temp; int sizeTemp;
    matrix a;
    matrix tempA;

    unsigned long *f;
    unsigned long *shift;
    unsigned long *res1;

    inline unsigned long multiply(unsigned long first, unsigned long second) const;
    inline unsigned long divide(unsigned long first, unsigned long second) const;
    void addPolynomials(const unsigned long* first, int size1, const unsigned long* second, int size2,
                                 unsigned long *result, int *sizeRes) const;
    void multiplyPolynomials(const unsigned long* first, int size1, const unsigned long* second, int size2,
                                      unsigned long *result, int *size) const;
    void dividePolynomial(unsigned long* first, int *size1, unsigned long* second, int *size2,
                                   unsigned long *result, int *sizeResult) const;
    unsigned long eval(const unsigned long *poly, int size, unsigned long elem) const;


    bool euclid(const unsigned char *word, unsigned long *lambda, int *size);
    bool locatorsAndRoots(const unsigned long *lambda, int size, unsigned long *locators) const;


public:

    Decoder(long pw, long n, long t, long k, unsigned long *antilogarithms, unsigned long *logarithms);
    ~Decoder();
    void findSyndromPoly(const unsigned char *word);
    void alterSyndromPoly(const unsigned char *word);
    bool decode(const unsigned char* word, unsigned char *answer);
    long getN() const;
    long getT() const;
    long getK() const;


    unsigned long *syndromPoly;
    long syndromPolySize;
};


#endif //BCH_DECODER_H
