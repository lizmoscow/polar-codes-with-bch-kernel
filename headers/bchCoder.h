//
// Created by Moskovskaya Elizaveta on 11.05.18.
//

#ifndef BCH_BCHCODER_H
#define BCH_BCHCODER_H

#include <algorithm>

void findMinimalPolynomial(int i, int power, const unsigned long * fieldElements, int *size, unsigned char* res);

bool comparePoly(const unsigned char* poly1, int size1, const unsigned char* poly2, int size2);

unsigned char* multiplyPolynomials(const unsigned char* first, int size1, const unsigned char* second, int size2, int *sizeRes = nullptr);

void multiplyPolynomials(const unsigned char* first, int size1, const unsigned char* second, int size2, unsigned char* res, int *sizeRes = nullptr);

unsigned char* dividePolynomial(const unsigned char* first, int size1, const unsigned char* second, int size2, int *size, bool needRemainder);

unsigned char* lcm(const unsigned char* first, int size1, const unsigned char* second, int size2, int *sizeRes);

inline unsigned long getPower(unsigned __int128 poly);

unsigned char* generateRandomPoly(long k);

void generateRandomPoly(unsigned char *res, long k);

unsigned char* generateError(const unsigned char *codeword, long t, long n);

double* addError(const unsigned char *codeword, long t, long n);

void addNoise(double standartDeviation, const unsigned char *codeword, double *wordWithNoise, unsigned long n);

void printVec(const unsigned char *poly, int size);

void printVec(const unsigned long* poly, int size);

void printVec(const double*poly, int size);

void printVec(std::ofstream& out, const unsigned char *poly, int size);

void printVec(std::ofstream& out, const double*poly, int size);

void printMatrix(unsigned char** const matrix, int sizeI, int sizeJ = -1);

void printMatrix(std::ofstream& out, unsigned char** const matrix, int sizeI, int sizeJ = -1);

void makeMatrix(int power, const unsigned long * fieldElements, unsigned char** matrix);

#endif //BCH_BCHCODER_H
