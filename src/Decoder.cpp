//
// Created by Moskovskaya Elizaveta on 12.05.18.
//

#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstring>
#include "../headers/Decoder.h"

Decoder::Decoder(long pw, long n, long t, long k, unsigned long *antilogarithms, unsigned long *logarithms):
        power(pw), n(n), t(t), k(k) {

    l = 2 * t;

    this -> antilogarithms = antilogarithms;
    this -> logarithms = logarithms;
    lambda = new unsigned long[t + 1];
    locators = new unsigned long[t];
    syndromPoly = new unsigned long[l];
    oldPoly = new unsigned char[n];

    s = new unsigned long[l + 1];
    p = new unsigned long[l + 1];
    q = new unsigned long[l + 1];
    temp = new unsigned long[l + 1];
    a = {new unsigned long[l + 1], new unsigned long[l + 1],
                       new unsigned long[l + 1], new unsigned long[l + 1],
                       1, 1, 1 ,1};
    tempA = {new unsigned long[l + 1], new unsigned long[l + 1],
                           new unsigned long[l + 1], new unsigned long[l + 1],
                           1, 1, 1 ,1};

    f = new unsigned long[l + 1];
    shift = new unsigned long[l + 1];
    res1 = new unsigned long[l + 1];
}

Decoder::~Decoder() {
    delete[] lambda;
    delete[] locators;
    delete[] syndromPoly;
    delete[] oldPoly;

    delete[] s;
    delete[] p;
    delete[] q;
    delete[] temp;
    delete[] a.ff;
    delete[] a.fs;
    delete[] a.sf;
    delete[] a.ss;
    delete[] tempA.ff;
    delete[] tempA.fs;
    delete[] tempA.sf;
    delete[] tempA.ss;

    delete[] f;
    delete[] shift;
    delete[] res1;
}

unsigned long Decoder::multiply(unsigned long first, unsigned long second) const {
    if (!first || !second) return 0;
    auto temp = logarithms[first] + logarithms[second];
    if (temp < n) return antilogarithms[temp];
    return antilogarithms[temp - n];
}

unsigned long Decoder::divide(unsigned long first, unsigned long second) const {
    if (!first) return 0;
    auto temp = n + logarithms[first] - logarithms[second];
    if (temp < n) return antilogarithms[temp];
    return antilogarithms[temp - n];
}

void Decoder::addPolynomials(const unsigned long* first, int size1, const unsigned long* second, int size2,
                                            unsigned long *result, int *sizeRes) const{
    int m = std::min(size1, size2);
    for (unsigned long i = 0; i < m; ++i) {
        result[i] = first[i] ^ second[i];
        if (result[i]) {
            *sizeRes = i + 1;
        }
    }
    if (size1 > m) {
        for (unsigned long i = m; i < size1; ++i) {
            result[i] = first[i];
        }
        *sizeRes = size1;
    }
    else if (size2 > m){
        for (unsigned long i = m; i < size2; ++i) {
            result[i] = second[i];
        }
        *sizeRes = size2;
    }
}

void Decoder::multiplyPolynomials(const unsigned long* first, int size1, const unsigned long* second, int size2,
                                   unsigned long *result, int *size) const{
    for (int i = 0; i < size1 + size2 - 1; ++i) {
        result[i] = 0;
    }
    for (int i = 0; i < size1; ++i) {
        for (int j = 0; j < size2; ++j) {
            result[j + i] ^= multiply(first[i], second[j]);
        }
    }
    *size = size1 + size2 - 1;
}

void Decoder::dividePolynomial(unsigned long* first, int *size1, unsigned long* second, int *size2,
                               unsigned long *result, int *sizeResult) const{


    for (unsigned int i = 0; i < *size1; ++i) {
        f[i] = first[i];
    }
    int s = *size1, s2 = *size2;
    unsigned long coeff;
    for (unsigned int i = 0; i < s - s2 + 1; ++i) {
        shift[i] = 0;
        res1[i] = 0;
    }

    while (s >= s2) {
        shift[s - s2] = 1;
        coeff = divide(f[s - 1], second[s2 - 1]);

        for (int i = 0; i < *size2; ++i) {
            f[s - s2 + i] ^= multiply(coeff, second[i]);
        }

        res1[s - s2] = coeff;

        shift[s - s2] = 0;
        while (s > 0 && f[s - 1] == 0) {
            --s;
        }
    }
    *sizeResult = *size1 - *size2 + 1;
    for (int i = 0; i < *sizeResult; ++i) {
        result[i] = res1[i];
    }
    *size1 = *size2;
    for (int i = 0; i < *size2; ++i) {
        first[i] = second[i];
    }
    *size2 = (s) ? s : 1;
    for (int i = 0; i < *size2; ++i) {
        second[i] = f[i];
    }
}

unsigned long Decoder::eval(const unsigned long *poly, int size, unsigned long elem) const {
    unsigned long value = poly[size - 1];
    if (size < 2 ) {
        return value;
    }
    for (int i = size - 2; i >= 0; --i) {
        value = multiply(value, elem) ^ poly[i];
    }
    return value;
}



void Decoder::findSyndromPoly(const unsigned char *word)  {
    long start = n - 1;
    while (!word[start] && start > 0) {
        start--;
    }
    syndromPolySize = 0;

    for (unsigned long i = 1; i <= l; ++i) {
        syndromPoly[i - 1] = word[start];
        auto antilogI = antilogarithms[i];
        for (long j = start - 1; j >=0; --j) {
            syndromPoly[i - 1] = ((unsigned long)word[j]) ^ multiply(antilogI, syndromPoly[i - 1]);
        }
        if (syndromPoly[i - 1]) {
            syndromPolySize = i;
        }
    }

    memcpy(oldPoly, word, n * sizeof(*word));
}

//TODO: simplify
void Decoder::alterSyndromPoly(const unsigned char *word) {
    for (int i = 0; i < n; ++i){
        if (oldPoly[i] != word[i]) {
            for (int j = 1; j <= l; ++j) {
                syndromPoly[j - 1] ^= multiply(oldPoly[i], antilogarithms[(j * i) % n]);
                syndromPoly[j - 1] ^= multiply(word[i], antilogarithms[(j * i) % n]);
            }
        }
    }

    syndromPolySize = l;
    for (int i = l - 1; i >= 0 && !syndromPoly[i]; --i) {
        syndromPolySize = i;
    }

    memcpy(oldPoly, word, n * sizeof(*word));
}


bool Decoder::euclid(const unsigned char *word, unsigned long *lambda, int *size) {
    sizeS = l + 1;
    for (int i = 0; i < l; ++i) {
        s[i] = 0;
    }
    s[l] = 1;


    memcpy(p, syndromPoly, syndromPolySize * sizeof(*syndromPoly));
    sizeP = syndromPolySize;

    for (int i = 0; i <= l; ++i) {
        a.ff[i] = 0; a.fs[i] = 0;
        a.sf[i] = 0; a.ss[i] = 0;
    }
    a.ff[0] = 1;
    a.ss[0] = 1;
    a.sizeFF = 1; a.sizeFS = 1;
    a.sizeSF = 1; a.sizeSS = 1;

    while (sizeP > t) {
        dividePolynomial(s, &sizeS, p, &sizeP, q, &sizeQ);

        addPolynomials(a.sf, a.sizeSF, nullptr, 0, tempA.ff, &tempA.sizeFF);

        addPolynomials(a.ss, a.sizeSS, nullptr, 0, tempA.fs, &tempA.sizeFS);

        multiplyPolynomials(q, sizeQ, a.sf, a.sizeSF, temp, &sizeTemp);
        addPolynomials(a.ff, a.sizeFF, temp, sizeTemp, tempA.sf, &tempA.sizeSF);

        multiplyPolynomials(q, sizeQ, a.ss, a.sizeSS, temp, &sizeTemp);
        addPolynomials(a.fs, a.sizeFS, temp, sizeTemp, tempA.ss, &tempA.sizeSS);

        std::swap(tempA, a);
    }


    unsigned long delta = a.ss[0];
    if (!delta) {
        return false;
    }
    *size = a.sizeSS;
    memcpy(lambda, a.ss, a.sizeSS * sizeof(*a.ss));
    return true;
}

bool Decoder::locatorsAndRoots(const unsigned long *lambda, int size, unsigned long *locators) const {
    int count = 0;
    unsigned long k = 0;
    for (unsigned long i = 0; i < size; ++i) {
        while (k < n && eval(lambda, size, antilogarithms[k]) != 0) {
            ++k;
        }
        if (k < n) {
            locators[i] = (n - k) % n;
            count = i;
            ++k;
        }
    }
    if (count != size - 2) {
        return false;
    }
    return true;
}

bool Decoder::decode(const unsigned char* word, unsigned char *answer){
    auto res = euclid(word, lambda, &sizeLambda);
    if (!res) {
        return false;
    }

    res = locatorsAndRoots(lambda, sizeLambda, locators);
    if (!res) {
        return false;
    }

    for (int i = 0; i < n; ++i) {
        answer[i] = 0;
    }

    for (int i = 0; i < sizeLambda - 1; ++i) {
        answer[ locators[i]] = 1;
    }

    for (int i = 0; i < n; ++i) {
        answer[i] ^= word[i];
    }
    return true;
}

long Decoder::getN() const {
    return n;
}

long Decoder::getT() const {
    return t;
}

long Decoder::getK() const {
    return k;
}

