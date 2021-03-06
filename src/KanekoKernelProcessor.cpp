//
// Created by Moskovskaya Elizaveta on 10/06/2018.
//

#include <cmath>
#include <algorithm>
#include <climits>
#include <cfloat>
#include <iostream>
#include <cstring>
#include <cassert>
#include "../headers/KanekoKernelProcessor.h"
#include "../headers/bchCoder.h"

#define COUNT

KanekoKernelProcessor::KanekoKernelProcessor(long pw, long n, long t, long k, unsigned long *antilogarithms, unsigned long *logarithms, double signalToNoiseRatio):
        decoder(pw, n, t, k, antilogarithms, logarithms), n(n), t(t) {

    sd = sqrt(1 / (pow(10, signalToNoiseRatio / 10) * 2 * k / n));
    alpha = new double[n];
    alphaSorted = new std::pair<double, int>[n];
    yH = new unsigned char[n];
    x = new unsigned char[n];
    err = new unsigned char[n];
}

KanekoKernelProcessor::~KanekoKernelProcessor() {
    delete[] alpha;
    delete[] alphaSorted;
    delete[] yH;
    delete[] x;
    delete[] err;
}

void KanekoKernelProcessor::calcError(long i) {
    long pos = 0;

    for (long k = 0; k < n; ++k) {
        err[k] = 0;
    }

    while (i > 0) {
        if (i & 1) {
            err[alphaSorted[pos].second] = 1;
        }

        ++pos;
        i >>= 1;
    }
}

//TODO: save alphaSoterd[j].second into a variable
double KanekoKernelProcessor::calcRightSide() const {
    long border = (2 * t + 1) - (m + m0) / 2;
    double l = 0;
    long i = 0;
    long j = 0;
    while (i < border && j < n) {
        if (yH[alphaSorted[j].second] == x[alphaSorted[j].second]) {
            l += alphaSorted[j].first;
            ++i;
        }
        ++j;
    }
    return l;
}

double KanekoKernelProcessor::calcL() const {
    double l = 0;
    for (long i = 0; i < n; ++i) {
        if (yH[i] != x[i]) {
            l += alpha[i];
        }
    }
    return l;
}

double KanekoKernelProcessor::calcL(const unsigned char *word) const {
    double l = 0;
    for (long i = 0; i < n; ++i) {
        if (yH[i] != word[i]) {
            l += alpha[i];
        }
    }
    return l;
}

long KanekoKernelProcessor::calcM() {
    long count = 0;
    for (int i = 0; i < n; ++i) {
        if (yH[i] != x[i]) {
            ++count;
        }
    }
    return count;
}

long KanekoKernelProcessor::calcM(const unsigned long *answer) {
    long count = 0;
    for (int i = 0; i < n; ++i) {
        if (yH[i] != answer[i]) {
            ++count;
        }
    }
    return count;
}


double KanekoKernelProcessor::calcT(uint64_t j) {
    long border = t - (m + m0) / 2;
    long i = 0;
    long k = 0;
    double l = 0;
    while (i < border) {
        if (yH[alphaSorted[k].second] == x[alphaSorted[k].second]) {
            l += alphaSorted[k].first;
            ++i;
        }
        ++k;
    }
    for (i = 0; i <= t; ++i) {
        l += alphaSorted[j + i].first;
    }
    return l;
}

double KanekoKernelProcessor::calcT(const unsigned long *answer, uint64_t j) {
    long mLocal = calcM(answer);
    long border = t - (mLocal + m0) / 2;
    long i = 0;
    long k = 0;
    double l = 0;
    while (i < border) {
        if (yH[alphaSorted[k].second] == answer[alphaSorted[k].second]) {
            l += alphaSorted[k].first;
            ++i;
        }
        ++k;
    }
    for (i = 0; i <= t; ++i) {
        l += alphaSorted[j + i].first;
    }
    return l;
}


bool myFunction(std::pair<double, int> i, std::pair<double, int> j) { return (i.first < j.first); }

void KanekoKernelProcessor::set(const double *word) const {
    for (long i = 0; i < n; ++i) {
        alpha[i] = 2 * word[i] / pow(sd, 2);;
        alphaSorted[i].first = fabs(alpha[i]);
        alphaSorted[i].second = i;
        yH[i] = (alpha[i] <= 0.0) ? (unsigned char)0 : (unsigned char)1;
        alpha[i] = fabs(alpha[i]);
    }
    std::sort(alphaSorted, alphaSorted + n, myFunction);
}

void KanekoKernelProcessor::decode(unsigned char *res) {
    uint64_t j = 0;
    long i = 0;
    //TODO: why T=LONG_MAX?
    uint64_t T = LONG_MAX;
    double l = 0;
    double l0 = DBL_MAX;
    bool success;
    bool firstDecodingSuccessful = true;
    decoder.findSyndromPoly(yH);
    while (i < ((T == LONG_MAX) ? LONG_MAX : (uint64_t(1) << T))) {
        calcError(i);
        for (long k = 0; k < n; ++k) {
            err[k] ^= yH[k];
        }
        decoder.alterSyndromPoly(err);
        if (!i && !success) firstDecodingSuccessful = false;
        if (success) {
            m = calcM();
            if (!i || !firstDecodingSuccessful) m0 = m;
            l = calcL();
            //if (l == answerL) {static int count = 0; ++count; std::cout << count << " ";}
            if (l < l0) {
                memcpy(res, x, n * sizeof(*res));
                l0 = l;
                if (l < calcRightSide()) {
                    return;
                }
                //TODO: j <=  n - 1 - t
                while (l >= calcT(j) && j <  n - 1 - t) {
                    ++j;
#ifdef COUNT
                    ++comparisonCount;
                    ++summCount;
#endif
                }
                T = (j > J) ? J : j;
                j = 0;
#ifdef COUNT
                ++comparisonCount;
#endif
            }
        }
        ++i;
#ifdef COUNT
        comparisonCount += n + 6;
        summCount += n + 1;
#endif
    }
}

void KanekoKernelProcessor::decode(const double *word, unsigned char *res) {
    for (int i = 0; i < n; ++i) {
        alpha[i] = 2 * word[i] / pow(sd, 2);
        alphaSorted[i].first = fabs(alpha[i]);
        alphaSorted[i].second = i;
        yH[i] = (alpha[i] <= 0.0) ? (unsigned char)0 : (unsigned char)1;
        alpha[i] = fabs(alpha[i]);
    }
    std::sort(alphaSorted, alphaSorted + n, myFunction);
#ifdef COUNT
    summCount += 2 * n + 1;
    comparisonCount += 2 * n + 1;
#endif

    uint64_t j = 0;
    long i = 0;
    //TODO: why T=LONG_MAX?
    uint64_t T = LONG_MAX;
    double l = 0;
    double l0 = DBL_MAX;
    bool success;
    bool firstDecodingSuccessful = true;
    decoder.findSyndromPoly(yH);
    //TODO: i < (uint64_t(1) << T) - 1
    while (i < ((T == LONG_MAX) ? LONG_MAX : (uint64_t(1) << T))) {
        calcError(i);
        for (long k = 0; k < n; ++k) {
            err[k] ^= yH[k];
        }
        decoder.alterSyndromPoly(err);
        success = decoder.decode(err, x);
#ifdef COUNT
        ++decodingCount;
#endif
        m = calcM();
        if (!i && success) m0 = m;
        else if (!i) firstDecodingSuccessful = false;
        if (!firstDecodingSuccessful) m0 = m;
        l = calcL();
        if (success && l < l0) {
            memcpy(res, x, n * sizeof(*res));
            l0 = l;
            if (l < calcRightSide()) {
                return;
            }
            while (l >= calcT(j)) {
                ++j;
#ifdef COUNT
                ++comparisonCount;
                ++summCount;
#endif
            }
            T = j;
            j = 0;
#ifdef COUNT
            ++comparisonCount;
#endif
        }
        ++i;
#ifdef COUNT
        comparisonCount += n + 6;
        summCount += n + 1;
#endif
    }
}

/*void KanekoKernelProcessor::decode(const unsigned char *answer, const double *word, unsigned char *res) {
    for (int i = 0; i < n; ++i) {
        alpha[i] = 2 * word[i] / pow(sd, 2);
        alphaSorted[i].first = fabs(alpha[i]);
        alphaSorted[i].second = i;
        yH[i] = (alpha[i] <= 0.0) ? (unsigned char)0 : (unsigned char)1;
        alpha[i] = fabs(alpha[i]);
    }
    std::sort(alphaSorted, alphaSorted + n, myFunction);

    double answerL = calcL(answer);
    uint64_t z = 0;
    while (answerL >= calcT(z)) {
        ++z;
    }

    uint64_t j = 0;
    long i = 0;
    uint64_t T = LONG_MAX;
    double l = 0;
    double l0 = DBL_MAX;
    bool success;
    bool firstDecodingSuccessful = true;
    decoder.findSyndromPoly(yH);
    while (i < ((T == LONG_MAX) ? LONG_MAX : (uint64_t(1) << T))) {
        calcError(i);
        for (long k = 0; k < n; ++k) {
            err[k] ^= yH[k];
        }
        decoder.alterSyndromPoly(err);
        success = decoder.decode(err, x);
        if (!i && !success) firstDecodingSuccessful = false;
        if (success) {
            m = calcM();
            if (!i || !firstDecodingSuccessful) m0 = m;
            l = calcL();
            //if (l == answerL) {static int count = 0; ++count; std::cout << count << " ";}
            if (l < l0) {
                memcpy(res, x, n * sizeof(*res));
                l0 = l;
                if (l < calcRightSide()) {
                    return;
                }
                while (l >= calcT(j)) {
                    ++j;
                    assert(j<n);
                }
                T = j;
                j = 0;
            }
        }
        ++i;
    }
    std::cout << i << "\n";
}*/

//TODO: delete answer
void KanekoKernelProcessor::decode(const unsigned char *answer, const double *word, unsigned char *res) {
    for (int i = 0; i < n; ++i) {
        alpha[i] = 2 * word[i] / pow(sd, 2);
        alphaSorted[i].first = fabs(alpha[i]);
        alphaSorted[i].second = i;
        yH[i] = (alpha[i] <= 0.0) ? (unsigned char)0 : (unsigned char)1;
        alpha[i] = fabs(alpha[i]);
    }
    std::sort(alphaSorted, alphaSorted + n, myFunction);

    double answerL = calcL(answer);
    /*uint64_t z = 0;
    while (answerL >= calcT(z)) {
        ++z;
    }*/

    long j = 0;
    long i = 0;
    //TODO: why T=10?
    long T = n;
    double l = 0;
    double l0 = DBL_MAX;
    bool success;
    bool firstDecodingSuccessful = true;
    decoder.findSyndromPoly(yH);
    //TODO: i < (1 << T) - 1
    while (i < (1 << T) - 1) {
        calcError(i);
        for (long k = 0; k < n; ++k) {
            err[k] ^= yH[k];
        }
        decoder.alterSyndromPoly(err);
#ifdef COUNT
        ++decodingCount;
#endif
        success = decoder.decode(err, x);
        if (!i && !success) firstDecodingSuccessful = false;
        if (success) {
            m = calcM();
            if (!i || !firstDecodingSuccessful) m0 = m;
            l = calcL();
            //if (l == answerL) {static int count = 0; ++count; std::cout << count << " ";}
            if (l < l0) {
                memcpy(res, x, n * sizeof(*res));
                l0 = l;
                if (l < calcRightSide()) {
                    return;
                }
                //TODO: j <=  n - 1 - t
                while (l >= calcT(j) && j <=  n - 1 - t) {
                    ++j;
#ifdef COUNT
                    ++comparisonCount;
                    ++summCount;
#endif
                }
                //TODO: why?
                //T = (j > J) ? J : j;
                T = j;
                j = 0;
#ifdef COUNT
                ++comparisonCount;
#endif
            }
        }
        ++i;
#ifdef COUNT
        comparisonCount += n + 6;
        summCount += n + 1;
#endif
    }
    //std::cout << i << "\n";
}

long KanekoKernelProcessor::getN() const { return decoder.getN(); }

long KanekoKernelProcessor::getT() const { return decoder.getT(); }

long KanekoKernelProcessor::getK() const { return decoder.getK(); }

unsigned long KanekoKernelProcessor::getComparisonCount() const { return comparisonCount; }

unsigned long KanekoKernelProcessor::getSummCount() const { return summCount; }

unsigned long KanekoKernelProcessor::getDecodingCount() const { return decodingCount; }

void KanekoKernelProcessor::setDecodingCount(unsigned long c) { decodingCount = c; }

void KanekoKernelProcessor::setComparisonCount(unsigned long comparisonCount) { this -> comparisonCount = comparisonCount; }

void KanekoKernelProcessor::setSummCount(unsigned long summCount) { this -> summCount = summCount; }

