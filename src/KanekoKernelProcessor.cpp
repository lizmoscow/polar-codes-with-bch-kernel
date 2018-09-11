//
// Created by Moskovskaya Elizaveta on 10/06/2018.
//

#include <cmath>
#include <algorithm>
#include <climits>
#include <cfloat>
#include <iostream>
#include <cstring>
#include "../headers/KanekoKernelProcessor.h"


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


double normalDistribution(double x, double sd, double expValue) {
    return 1 / (sqrt(2 * M_PI) * sd) * exp(-pow((x - expValue) / sd, 2) / 2);
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

//TODO: split into two functions to reduce the number of calls to memory
double KanekoKernelProcessor::calcT(long j) {
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

bool myFunction(std::pair<double, int> i, std::pair<double, int> j) { return (i.first < j.first); }

void KanekoKernelProcessor::set(const double *word) const {
    for (long i = 0; i < n; ++i) {
        alpha[i] = log(normalDistribution(word[i], sd, 1.0) / normalDistribution(word[i], sd, -1.0));
        alphaSorted[i].first = fabs(alpha[i]);
        alphaSorted[i].second = i;
        yH[i] = (alpha[i] <= 0.0) ? (unsigned char)0 : (unsigned char)1;
        alpha[i] = fabs(alpha[i]);
    }
    std::sort(alphaSorted, alphaSorted + n, myFunction);
}

void KanekoKernelProcessor::decode(unsigned char *res) {
    long j = 0;
    long i = 0;
    long T = LONG_MAX;
    double l = 0;
    double l0 = DBL_MAX;
    bool success;
    bool firstDecodingSuccessful = true;
    decoder.findSyndromPoly(yH);
    while (i < ((T == LONG_MAX) ? LONG_MAX : (1 << T))) {
        calcError(i);
        for (long k = 0; k < n; ++k) {
            err[k] ^= yH[k];
        }
        decoder.alterSyndromPoly(err);
        success = decoder.decode(err, x);
        ++decodingCount;
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
                ++comparisonCount;
                ++summCount;
            }
            T = j;
            j = 0;
            ++comparisonCount;
        }
        ++i;
        comparisonCount += n + 6;
        summCount += n + 1;
    }
    state = calcL();
}


void KanekoKernelProcessor::decode(const double *word, unsigned char *res) {
    for (int i = 0; i < n; ++i) {
        alpha[i] = log(normalDistribution(word[i], sd, 1.0) / normalDistribution(word[i], sd, -1.0));
        alphaSorted[i].first = fabs(alpha[i]);
        alphaSorted[i].second = i;
        yH[i] = (alpha[i] <= 0.0) ? (unsigned char)0 : (unsigned char)1;
        alpha[i] = fabs(alpha[i]);
    }
    std::sort(alphaSorted, alphaSorted + n, myFunction);
    summCount += 2 * n + 1;
    comparisonCount += 2 * n + 1;

    long j = 0;
    long i = 0;
    long T = LONG_MAX;
    double l = 0;
    double l0 = DBL_MAX;
    bool success;
    bool firstDecodingSuccessful = true;
    decoder.findSyndromPoly(yH);
    while (i < ((T == LONG_MAX) ? LONG_MAX : (1 << T))) {
        calcError(i);
        for (long k = 0; k < n; ++k) {
            err[k] ^= yH[k];
        }
        decoder.alterSyndromPoly(err);
        success = decoder.decode(err, x);
        ++decodingCount;
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
                ++comparisonCount;
                ++summCount;
            }
            T = j;
            j = 0;
            ++comparisonCount;
        }
        ++i;
        comparisonCount += n + 6;
        summCount += n + 1;
    }
    state = calcL();
}

long KanekoKernelProcessor::getN() const {
    return decoder.getN();
}

long KanekoKernelProcessor::getT() const {
    return decoder.getT();
}

long KanekoKernelProcessor::getK() const {
    return decoder.getK();
}

unsigned long KanekoKernelProcessor::getComparisonCount() const {
    return comparisonCount;
}

unsigned long KanekoKernelProcessor::getSummCount() const {
    return summCount;
}

unsigned long KanekoKernelProcessor::getDecodingCount() const {
    return decodingCount;
}

void KanekoKernelProcessor::setDecodingCount(unsigned long c) {
    decodingCount = c;
}

/*void* KanekoKernelProcessor::GetState(unsigned Stride) const {
    return state;
}*/
