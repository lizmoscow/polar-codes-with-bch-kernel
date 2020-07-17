//
// Created by Moskovskaya Elizaveta on 10.05.18.
//

#include <random>
#include <iostream>
#include <fstream>
#include <cstring>
#include <chrono>
#include <stdlib.h>
#include <LinAlg.h>
#include <TrellisKernelProcessor.h>
#include <SectionedTrellisKernelProcessor.h>
#include "bchCoder.h"

#define RANDOM
//#define DEBUG

#ifdef RANDOM
auto seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
#endif

#ifndef RANDOM
std::default_random_engine generator(6);
#endif
std::uniform_int_distribution<unsigned short> distribution(0, 1);


void findMinimalPolynomial(int i, int power, const unsigned long * fieldElements, int *size, unsigned char * res) {
    //basic definition and initialization
    auto minimalPolynomial = new unsigned long[power + 1];
    int sizeMP = 1;
    minimalPolynomial[0] = 1;
    for (int ind = 1; ind <= power; ++ind) {
        minimalPolynomial[ind] = 0;
    }
    auto binom = new unsigned long[2];
    binom[1] = 1;
    auto temp = new unsigned long[power + 1];
    for (int ind = 0; ind <= power; ++ind) {
        temp[ind] = 0;
    }
    unsigned long element = fieldElements[i];
    int index = 2;
    int i1, i2;

    //calculating minimal polynomial
    do {
        binom[0] = element;
        for (int ind1 = 0; ind1 < sizeMP; ++ind1) {
            for (int ind2 = 0; ind2 < 2; ++ind2) {
                if (minimalPolynomial[ind1]) {
                    //finding orders of elements which will be multiplied
                    for (int ind = 0; ind < (1 << power) - 1; ++ind) {
                        if (minimalPolynomial[ind1] == fieldElements[ind]) {
                            i1 = ind;
                            break;
                        }
                    }
                    for (int ind = 0; ind < (1 << power) - 1; ++ind) {
                        if (binom[ind2] == fieldElements[ind]) {
                            i2 = ind;
                            break;
                        }
                    }
                    //adding the result of multiplication to the coefficient of polynomial
                    temp[ind1 + ind2] ^= fieldElements[(i1 + i2) % ((1 << power) - 1)];
                }

            }
        }
        //cleaning up
        std::swap(temp, minimalPolynomial);
        ++sizeMP;
        for (unsigned long ind = 0; ind <= power; ++ind) {
            temp[ind] = 0;
        }
        //getting next root
        element = fieldElements[index * i % ((1 << power) - 1)];
        index *= 2;
    } while (element != fieldElements[i]);

    //forming the result
    *size = sizeMP;
    //auto res = new unsigned char[*size];
    for (int ind = sizeMP - 1; ind >= 0; --ind) {
        res[ind] = ((minimalPolynomial[ind] % 2) ? 1 : 0);
    }

    delete[] minimalPolynomial;
    delete[] binom;
    delete[] temp;

    //return res;
}
bool comparePoly(const unsigned char* poly1, int size1, const unsigned char* poly2, int size2) {
    if (size1 != size2) {
        return false;
    }
    for (int i = 0; i < size1; ++i) {
        if (poly1[i] != poly2[i]) {
            return false;
        }
    }
    return true;
}

unsigned char* multiplyPolynomials(const unsigned char* first, int size1, const unsigned char* second, int size2, int *sizeRes) {
    auto res = new unsigned char[size1 + size2 - 1];
    for (int i = 0; i < size1 + size2 - 1; ++i) {
        res[i] = 0;
    }
    for (int i = 0; i < size1; ++i) {
        for (int j = 0; j < size2; ++j) {
            res[j + i] ^= first[i] & second[j];
        }
    }
    if (sizeRes != nullptr) {
        *sizeRes = size1 + size2 - 1;
    }
    return res;
}

void multiplyPolynomials(const unsigned char* first, int size1,
                        const unsigned char* second, int size2,
                        unsigned char *res, int *sizeRes, int offset) {
    for (int i = 0; i < size1 + size2 - 1; ++i) {
        res[i + offset] = 0;
    }
    for (int i = 0; i < size1; ++i) {
        for (int j = 0; j < size2; ++j) {
            res[j + i + offset] ^= first[i] & second[j];
        }
    }
    if (sizeRes != nullptr) {
        *sizeRes = size1 + size2 - 1;
    }
}

unsigned char* dividePolynomial(const unsigned char* first, int size1, const unsigned char* second, int size2, int *size, bool needRemainder) {
    auto *a = new unsigned char[size1];
    for (unsigned int i = 0; i < size1; ++i) {
        a[i] = first[i];
    }
    int s = size1, s2 = size2;
    auto *shift = new unsigned char[s - s2 + 1];
    auto *res1 = new unsigned char[s - s2 + 1];
    for (unsigned int i = 0; i < s - s2 + 1; ++i) {
        shift[i] = 0;
        res1[i] = 0;
    }

    while (s >= s2) {
        shift[s - s2] = 1;

        for (int i = 0; i < size2; ++i) {
            a[s - s2 + i] = a[s - s2 + i] ^ second[i];
        }

        res1[s - s2] = 1;

        shift[s - s2] = 0;
        while (s > 0 && a[s - 1] == 0) {
            --s;
        }
    }

    if (needRemainder) {
        *size = (s) ? s : 1;
        auto remainder = new unsigned char[*size];
        for (int i = 0; i < *size; ++i) {
            remainder[i] = (a[i]) ? 1 : 0;
        }
        delete[] a;
        delete[] res1;
        delete[] shift;
        return remainder;
    }
    else {
        *size = size1 - size2 + 1;
        auto quotient = new unsigned char[*size];
        for (int i = 0; i < size1 - size2 + 1; ++i) {
            quotient[i] = (res1[i]) ? 1 : 0;
        }
        delete[] a;
        delete[] res1;
        delete[] shift;
        return quotient;
    }
}

unsigned char* gcd(const unsigned char* first, int size1, const unsigned char* second, int size2, int *sizeRes) {
    auto *zero = new unsigned char[1];
    zero[0] = 0;
    unsigned char *a, *c;
    int sa, sc, st;

    auto res = new unsigned char[size1];
    *sizeRes = size1;
    for (unsigned int i = 0; i < size1; ++i) {
        res[i] = first[i];
    }

    c = new unsigned char[size2];
    sc = size2;
    for (unsigned int i = 0; i < size2; ++i) {
        c[i] = second[i];
    }

    while (!comparePoly(c, sc, zero, 1)) {
        a = res;
        sa = *sizeRes;
        res = c;
        *sizeRes = sc;
        c = dividePolynomial(a, sa, res, *sizeRes, &sc, true);
        delete[] a;
    }
    delete[] zero;
    delete[] c;
    return res;
}

unsigned char* lcm(const unsigned char* first, int size1, const unsigned char* second, int size2, int *sizeRes) {
    unsigned char *temp1, *temp2, *res;
    int st1, st2;
    temp1 = gcd(first, size1, second, size2, &st1);
    temp2 = multiplyPolynomials(first, size1, second, size2, &st2);
    res = dividePolynomial(temp2, st2, temp1, st1, sizeRes, false);
    delete[] temp1;
    delete[] temp2;
    return res;
}

unsigned char* generateRandomPoly(long k) {
    auto res = new unsigned char[k];
    for (int i = 0; i < k; ++i) {
        res[i] = (unsigned char)distribution(generator);
    }
    return res;
}

void generateRandomPoly(unsigned char *res, long k) {
    for (int i = 0; i < k; ++i) {
        res[i] = (unsigned char)distribution(generator);
    }
}


void addNoise(double standartDeviation, const unsigned char *codeword, double *wordWithNoise, unsigned long n) {
    std::normal_distribution<double> distribution(0.0, standartDeviation);
    double g;
    for (int i = 0; i < n; ++i) {
        g = distribution(generator);
        wordWithNoise[i] = (codeword[i] ? 1 : -1) + g;
    }
}


void printVec(const double *poly, int size) {
    for (int i = 0; i < size; ++i) {
        std::cout << (poly[i]) << ' ';
    }
    std::cout << std::endl;
}


void printVec(const unsigned char* poly, int size) {
    for (int i = 0; i < size; ++i) {
        std::cout << ((poly[i]) ? 1 : 0) << ' ';
    }
    std::cout << std::endl;
}

void printVec(const unsigned long* poly, int size) {
    for (int i = 0; i < size; ++i) {
        std::cout << poly[i] << ' ';
    }
    std::cout << std::endl;
}

void printVec(const int* poly, int size) {
    for (int i = 0; i < size; ++i) {
        std::cout << poly[i] << ' ';
    }
    std::cout << std::endl;
}


void printVec(std::ofstream& out, const double *poly, int size) {
    for (int i = 0; i < size; ++i) {
        out << (poly[i]) << ' ';
    }
    out << std::endl;
}


void printVec(std::ofstream& out, const unsigned char* poly, int size) {
    for (int i = 0; i < size; ++i) {
        out << ((poly[i]) ? 1 : 0) << ' ';
    }
    out << std::endl;
}


/*void printMatrix(unsigned char** const matrix, int sizeI, int sizeJ) {
    if (sizeJ == -1) sizeJ = sizeI;
    for (auto i = 0; i < sizeI; ++i) {
        for (auto j = 0; j < sizeJ; ++j) {
            std::cout << ((matrix[i][j]) ? 1 : 0) << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}*/


void printMatrix(const unsigned char* const matrix, int sizeI, int sizeJ) {
    if (sizeJ == -1) sizeJ = sizeI;
    for (auto i = 0; i < sizeI; ++i) {
        for (auto j = 0; j < sizeJ; ++j) {
            std::cout << ((matrix[i * sizeI + j]) ? 1 : 0) << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


void printMatrix(std::ofstream& out, unsigned char** const matrix, int sizeI, int sizeJ) {
    if (sizeJ == -1) sizeJ = sizeI;
    for (auto i = 0; i < sizeI; ++i) {
        for (auto j = 0; j < sizeJ; ++j) {
            out << ((matrix[i][j]) ? 1 : 0) << ' ';
        }
        out << std::endl;
    }
    out << std::endl;

}


void printMatrix(std::ofstream& out, const unsigned char* const matrix, int sizeI, int sizeJ) {
    if (sizeJ == -1) sizeJ = sizeI;
    if (sizeJ == -1) sizeJ = sizeI;
    for (auto i = 0; i < sizeI; ++i) {
        for (auto j = 0; j < sizeJ; ++j) {
            out << ((matrix[i * sizeI + j]) ? 1 : 0) << ' ';
        }
        out << std::endl;
    }
    out << std::endl;
}


void makeMatrix(int power, const unsigned long * fieldElements, unsigned char* matrix) {
    const auto amount = (power != 2) ? ((1 << power) - 2) / 2 : 2;
    const auto len = (1 << power) - 1;
    int polySize, gSizeOld = 1, gSizeNew, size;
    auto poly = new unsigned char[power + 1];
    auto g = new unsigned char[len];
    g[0] = 1;
    for (auto i = 0; i < len + 1; ++i) {
        matrix[i * (len + 1)] = 1;
        for (auto j = 1; j < len + 1; ++j) {
            matrix[i * (len + 1) + j] = 0;
        }
    }
    matrix[len + 2] = 1;
    for (unsigned long i = 2; i <= amount; i++) {
        findMinimalPolynomial(i, power, fieldElements, &polySize, poly);
        if (gSizeOld >= polySize && !dividePolynomial(g, gSizeOld, poly, polySize, &size, true)[0] && size == 1) { continue; }
        gSizeNew = polySize + gSizeOld - 1;
        multiplyPolynomials(poly, polySize, g, gSizeOld, matrix + gSizeNew * (len + 1), nullptr, 1);
        auto count = 1;
        for (auto j = gSizeOld; j < gSizeNew - 1; ++j) {
            for (auto k = count; k < count + gSizeOld; ++k) {
                matrix[(j + 1) * (len + 1) + k + 1] = g[k - count];
            }
            ++count;
        }
        for (auto j = 1; j <= gSizeNew; ++j) {
            g[j - 1] = matrix[gSizeNew * (len + 1) + j];
        }
        gSizeOld = gSizeNew;
    }
    delete[] g;
    delete[] poly;
}


/*void swapColumns(int power, const unsigned long * fieldElements, unsigned char** matrix, unsigned char** newMatrix) {
    int n = 1 << power;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= 2; ++j) {
            newMatrix[i][j] = matrix[i][j];
        }
    }
    for (int i = 3; i < n; ++i) {
        int j = 2;
        for (; j < n; ++j) {
            if (fieldElements[j] == i) {
                break;
            }
        }
        for (int k = 0; k < n; ++k) {
            newMatrix[k][i] = matrix[k][j + 1];
        }
    }
}*/


/*void swapColumns(int power, const unsigned long* fieldElements, unsigned char* matrix, CMatrixBinaryKernel* newKernel, unsigned int length) {
    int n = 1 << power;
    auto oldBasis = new unsigned int[power];
    auto newBasis = new unsigned int[power];
    auto b = new unsigned char[power * power];
    unsigned char* tempMatrix = new unsigned char[n * n];
    CMatrixBinaryKernel* tempKernel;
    unsigned long long int minSum = 0, minCmp = 0, sum = 0, cmp = 0;
    for (int i = 0; i < power; ++i) {
        oldBasis[i] = 1 << i;
    }
    newBasis[0] = 0;
    for (int i = 0; i < power * power; ++i) {
        b[i] = 0;
    }
    unsigned long border = (unsigned long)1 << (power * power);
    for (int i = 0; i < border; ++i) {
        if (increase(power, b)) {
            for (int j = 0; j < power; ++j) {
                newBasis[j] = 0;
                for (int k = 0; k < power; ++k) {
                    newBasis[j] ^= oldBasis[k] * b[power * k + j];
                }
            }
            for (int j = 1; j < n; ++j) {
                int temp = 0;
                for (int k = 0; k < power; ++k) {
                    temp += newBasis[k] * (j & (1 << k));
                }
                for (int k = 0; k < n; ++k) {
                    tempMatrix[k * n + temp] = matrix[k * n + j];
                }
            }
            tempKernel = new CMatrixBinaryKernel(length, tempMatrix, "name", nullptr);
            CKernProcLLR* processor = tempKernel->GetProcessor(0);
            auto processorState = processor->GetState(1);
            auto processorTemp = processor->GetTemp(1);
            sum = SumCount;
            cmp = CmpCount;
            tBit* pKnownInputSymbols = new tBit[processor->Size()];
            MType* pChannelLLRs = new MType[processor->Size()];
            MType* pLLRs = new MType[1];
            for (int j = 0; j < processor->Size(); ++j) {
                pKnownInputSymbols[j] = 0;
            }
            for (int j = 0; j < n; ++j) {
                processor->GetLLRs(1, j, pKnownInputSymbols, pChannelLLRs, pLLRs, processorState, processorTemp);
            }
            if (SumCount - sum < minSum && CmpCount - cmp < minCmp) {
                minSum = SumCount - sum;
                minCmp = CmpCount - cmp;
                std::memcpy(newKernel, tempKernel, sizeof(tempKernel));
            }
            processor->FreeState(processorState);
            processor->FreeTemp(processorTemp);
            delete tempKernel;
            delete processor;
            delete[] pKnownInputSymbols;
            delete[] pChannelLLRs;
            delete[] pLLRs;
        }
    }
}*/


void swapColumns(int power, const unsigned long * fieldElements, unsigned char* matrix, tBit* newMatrix, unsigned int length, unsigned int sectionSize, unsigned int codeLength) {
    unsigned long long int minSum = ~(unsigned long long int)0, minCmp = ~(unsigned long long int)0, sum = 0, cmp = 0, prevSum = SumCount, prevCmp = CmpCount;
    int n = 1 << power;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= 2; ++j) {
            newMatrix[i * n + j] = matrix[i * n + j];
        }
    }
    for (int i = 3; i < n; ++i) {
        int j = 2;
        for (; j < n; ++j) {
            if (fieldElements[j] == i) {
                break;
            }
        }
        for (int k = 0; k < n; ++k) {
            newMatrix[k * n + i] = matrix[k * n + j + 1];
        }
    }
    auto tempKernel = new CMatrixBinaryKernel(length, newMatrix, "name", nullptr);
    auto processor = (SectionedTrellisKernelProcessor*)(tempKernel->GetProcessor(0, sectionSize, codeLength));
    auto processorState = processor->GetState(1);
    auto processorTemp = processor->GetTemp(1);
    auto processorMetrics = processor->GetMetrics();
    sum = SumCount;
    cmp = CmpCount;
    tBit* pKnownInputSymbols = new tBit[processor->Size()];
    MType* pChannelLLRs = new MType[processor->Size()];
    MType* pLLRs = new MType[1];
    for (int j = 0; j < processor->Size(); ++j) {
        pKnownInputSymbols[j] = 0;
    }
    for (int j = 0; j < processor -> Size(); ++j) {
        pChannelLLRs[j] = (rand() % 10) - 5;
    }
    for (int j = 0; j < n; ++j) {
        processor->GetLLRs(1, j, pKnownInputSymbols, pChannelLLRs, pLLRs, processorState, processorTemp, processorMetrics);
    }
    processor->FreeState(processorState);
    processor->FreeTemp(processorTemp);
    processor->FreeMetrics(processorMetrics);
    delete processor;
    delete tempKernel;
    delete[] pKnownInputSymbols;
    delete[] pChannelLLRs;
    delete[] pLLRs;
#ifdef DEBUG
    printMatrix(newMatrix, n, n);
#endif
    std::cout << (SumCount - sum) << " " << (CmpCount - cmp) << "\n";
}


/*unsigned int logarithm(unsigned int base, unsigned int val) {
    unsigned int i = 0;
    while (val / base > 0) {
        val /= base;
        i++;
    }
    return i;
}*/


void randomSwapColumns(int power, unsigned char* matrix, tBit* newKernel, unsigned int length, unsigned int sectionSize, unsigned int codeLength) {
    int n = 1 << power;
    auto oldBasis = new unsigned long[power];
    auto newBasis = new unsigned long[power];
    auto b = new unsigned char[power * power];
    auto l = new unsigned char[power * power];
    auto u = new unsigned char[power * power];
    auto tempMatrix = new tBit[n * n];
    auto tempArray = new int[n];
    auto finalTempArray = new int[n];
    auto finalBasis = new unsigned long[power];
    //CMatrixBinaryKernel tempKernel;
    unsigned long long int minSum = ~(unsigned long long int)0, minCmp = ~(unsigned long long int)0, sum = 0, cmp = 0, prevSum = SumCount, prevCmp = CmpCount;
#ifdef DEBUG
    //printVec(oldBasis, power);
#endif
    MType* pChannelLLRs = new MType[n];
    tBit* pKnownInputSymbols = new tBit[n];
    MType* pLLRs = new MType[1];
    for (int j = 0; j < n; ++j) {
        pKnownInputSymbols[j] = 0;
    }
    for (int j = 0; j < n; ++j) {
        pChannelLLRs[j] = (rand() % 10) - 5;
    }
    auto processorState = AlignedAlloc(n);
    auto processorTemp = malloc(sizeof(MType) * 2  << n);
    unsigned int size1 = logarithm(n, codeLength);
    unsigned int size2;
    unsigned int N = n + 1;
    unsigned int L = (N % sectionSize) ? (N / sectionSize + (N % sectionSize)) : (N / sectionSize - 1  + (N % sectionSize) + sectionSize);
    KType**** processorMetrics = new KType***[size1];
    for (int i = 0; i < size1; ++i) {
        size2 = (unsigned int)pow(n, i);
        processorMetrics[i] = new KType**[size2];
        for (int j = 0; j < size2; ++j) {
            processorMetrics[i][j] = new KType*[L - 1];
            for (int k = 0; k < L - 1; ++k) {
                processorMetrics[i][j][k] = new KType[1 << sectionSize];
            }
        }
    }

    //filling in upper- and lower-triangular matrices with zeroes and ones in diagonals
    for (int i = 0; i < power; ++i) {
        l[i * power + i] = u[i * power + i] = 1;
        for (int j = i + 1; j < power; ++j) {
            l[i * power + j] = 0;
        }
        for (int j = 0; j < i; ++j) {
            u[i * power + j] = 0;
        }
    }
    unsigned long border = 20000000;

    //trying different permutations
    for (int i = 0; i < border; ++i) {
        //creating a standard basis
        for (int i = 0; i < power; ++i) {
            oldBasis[i] = (unsigned long)1 << i;
        }
        //generating random permutation
        randomInvertibleMatrix(power, l, u, b);
#ifdef DEBUG1
        printMatrix(b, power, power);
#endif
        //applying it to old basis
        for (int j = 0; j < power; ++j) {
            newBasis[j] = 0;
            for (int k = 0; k < power; ++k) {
                newBasis[j] ^= oldBasis[k] * b[power * k + j];
            }
        }
#ifdef DEBUG1
        printVec(newBasis, power);
#endif
        for (int j = 0; j < n; ++j) {
            //calculating the new order of the columns
            int temp = 0;
            for (int k = 0; k < power; ++k) {
                temp ^= (j & (1 << k)) ? newBasis[k] : 0;
            }
            tempArray[j] = temp;
#ifdef DEBUG1
            std::cout << temp << ' ';
#endif
            for (int k = 0; k < n; ++k) {
                tempMatrix[k * n + j] = matrix[k * n + temp];
            }
        }
#ifdef DEBUG1
        std::cout << std::endl;
        printMatrix(tempMatrix, n, n);
        std::cout << std::endl;
#endif
        auto tempKernel = new CMatrixBinaryKernel(length, tempMatrix, "name", nullptr);
        auto processor = (SectionedTrellisKernelProcessor*)(tempKernel->GetProcessor(0, sectionSize, codeLength));
        //auto processorState = processor->GetState(1);
        //auto processorTemp = processor->GetTemp(1);
        //auto processorMetrics = processor->GetMetrics();
        sum = SumCount;
        cmp = CmpCount;


        for (int j = 0; j < n; ++j) {
            processor->GetLLRs(1, j, pKnownInputSymbols, pChannelLLRs, pLLRs, processorState, processorTemp, processorMetrics);
        }
#ifdef DEBUG1
        printVec(oldBasis, power);
#endif
        if ((SumCount - sum) < minSum && (CmpCount - cmp) < minCmp) {
            minSum = SumCount - sum;
            minCmp = CmpCount - cmp;
            std::memcpy(newKernel, tempMatrix, n * n * sizeof(tBit));
            std::memcpy(finalTempArray, tempArray, n * sizeof(int));
            std::memcpy(finalBasis, newBasis, power * sizeof(unsigned long));
        }
#ifdef DEBUG1
        printVec(oldBasis, power);
#endif
        //processor->FreeState(processorState);
        //processor->FreeTemp(processorTemp);
        //processor->FreeMetrics(processorMetrics);
        delete processor;
        delete tempKernel;
    }
    delete[] pKnownInputSymbols;
    delete[] pLLRs;
    delete[] pChannelLLRs;
    tBit* S = (tBit*)processorState;
    AlignedFree(S);
    free(processorTemp);
    for (int i = 0; i < size1; ++i) {
        for (int j = 0; j < pow(n, i); ++j) {
            for (int k = 0; k < L - 1; ++k) {
                delete[] processorMetrics[i][j][k];
            }
            delete[] processorMetrics[i][j];
        }
        delete[] processorMetrics[i];
    }
    delete[] processorMetrics;
#ifdef DEBUG
    printVec(finalBasis, power);
    printVec(finalTempArray, n);
    printMatrix(newKernel, n, n);
#endif
    std::cout << minSum << " " << minCmp << "\n";
    delete[] oldBasis;
    delete[] newBasis;
    delete[] b;
    delete[] l;
    delete[] u;
    delete[] tempMatrix;
    delete[] tempArray;
    delete[] finalTempArray;
    SumCount = prevSum;
    CmpCount = prevCmp;
}


/*void getMatr(unsigned char *mas, unsigned char *p, int i, int j, int m) {
    int ki, kj, di, dj;
    di = 0;
    for (ki = 0; ki < m - 1; ki++) {
        if (ki == i) di = 1;
        dj = 0;
        for (kj = 0; kj < m - 1; kj++) {
            if (kj == j) dj = 1;
            p[ki * m + kj] = mas[(ki + di) * m + j + dj];
        }
    }
}


int determinant(unsigned char *mas, int m) {
    int i, j, d, k, n;
    unsigned char *p;
    p = new unsigned char[m];
    j = 0; d = 0;
    k = 1;
    n = m - 1;
    if (m == 1) {
        d = mas[0];
        return(d);
    }
    if (m == 2) {
        d = mas[0] * mas[m + 1] - (mas[m] * mas[1]);
        return(d);
    }
    if (m>2) {
        for (i = 0; i<m; i++) {
            getMatr(mas, p, i, 0, m);
            d = d + k * mas[i * m] * determinant(p, n);
            k = -k;
        }
    }
    return(d);
}


bool increase(int power, unsigned char* b) {
    int n = 1 << power;
    char transfer, old;
    if (b[0]) {
        transfer = 1;
        b[0] = 0;
    }
    else {
        transfer = 0;
        b[0] = 1;
    }
    for (int i = 0; i < power; ++i) {
        for (int j = 0; j < power; ++j) {
            if (!i && !j) {
                old = b[i * n + j];
                b[i * n + j] ^= transfer;
                transfer &= old;
            }
        }
    }
    return determinant(b, power) != 0;
}*/


void randomInvertibleMatrix(int power, unsigned char* l, unsigned char* u, unsigned char* b) {
    //filling in matrices l and u with random values
    for (int i = 0; i < power; ++i) {
        for (int j = 0; j < i; ++j) {
            l[i * power + j] = distribution(generator);
        }
        for (int j = i + 1; j < power; ++j) {
            u[i * power + j] = distribution(generator);
        }
    }
    //multiplying matrices l and u and writing the result in b
    for (int i = 0; i < power; ++i) {
        for (int j = 0; j < power; ++j) {
            b[i * power + j] = 0;
            for (int k = 0; k < power; ++k) {
                b[i * power + j] ^= l[i * power + k] & u[k * power + j];
            }
        }
    }
}