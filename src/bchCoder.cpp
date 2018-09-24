//
// Created by Moskovskaya Elizaveta on 10.05.18.
//

#include <random>
#include <iostream>
#include <fstream>
#include <cstring>
#include <chrono>
#include "../headers/bchCoder.h"

auto seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
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

void multiplyPolynomials(const unsigned char* first, int size1, const unsigned char* second, int size2, unsigned char *res, int *sizeRes) {
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


void printMatrix(unsigned char** const matrix, int sizeI, int sizeJ) {
    if (sizeJ == -1) sizeJ = sizeI;
    for (auto i = 0; i < sizeI; ++i) {
        for (auto j = 0; j < sizeJ; ++j) {
            std::cout << ((matrix[i][j]) ? 1 : 0) << ' ';
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


void makeMatrix(int power, const unsigned long * fieldElements, unsigned char** matrix) {
    const auto amount = ((1 << power) - 2) / 2;
    const auto len = (1 << power) - 1;
    int polySize, gSizeOld = 1, gSizeNew, size;
    auto poly = new unsigned char[power + 1];
    auto g = new unsigned char[len];
    g[0] = 1;
    for (auto i = 0; i < len; ++i) {
        for (auto j = 0; j < len; ++j) {
            matrix[i][j] = 0;
        }
    }
    matrix[0][0] = 1;
    for (unsigned long i = 2; i <= amount; i++) {
        findMinimalPolynomial(i, power, fieldElements, &polySize, poly);
        if (gSizeOld >= polySize && !dividePolynomial(g, gSizeOld, poly, polySize, &size, true)[0] && size == 1) { continue; }
        gSizeNew = polySize + gSizeOld - 1;
        multiplyPolynomials(poly, polySize, g, gSizeOld, matrix[gSizeNew - 1]);
        auto count = 1;
        for (auto j = gSizeOld; j < gSizeNew - 1; ++j) {
            for (auto k = count; k < count + gSizeOld; ++k) {
                matrix[j][k] = g[k - count];
            }
            ++count;
        }
        memcpy(g, matrix[gSizeNew - 1], gSizeNew);
        gSizeOld = gSizeNew;
    }
}