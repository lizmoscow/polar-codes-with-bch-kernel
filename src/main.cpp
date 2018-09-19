//
// Created by Moskovskaya Elizaveta on 11.05.18.
//


#include <iostream>
#include <cmath>
#include <climits>
#include <fstream>
#include "../headers/bchCoder.h"
#include "../headers/KanekoKernelProcessor.h"
#include "../headers/dataForPlot.h"

const unsigned long fieldElements[16] =
        {3, 7, 11, 19, 37, 67, 137, 285, 529, 1033, 2053, 4179, 8219, 17475, 32771, 69643};

int main (int argc, char* argv[]) {

    try {

        long power;
        long t;
        std::string filename;
        long p = 1, e = 1;
        double signalToNoiseRatio = 0;

        //random words with particular noise
        if (argc == 4) {
            power = std::atoi(argv[1]);
            t = std::atoi(argv[2]);
            signalToNoiseRatio = std::atof(argv[3]);
        }

        //word from file with particular noise
        else if (argc == 5) {
            power = std::atoi(argv[1]);
            t = std::atoi(argv[2]);
            signalToNoiseRatio = std::atof(argv[3]);
            filename = argv[4];
        }

        //plot computing
        else if (argc == 6) {
            power = std::atoi(argv[1]);
            t = std::atoi(argv[2]);
            filename = argv[3];
            p = std::atoi(argv[4]);
            e = std::atoi(argv[5]);
        }

        else {
            throw "Invalid arguments\n";
        }

        if (t <= 0 || t >= (1 << power - 1)|| power <= 1 || p <= 0 || e <= 0 || signalToNoiseRatio < 0) {
            throw "Invalid values of arguments\n";
        }

        //setting up the field
        int n = (1 << power) - 1;
        unsigned long primitivePolynomial = fieldElements[power - 1];

        auto *antilogarithms = new unsigned long[n];
        antilogarithms[0] = 1;
        antilogarithms[1] = 1 << 1;

        auto *logarithms = new unsigned long[n + 1];
        logarithms[0] = LONG_MAX;
        logarithms[1] = 0;
        logarithms[2] = 1;

        for (unsigned long i = 2; i < n; ++i) {
            antilogarithms[i] = antilogarithms[i - 1] << 1;
            if (antilogarithms[i] >> power == 1) {
                antilogarithms[i] ^= primitivePolynomial;
            }
            logarithms[antilogarithms[i]] = i;
        }

        //setting up the coder
        int gSize;
        auto g = new unsigned char[(1 << power) - 1];
        unsigned char *temp;
        findMinimalPolynomial(1, power, antilogarithms, &gSize, g);
        auto minimalPolynomial = new unsigned char[power + 1];
        int mpSize;
        for (unsigned long i = 2; i < 2 * t; ++i) {
            findMinimalPolynomial(i, power, antilogarithms, &mpSize, minimalPolynomial);
            temp = lcm(g, gSize, minimalPolynomial, mpSize, &gSize);
            delete[] g;
            g = temp;
        }
        int k = n - gSize + 1;
        int d = 2 * t + 1;
        printVec(g, gSize);

        std::cout << "(" << n << ", " << k << ", " << d << ")\n";

        //random words with particular noise
        if (argc == 4) {
            //generating a word
            auto info = generateRandomPoly(k);
            auto res = new unsigned char[n];
            multiplyPolynomials(info, k, g, gSize, res);
            printVec(res, n);

            //adding up noise
            double sd = sqrt(1 / (pow(10, signalToNoiseRatio / 10) * 2 * k / n));
            auto err = new double[n];
            addNoise(sd, res, err, n);
            printVec(err, n);

            //decoding
            KanekoKernelProcessor decoder(power, n, t, k, antilogarithms, logarithms, 0.5);
            auto decoded = new unsigned char[n];
            decoder.decode(err, decoded);

            //checking the result
            printVec(decoded, n);
            if (!comparePoly(res, n, decoded, n)) {
                std::cout << "Errors were not corrected!\n";
            }
            else {
                std::cout <<"Ok\n";
            }

            delete[] info;
            delete[] res;
            delete[] err;
            delete[] decoded;
        }

        //word from file with particular noise
        else if (argc == 5) {
            //reading a word from file
            std::ifstream in(filename);
            auto res = new unsigned char[n];
            auto err = new double[n];
            if (in.is_open()) {
                for (int i = 0; i < n; ++i) {
                    in >> res[i];
                    res[i] = (res[i] == '1') ? res[i] = 1 : res[i] = 0;
                }
                printVec(res, n);
                for (int i = 0; i < n; ++i) {
                    in >> err[i];
                }
                printVec(err, n);

            }
            else {
                throw("File does not exsist!");
            }

            //decoding
            KanekoKernelProcessor decoder(power, n, t, k, antilogarithms, logarithms, 0.5);
            auto decoded = new unsigned char[n];
            decoder.decode(err, decoded);

            //checking the results
            printVec(decoded, n);
            if (!comparePoly(res, n, decoded, n)) {
                std::cout << "Errors were not corrected!\n";
            }
            else {
                std::cout <<"Ok\n";
            }

            delete[] res;
            delete[] err;
            delete[] decoded;
        }

        //plot computing
        else if (argc == 6) {
            KanekoKernelProcessor decoder(power, n, t, k, antilogarithms, logarithms, 0.5);
            fun(filename, decoder, g, gSize, p, e, 5.0);
        }

        delete[] antilogarithms;
        delete[] g;
        delete[] minimalPolynomial;


    }
    catch (const char *err) {
        std::cerr << err;
    }
    return 0;
}


//-1.54085 3.43615 2.65674 -0.793389 -1.91762 0.428461 -0.893624 0.454996 -0.464513 -1.48376 0.67617 -1.75345 1.8906 -0.252364 1.40103
