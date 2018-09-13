//
// Created by Moskovskaya Elizaveta on 14/06/2018.
//

#include <string>
#include <fstream>
#include <cmath>
#include <iostream>
#include <time.h>
#include "../headers/KanekoKernelProcessor.h"
#include "../headers/bchCoder.h"

#define DEBUG

void fun(const std::string& file, KanekoKernelProcessor& decoder, const unsigned char *g, unsigned long gSize, long p, long e, double maxSTNR = 5.0) {

    setlocale(LC_ALL, "Russian");

    int count = 0, countErr = 0;
    double standardDiviation;
    auto info = new unsigned char[decoder.getK()];
    auto res = new unsigned char[decoder.getN()];
    auto err = new double[decoder.getN()];
    auto decoded = new unsigned char[decoder.getN()];
    unsigned long wordCount = 0;

    std::ofstream fout;
    fout.open(file + ".csv");

    /*std::ofstream foutT;
    foutT.open(file + "_t.csv");

    std::ofstream foutC;
    foutC.open(file + "_c.csv");

    std::ofstream foutS;
    foutS.open(file + "_s.csv");*/

    clock_t start, end;
    start = clock();

#ifdef DEBUG
    std::ofstream ferr;
    ferr.open("out/errWords.txt");
#endif

    for (double stnr = 0.0; stnr <= maxSTNR; stnr += 0.5) {

        while (count < p && countErr < e) {

            standardDiviation = sqrt(1 / (pow(10, stnr / 10) * 2 * decoder.getK() / decoder.getN()));
            generateRandomPoly(info, decoder.getK());
            multiplyPolynomials(info, decoder.getK(), g, gSize, res);
            addNoise(standardDiviation, res, err, decoder.getN());


#ifndef DEBUG
            decoder.decode(err, decoded);
#endif

#ifdef DEBUG
            decoder.set(err);
            decoder.decode(decoded);
            if (decoder.calcL(res) < decoder.calcL(decoded)) {
                ferr << stnr << "\n";
                printVec(ferr, res, decoder.getN());
                printVec(ferr, err, decoder.getN());
                ferr << "\n";
            }
#endif

            if (!comparePoly(res, decoder.getN(), decoded, decoder.getN())) {
                ++countErr;
            }
            ++count;
            ++wordCount;
        }

        fout << stnr << "," << ((double)countErr)/count
                     << "," << ((double) decoder.getDecodingCount()) / wordCount
                     << "," << ((double) decoder.getComparisonCount()) / wordCount
                     << "," << ((double) decoder.getSummCount()) / wordCount << "\n";

        /*fout << stnr << "," << ((double)countErr)/count << "\n";
        foutT << stnr << "," << ((double) decoder.getDecodingCount()) / wordCount << "\n";
        foutC << stnr << "," << ((double) decoder.getComparisonCount()) / wordCount << "\n";
        foutS << stnr << "," << ((double) decoder.getSummCount()) / wordCount << "\n";*/

        decoder.setDecodingCount();
        decoder.setComparisonCount();
        decoder.setSummCount();
        wordCount = 0;
        count = 0, countErr = 0;

        std::cout << stnr << "\n";

#ifdef DEBUG
        ferr << "\n";
#endif

    }

    end = clock();

    std::cout << "Общее время: " << ((double) end - start) / (double) CLOCKS_PER_SEC << " секунд\n";

    fout.close();
    /*foutT.close();
    foutC.close();
    foutS.close();*/

    delete[] info;
    delete[] res;
    delete[] err;
    delete[] decoded;

}