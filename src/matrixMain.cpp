//
// Created by Moskovskaya Elizaveta on 22/07/2018.
//


#include <iostream>
#include <cmath>
#include <fstream>
#include "../headers/bchCoder.h"

unsigned long primitivePolynomials[16] =
        {3, 7, 11, 19, 37, 67, 137, 285, 529, 1033, 2053, 4179, 8219, 17475, 32771, 69643};

int main (int argc, char* argv[]) {
    try {
        long power;
        std::string filename;
        if (argc == 2) {
            power = std::atoi(argv[1]);
        }
        else if (argc == 3) {
            power = std::atoi(argv[1]);
            filename = argv[2];
        }
        else {
            throw "Invalid argument(s)\n";
        }
        if (power <= 1) {
            throw "Invalid values of argument(s)\n";
        }

        //setting up the field
        int n = (1 << power) - 1;
        unsigned long primitivePolynomial = primitivePolynomials[power - 1];
        auto *fieldElements = new unsigned long[n];
        fieldElements[0] = 1;
        fieldElements[1] = 1 << 1;
        for (int i = 2; i < n; ++i) {
            fieldElements[i] = fieldElements[i - 1] << 1;
            if (fieldElements[i] >> power == 1) {
                fieldElements[i] ^= primitivePolynomial;
            }
        }

        //setting up the coder
        auto matrix = new unsigned char*[n];
        for (int i = 0; i < n; ++i) {
            matrix[i] = new unsigned char[n];
        }
        makeMatrix(power, fieldElements, matrix);
        if (filename != "") {
            std::ofstream out(filename);
            printMatrix(out, matrix, n);
            out.close();
        }
        else {
            printMatrix(matrix, n);
        }
    }
    catch (const char *err) {
        std::cerr << err;
    }
    return 0;
}
