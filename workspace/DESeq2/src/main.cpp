//============================================================================
// Name        : DESeq2.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <armadillo>
using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
    mat A = randu<mat>(4,5);
    mat B = randu<mat>(4,5);

    cout << A*trans(B) << endl;

    return 0;
}

