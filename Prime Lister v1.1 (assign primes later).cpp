#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <sstream>
#include <cstdio>
#include <algorithm>
#include <string>
#include <iterator>
#include <cstdlib> 
#include <ctime> 
#include <list>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <cctype>
#include <iomanip>
#include "omp.h"


using std::cerr;
using std::endl;
using std::ofstream;
using namespace std;

// Include all that common stuff **************

int main()
{

    unsigned long long int i;
    unsigned long long int j;
    unsigned long long int count;
    unsigned long long int obergrenze;
    unsigned long long int bereich;
    unsigned long long int breite;
    ofstream outdata; // outdata is like cin
    unsigned long long int limiter;
    unsigned long long int step;

    cout << "*****************************************\n";
    cout << "  PrimeLister v1.1 (2020)                *\n";
    cout << "******************************************* \n\n";
    cout << " Please specify the upper limit as number [N]. \n\n";
    cout << " With 16 GB of RAM a limit of 4.500.000.000 is recommended for [N]. \n\n";
    cout << " With 32 GB of RAM you may try 10.000.000.000, \n\n";
    cout << " With 64 GB of RAM even 20.000.000.000 can be possible, \n\n";
    cout << " With 128 GB of RAM you need to get laid soon... \n\n";
    cout << " Software may become unstable at the memory limit of your machine.\n\n";

    cout << " Please enter an integer value to calculate Primes below [N, should start from 1 million...]: ";
    cin >> obergrenze;

    cout << "\n Please enter the number of columns (integer) for the output table (PrimeList.txt): ";
    cin >> breite;
    if (breite < 1) { breite = 1; }


    const clock_t t1 = clock();
    bereich = roundl(obergrenze / 3) + (15 * breite);

    cout << "\n Fill a Vector with ones...\n";

    std::vector<unsigned long long int>A(bereich, 1);

    cout << " Delete all the multiples of the single elements (E0, E2, E4...)\n";

    for (i = 0; (7 * i + 10) < A.size(); i = i + 2) {
        while (A[i] == 0) { i = i + 2; }
        step = 6 * i + 10;
        for (j = 7 * i + 10; j < A.size(); j = j + step) { A[j] = 0; }
    }

    cout << " Delete all the multiples of the single elements (E1, E3, E5...)\n";


    for (i = 1; (3 * i * i + 8 * i + 4) < A.size(); i = i + 2) {
        while (A[i] == 0) { i = i + 2; }
        step = 6 * i + 8;
        for (j = 3 * i * i + 8 * i + 4; j < A.size(); j = j + step) { A[j] = 0; }
    }

    cout << " Now delete all the remaining squares of the single elements (E0, E2, E4...)\n";


    for (i = 0; (3 * i * i + 10 * i + 7) < A.size(); i = i + 2) {
        if (A[3 * i * i + 10 * i + 7] == 1) { step = 6 * i + 10; A[3 * i * i + 10 * i + 7] = 0; for (j = (3 * i * i + 10 * i + 7); j < A.size(); j = j + step) { A[j] = 0; } }
    }

    cout << " Replace all the ones by the according Primes\n";

    for (i = 0; i < A.size(); i = i + 2) {
        while (A[i] == 0 && i < A.size()) { i = i + 2; }
        A[i] = 3 * i + 5;
    }

    for (i = 1; i < A.size(); i = i + 2) {
        while (A[i] == 0 && i < A.size()) { i = i + 2; }
        A[i] = 3 * i + 4;
    }

    cout << " And delete all zeros\n";
    A.erase(remove(A.begin(), A.end(), 0), A.end());

    const clock_t t2 = clock();

    cout << " Done!\n";

    cout << " Insert the missing Primes 2 and 3 \n";

    A.insert(A.begin(), 3);
    A.insert(A.begin(), 2);

    cout << "\n That's the first 100 Primes generated: \n \n";

    for (i = 0; i < 100; i++) {
        std::cout << A[i] << ' ';
    }

    cout << "\n \n That's the last 100 Primes generated: \n \n";

    for (i = (A.size() - 100); i < A.size(); i++) {
        std::cout << A[i] << ' ';
    }

    cout << "\n \n Overall Primes found: ";
    std::cout << A.size() << '\n';

    float TO = (t2 - t1);

    cout << "\n Overall calculation time [Sec] : " << TO * 0.001 << "\n";

    int power = (A.size() / (TO * 0.001));

    cout << " => About " << power << " Primes per second generated \n\n";

	cout << "\n The Primes (all!) are now written to: 'PrimeList.txt' (that may take a few secons) \n";

    outdata.open("PrimeList.txt");
    if (!outdata) {
        cerr << " Error: File could not be opened!" << endl;
        exit(1);
    }

    if (breite == 1) {

        for (i = 0; i < A.size(); i++) {
            outdata << (i + 1) << ":\t" << A[i] << "\n";
        }
    }

    if (breite > 1) {

        count = 0;

        limiter = A.size() % breite;

        outdata << "1" << "-" << breite << "\t";

        for (i = 0; i < (A.size() - limiter); i++) {
            outdata << A[i] << "\t";
            count++;
            if (count == breite && i != (A.size() - limiter - 1)) { count = 0; outdata << "\n" << (i + 2) << "-" << (i + 1 + breite) << "\t"; }
        }
    }


    cout << "\n\n Close this Window to EXIT... ";
    outdata.close();

    int end;

    cin >> end;

}

