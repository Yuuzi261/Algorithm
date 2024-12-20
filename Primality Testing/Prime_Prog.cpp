#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>
#include <fstream>

using namespace std;

// const
int FERMAT_K = 5, MILLER_RABIN_K = 5;

// algorithms
bool BasicPrimeTesting(long long);
bool FermatPrimalityTesting(long long);
bool MillerRabin_PrimalityTesting(long long);
vector<bool> SieveOfEratosthenes(long long);
vector<bool> BasicPrimeTestingToN(long long);
vector<bool> FermatPrimalityTestingToN(long long);
vector<bool> MillerRabin_PrimalityTestingToN(long long);

// tools
long long fast_pow_mod(long long, long long, long long);
void recordTime(long long, long double, string);
void outputTime(long long, long double);

// for time calculation
void FunctionExecutionTimeTesting(long long, bool (*func) (long long), string);
void FunctionExecutionTimeTesting(long long, bool (*func) (long long), int, string);
void FunctionExecutionTimeTesting(long long, vector<bool> (*func) (long long), string);

// for testing
void SingleNumberTestToN(long long);
void SegmentalTesting(long long, int);

// output
void outputIsPrime(long long, bool);

int main(void) {

    long long natural;

    srand(time(NULL));

    // test the functions you want

    return 0;
}

// algorithms
bool BasicPrimeTesting(long long n) {

    for (int i = 2;i*i <= n;i++)
        if (n % i == 0) return false;
    return true;

}

bool FermatPrimalityTesting(long long n) {

    long long x;
    
    for(int i = 0;i < FERMAT_K;i++) {
        x = rand()%(n-1) + 1;
        if(fast_pow_mod(x, n-1, n) != 1) return false;
    }

    return true;

}

bool MillerRabin_PrimalityTesting(long long n) {

    if(n < 2) return false;
    if(n <= 3) return true;
    if(n % 2 == 0) return false;

    long long a, x;
    long long r = 0, d = n - 1;

    while(d % 2 == 0) {
        r++;
        d /= 2;
    }

    for(int i = 0;i < MILLER_RABIN_K;i++) {

        a = rand()%(n-3) + 2;
        x = fast_pow_mod(a, d, n);

        if(x == 1 || x == n-1) continue;

        for(int j=1;j < r;j++) {
            x = fast_pow_mod(x, 2, n);
            if(x == 1) return false;
            if(x == n-1) break;
        }
        if (x == n - 1) continue;
        return false;

    }

    return true;

}

vector<bool> SieveOfEratosthenes(long long n) {

    vector<bool> isPrimeVector;
    isPrimeVector.resize(n+1);

    isPrimeVector[0] = true, isPrimeVector[1] = true;
    for(int i=2;i<sqrt(n);i++) {
        if(isPrimeVector[i] == false) {
            for(int j=2;j<n;j++) {
                if(i * j > n) break;
                isPrimeVector[i*j] = true;
            }
        }
    }

    return isPrimeVector;

}

vector<bool> BasicPrimeTestingToN(long long n) {

    vector<bool> isPrimeVector;
    isPrimeVector.resize(n+1);

    isPrimeVector[0] = false, isPrimeVector[1] = false;
    for(int i = 2;i < n;i++)
        isPrimeVector[i] = BasicPrimeTesting(i);

    return isPrimeVector;

}

vector<bool> FermatPrimalityTestingToN(long long n) {

    vector<bool> isPrimeVector;
    isPrimeVector.resize(n+1);

    isPrimeVector[0] = false, isPrimeVector[1] = false;
    for(int i = 2;i < n;i++)
        isPrimeVector[i] = FermatPrimalityTesting(i);

    return isPrimeVector;

}

vector<bool> MillerRabin_PrimalityTestingToN(long long n) {

    vector<bool> isPrimeVector;
    isPrimeVector.resize(n+1);

    isPrimeVector[0] = false, isPrimeVector[1] = false;
    for(int i = 2;i < n;i++)
        isPrimeVector[i] = MillerRabin_PrimalityTesting(i);

    return isPrimeVector;

}

// tools
long long fast_pow_mod(long long a, long long b, long long n) {
    long long ans = 1;
    while (b > 0) {
        if (b & 1) ans = ans * a % n;
        a = a * a % n;
        b >>= 1;
    }
    return ans;
}

void recordTime(long long number, long double elapsed, string outputFile) {

    ofstream ofs;

    ofs.open(outputFile, std::ios_base::app);
    ofs << number << ' ' << elapsed / 1000000.0 << ' ' << elapsed << endl;
    ofs.close();

}

void outputTime(long long n, long double elapsed) {

    cout.setf(ios::fixed);
    cout.precision(3);
    cout << "Time: " << elapsed / 1000000.0 << 's';
    cout.precision(2);
    cout << " (" << elapsed << "ms)" << endl;
    // cout << n << endl;

}

// for time calculation
void FunctionExecutionTimeTesting(long long n, bool (*func) (long long), string outputFile = NULL) {

    long double elapsed;
    struct timespec start, end;
    bool isPrime;

    clock_gettime(CLOCK_REALTIME, &start);

    isPrime = (*func) (n);

    clock_gettime(CLOCK_REALTIME, &end);
    elapsed = (end.tv_sec - start.tv_sec) * 1000000.0 + (end.tv_nsec - start.tv_nsec) / 1000.0;

    if(!outputFile.empty()) recordTime(n, elapsed, outputFile);
    else outputTime(n, elapsed);

}

void FunctionExecutionTimeTesting(long long n, bool (*func) (long long), int interval, string outputFile = NULL) {

    long double elapsed;
    struct timespec start, end;
    bool isPrime;

    clock_gettime(CLOCK_REALTIME, &start);

    for(long long i = n;i < n + interval;i++) 
        isPrime = (*func) (n);

    clock_gettime(CLOCK_REALTIME, &end);
    elapsed = (end.tv_sec - start.tv_sec) * 1000000.0 + (end.tv_nsec - start.tv_nsec) / 1000.0;

    if(!outputFile.empty()) recordTime(n, elapsed, outputFile);
    else outputTime(n, elapsed);

}

void FunctionExecutionTimeTesting(long long n, vector<bool> (*func) (long long), string outputFile = NULL) {

    long double elapsed;
    struct timespec start, end;
    vector<bool> isPrimeVector;
    isPrimeVector.resize(n);

    clock_gettime(CLOCK_REALTIME, &start);

    isPrimeVector = (*func) (n);

    clock_gettime(CLOCK_REALTIME, &end);
    elapsed = (end.tv_sec - start.tv_sec) * 1000000.0 + (end.tv_nsec - start.tv_nsec) / 1000.0;

    if(!outputFile.empty()) recordTime(n, elapsed, outputFile);
    else outputTime(n, elapsed);

}

// for testing
void SingleNumberTestToN(long long n) {

    for(int i=2;i < n;i++) {
        FunctionExecutionTimeTesting(i, BasicPrimeTesting, "test.txt");
        FunctionExecutionTimeTesting(i, SieveOfEratosthenes, "test.txt");
        FunctionExecutionTimeTesting(i, FermatPrimalityTesting, "test.txt");
        if(i % 10000 == 0) cout << "dealing with " << i << "..." << endl;
    }

}

void SegmentalTesting(long long n, int interval) {

    long long index = 2;

    while(index < n) {
        cout << "dealing with " << index << "..." << endl;
        FunctionExecutionTimeTesting(index, BasicPrimeTesting, interval, "speed_interval.txt");
        FunctionExecutionTimeTesting(index, FermatPrimalityTesting, interval, "speed_interval.txt");
        index += interval;
    }

}

// output
void outputIsPrime(long long n, bool isPrime) {

    if(isPrime) cout << n << " is prime." << endl;
    else cout << n << " is not prime." << endl;

}