/* Sieve of Pritchard in C++, as described at https://en.wikipedia.org/wiki/Sieve_of_Pritchard
   arguments: N [-p]
      N: finds primes up to N
     -p: (optional) print the primes found
   simple implementation using a simple array of integers and a bit array (using bytes for speed)
   2 <= N <= 1000000000
   (like the standard Sieve of Eratosthenes, this algorithm is not suitable for very large N due to memory requirements) */

// compile with -DLOGGING to print array after each main step

#include <cstring>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctime>

double e_time = 0.0; // Extend
double d_time = 0.0; // Delete
double c_time = 0.0; // Compress
double p_time = 0.0; // Primes (counting)

void check(uint32_t w[], uint32_t from, uint32_t to, bool d[]) { // for logging array contents
    for (uint32_t i=from; i <= to; i++) printf(" w[%u]=%u%s", i, w[i], (d[w[i]] ? "*" : ""));
    printf("\n");
}

void Extend (uint32_t w[], uint32_t &w_end, uint32_t &length, uint32_t n, bool d[]) {
    // Rolls full wheel W up to n, and sets length=n
    uint32_t i, j, x;
    int start_s = clock();
    i = 0; j = w_end;
    x = length + 1; // length+w[0]
    while (x <= n) {
        w[++j] = x; // Append x to the ordered set W
        d[x] = false;
        x = length + w[++i];
    }
    length = n; w_end = j;
    int stop_s=clock();
    e_time += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    #if LOGGING
    printf("After Extend(%u):", n); check(w, 0, w_end, d);
    #endif
}

void Delete (uint32_t w[], uint32_t length, uint32_t p, bool d[], uint32_t &imaxf) {
    // Deletes multiples p*w[i] of p from W, and sets imaxf to last i deleted
    uint32_t i, x;
    int start_s = clock();
    i = 0;
    x = p; // p*w[0]=p*1
    while (x <= length) {
        d[x] = true; // Remove x from W;
        x = p*w[++i];
    }
    imaxf = i-1;
    int stop_s=clock();
    d_time += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    #if LOGGING
    printf("After Delete(%u): (partial)", p); check(w, 0, imaxf, d);
    #endif
}

void Compress(uint32_t w[], bool d[], uint32_t to, uint32_t &w_end) {
    // Removes deleted values in w[0..to], and if to=w_end, updates w_end, otherwise pads with zeros on right
    uint32_t i, j;
    int start_s = clock();
    j = 0;
    for (i=1; i <= to; i++) if (!d[w[i]]) w[++j] = w[i];
    if (to == w_end) w_end = j; else for (uint32_t k=j+1; k <= to; k++) w[k] = 0;
    int stop_s=clock();
    c_time += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    #if LOGGING
    printf("After Compress(%u):", to); check(w, 0, w_end, d);
    #endif
}

void Sift(uint32_t N, bool printPrimes, uint32_t &nrPrimes) {
    // finds the nrPrimes primes up to N, printing them if printPrimes
    uint32_t *w = new uint32_t[N/4+5];
    bool *d = new bool[N+1];
    uint32_t w_end, length;
    /* representation invariant (for the main loop):
       if length < N (so W is a complete wheel), w[0..w_end] is the ordered set W;
       otherwise, w[0..w_end], omitting zeros and values w with d[w] true, is the ordered set W,
       and no values are omitted before a value between 1 and N/p */
    uint32_t p, imaxf;
    // W,k,length = {1},1,2:
    w_end = 0; w[0] = 1; length = 2;
    // Pr = {2}:
    nrPrimes = 1; if (printPrimes) printf("%u\n", 2);
    p = 3;
    /* invariant: p = p_(k+1) and W = W_k inter {1,...,N} and length = min(P_k,N) and Pr = the first k primes
       (where p_i denotes the i'th prime, W_i denotes the i'th wheel, P_i denotes the product of the first i primes) */
    while (p*p <= N) {
        // Append p to Pr:
        nrPrimes++; if (printPrimes) printf("%u\n", p);
        if (length < N) Extend (w, w_end, length, std::min(p*length,N), d);
        Delete(w, length, p, d, imaxf);
        if (length < N) {
            Compress(w, d, w_end, w_end);
        } else {
            while (d[w[imaxf]]) imaxf++; // get sentinel for Delete loop
            Compress(w, d, imaxf, w_end);
        }
        // p = next(W, 1):
        p = w[1];
        // k++
    }
    if (length < N) {
        // Extend full wheel W,length to N:
        Extend (w, w_end, length, N, d);
    }
    // gather remaining primes:
    int start_s = clock();
    for (uint32_t i=1; i <= w_end; i++) {
        if (w[i] == 0 || d[w[i]]) continue;
        nrPrimes++; if (printPrimes) printf("%u\n", w[i]);
    }
    delete[] w; delete[] d;
    int stop_s=clock();
    p_time += (stop_s-start_s)/double(CLOCKS_PER_SEC);
}

int main (int argc, char *argw[]) {
    bool error = false, printPrimes = false;
    uint64_t max = 1000000000, N;
    uint32_t nrPrimes;
    if (argc == 3) {
        if (strcmp(argw[2], "-p") == 0) {
            printPrimes = true;
            argc--;
        } else {
            error = true;
        }
    }
    if (argc == 2) { N = strtol(argw[1], NULL, 10); if (N < 2 || N > max) error = true; }
    if (error) { printf("call with: %s N -p where 2 <= N <= %lu and -p to print the primes is optional \n", argw[0], max); exit(1); }
    int start_s = clock();
    Sift(N, printPrimes, nrPrimes);
    int stop_s=clock();
    printf("%u primes up to %lu found in %.3f ms\n", nrPrimes, N, (stop_s-start_s)*1E3/double(CLOCKS_PER_SEC));
    printf("Time for: Extend=%.2f ms, Delete=%.2f ms, Compress=%.2f ms, Primes=%.2f ms\n", e_time*1E3, d_time*1E3, c_time*1E3, p_time*1E3);
}