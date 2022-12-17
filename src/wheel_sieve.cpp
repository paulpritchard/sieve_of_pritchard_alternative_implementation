/* Sieve of Pritchard in C++, as described at https://en.wikipedia.org/wiki/Sieve_of_Pritchard
   arguments: N [-p]
      N: finds primes up to N
     -p: (optional) print the primes found
   efficient implementation using a simple array of integers and a bit array
   2 <= N <= 1000000000
   (like the standard Sieve of Eratosthenes, this algorithm is not suitable for very large N due to memory requirements) */

// compile with -O3 -march=native -lm
// add -DLOGGING to print array after each main step

#include <cstring>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctime>
#include <math.h>
#include <x86intrin.h>

#define delete(x,d) d[x / 8] |= 1 << (x % 8)
#define deleted(x,d) (d[x / 8] & 1 << (x % 8))

double e_time = 0.0; // Extend
double d_time = 0.0; // Delete
double c_time = 0.0; // Compress
double p_time = 0.0; // Primes (counting)

void check (uint32_t* w, uint32_t from, uint32_t to, char* d) { // for logging array contents
    for (uint32_t i=from; i <= to; i++) printf(" w[%u]=%u%s", i, w[i], (deleted(w[i],d) ? "*" : ""));
    printf("\n");
}

uint32_t find (uint32_t* w, uint32_t bound, uint32_t low, uint32_t high) {
    // returns max{i: w[i] <= bound} given w[low] <= bound < w[high]
    uint32_t i = low, j = high, mid;
    /* invariant: w[i] <= bound < w[j] */
    while (i+1 < j) if (w[mid = i+(j-i)/2] <= bound) i = mid; else j = mid;
    return(i);
}

void Extend (uint32_t* &w, uint32_t &w_end, uint32_t &length, uint32_t n, char* d) {
    // Rolls full wheel W up to n, and sets length=n
    int start_s = clock();
    // get bigger array
    uint32_t nr_complete = (n/length)*(w_end+1); // from full copies of current wheel
    uint32_t rem = n%length;
    uint32_t nr_extra = 0; // from partial current wheel
    if (rem > 0) nr_extra = find(w, rem, 0, w_end+1)+1;
    uint32_t *w_next = new uint32_t[nr_complete+nr_extra+1];
    memcpy(w_next, w, (w_end+1)*sizeof(uint32_t)); // full copy of current wheel
    delete[] w; w = w_next;
    uint32_t offset = 0; uint32_t* dest = w + w_end+1;
    for (uint32_t i=1; i < n/length; i++) {
        memcpy(dest, w, (w_end+1)*sizeof(uint32_t)); // full copy of current wheel
        offset += length;
        for (uint32_t* end = dest + w_end+1; dest < end; dest++) *dest += offset; // (should vectorize)
    }
    memcpy(dest, w, nr_extra*sizeof(uint32_t)); // partial copy of current wheel
    offset += length;
    for (uint32_t* end = dest + nr_extra; dest < end; dest++) *dest += offset; // (should vectorize)
    w[nr_complete+nr_extra] = n; // sentinel
    length = n; w_end = nr_complete+nr_extra-1;
    int stop_s=clock();
    e_time += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    #if LOGGING
    printf("After Extend(%u):", n); check(w, 0, w_end, d);
    #endif
}

void Delete (uint32_t* w, uint32_t p, uint32_t from, uint32_t to, char *d, uint32_t* f_deleted, uint32_t fdbound, uint32_t &nfd, uint32_t N) {
    // Deletes multiples p*w[i] from W for from <= i <= to
    int start_s = clock();
    // for (i=from; i <= to; i++) delete(p*w[i],d);
    alignas(32) uint32_t buf[8]; alignas(32) uint32_t buf2[8]; // 16 is better if processor supports
    alignas(32) uint32_t ps[8]; for (uint32_t k=0; k < 8; k++) ps[k] = p;
    uint32_t i;
    __m256i Vec1 = _mm256_set1_epi32(1);
    __m256i Vec3 = _mm256_set1_epi32(3);
    __m256i Vec7 = _mm256_set1_epi32(7);
    __m256i Vecp = _mm256_loadu_si256((__m256i*)ps);
    for (i=from; i+7 <= to; i+=8) {
        //memcpy(buf, w+i, 8*sizeof(uint32_t));
        __m256i Vecb = _mm256_loadu_si256((__m256i*)(w+i)); // Vecb = buf
        //for (uint32_t k=0; k < 8; k++) buf[k] *= p;
        Vecb = _mm256_mullo_epi32(Vecb, Vecp);
        //for (uint32_t k=0; k < 8; k++) buf2[k] = buf[k] & 7; // % 8
        __m256i Vecb2 = (__m256i)_mm256_and_ps((__m256)Vecb, (__m256)Vec7); // Vecb2 = buf2
        //for (uint32_t k=0; k < 8; k++) buf[k] /= 8;
        Vecb = _mm256_srav_epi32(Vecb, Vec3);
        //for (uint32_t k=0; k < 8; k++) d[buf[k]] |= 1 << buf2[k];
        Vecb2 = _mm256_sllv_epi32(Vec1, Vecb2); // Vecb2 = 1 << buf2
        _mm256_storeu_si256((__m256i *) buf, Vecb);
        _mm256_storeu_si256((__m256i *) buf2, Vecb2);
        for (uint32_t k=0; k < 8; k++) d[buf[k]] |= buf2[k];
    }
    for (i=i; i <= to; i++) delete(p*w[i],d); // do < 8 remaining
    // remember deleted factors for Compress
    i = from;
    uint32_t j = 0;
    while (1) {
        uint32_t c = p*w[i++];
        if (c > fdbound) break;
        f_deleted[j++] = c;
    }
    nfd = j;
    int stop_s=clock();
    d_time += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    #if LOGGING
    printf("After Delete(%u, %u->%u): (partial)", p, from, to); check(w, 0, to, d);
    #endif
}

void Compress(uint32_t* w, char* d, uint32_t from, uint32_t &to, uint32_t w_end, uint32_t* f_deleted, uint32_t nfd, bool isFullWheel) {
    // Removes deleted values in w[from..to], updates to, and pads with evens if to < w_end
    int start_s = clock();
    uint32_t padding = w[to]+1; // even
    uint32_t to_init = to, dest = from, source = from;
    uint32_t i = nfd-1;
    while (f_deleted[i] > w[to]) i--;
    f_deleted[i+1] = w[to+1]; nfd = i+2; // sentinel
    uint32_t ifd = 0;
    i = from;
    if (isFullWheel) {
        while (ifd < nfd) {
            //while (w[i] < f_deleted[ifd]) i++;
            i = find(w, f_deleted[ifd], i, to+2); // (not safe to use doubling search)
            memcpy(w+dest, w+source, (i-source)*sizeof(uint32_t));
            dest += (i-source);
            source = ++i;
            ifd++;
        }
    } else {
        while (ifd < nfd) {
            //while (w[i] < f_deleted[ifd]) i++;
            uint32_t delta = 1;
            while (w[i+delta] <= f_deleted[ifd]) delta *= 2; // safe?!
            i = find(w, f_deleted[ifd], i, i+delta);
            memcpy(w+dest, w+source, (i-source)*sizeof(uint32_t));
            dest += (i-source);
            source = ++i;
            ifd++;
        }
    }
    if (to_init < w_end) { // pad with evens to support search in Compress
        uint32_t* at = w+dest;
        for (uint32_t* end = w+to_init; at <= end; at++) *at = padding; // (should vectorize)
    }
    to = dest-1; // return value
    int stop_s=clock();
    c_time += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    #if LOGGING
    printf("After Compress(%u,%u): (partial)", from, to_init); check(w, from, to, d);
    #endif
}

uint32_t Sift(uint32_t N, bool printPrimes) {
    // returns the number of primes up to N, printing them if printPrimes
    uint32_t *w = new uint32_t[2];
    char* d = (char*) calloc((N+1)/8+1, 1); // initialized to zero=false
    /* representation invariant (for the main loop):
       if length < N (so W is a complete wheel), w[0..w_end] is the ordered set W;
       otherwise, w[0..w_end] is ordered, and omitting evens and values w with d[w] true, is the set W \union Pr - {p_i|i<=nrPrimes},
       and no even or deleted values occur in w before a value between 1 and N/p' (where p' is the prime after p) */
    uint32_t p, p_index, p2, p_next, w_end_prev, low, high, delta;
    // W,k,length = {1},1,2:
    uint32_t w_end = 0; w[0] = 1;
    w[1] = 2; // sentinel
    uint32_t length = 2;
    uint32_t nr_appended = 0, nr_deleted = 0;
    uint32_t* f_deleted; uint32_t nfd;
    // Pr = {2}:
    if (printPrimes) printf("%u\n", 2); /*nrPrimes = 1;*/
    p = 3; p_index = 1; p2 = p*p;
    uint32_t c_from = 1, imaxf = 0;
    /* invariant: p = p_(k+1) and W = W_k inter {1,...,N} and length = min(P_k,N) and Pr = the first k primes
       (where p_i denotes the i'th prime, W_i denotes the i'th wheel, P_i denotes the product of the first i primes) */
    while (p2 <= N) { // p^2 <= N
        if (length < N) {
            w_end_prev = w_end;
            Extend (w, w_end, length, std::min((uint64_t)p*length,(uint64_t)N), d);
            nr_appended += w_end-w_end_prev;
        }
        if (length < N) {
            f_deleted = new uint32_t[w_end+1];
            if (printPrimes) printf("%u\n", p);/* nrPrimes++;*/
            Delete(w, p, 0, w_end_prev, d, f_deleted, length, nfd, N); // starts from w[0] to delete p*w[0]=p*1=p
            nr_deleted += w_end_prev; // don't count prime deleted
            Compress(w, d, c_from, w_end, w_end, f_deleted, nfd, true); // c_from=1 since p=w[1] is deleted
        } else {
            p_next = w[p_index+1];
            if (imaxf == 0) { imaxf = find(w, N/p, 1, w_end); nfd = find(w, N/p2, 1, imaxf); } // first time
            f_deleted = new uint32_t[nfd];
            Delete(w, p, p_index, imaxf, d, f_deleted, N/p_next, nfd, N); // starts from w[p_index]=p so p^2 is first deletion
            nr_deleted += imaxf-p_index+1;
            while (w[imaxf] > N/p_next) imaxf--; // doubling then binary search could be used here
            if (p2 <= N/p_next) {
                while (w[c_from] < p2) c_from++; // doubling then binary search could be used here
                Compress(w, d, c_from, imaxf, w_end, f_deleted, nfd, false);
            }
            p_index++;
        }
        // p = next(W, 1):
        p = w[p_index]; p2 = p*p;
        // k++
        delete[] f_deleted;
    }
    if (length < N) {
        w_end_prev = w_end;
        Extend (w, w_end, length, N, d);
        nr_appended += w_end-w_end_prev;
    }
    // print remaining primes:
    if (printPrimes) {
        int start_s = clock();
        for (uint32_t i=1; i <= w_end; i++) {
            if (w[i]%2 == 0 || deleted(w[i],d)) continue;
            printf("%u\n", w[i]);
        }
        int stop_s=clock();
        p_time += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    }
    delete[] w; free(d);
    return (1+nr_appended-nr_deleted);
}

int main (int argc, char *argw[]) {
    bool error = false, printPrimes = false;
    uint64_t max = 4000000000, N;
    uint32_t PiN;
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
    PiN = Sift(N, printPrimes);
    int stop_s=clock();
    float duration = (stop_s-start_s)*1E3/double(CLOCKS_PER_SEC);
    printf("%u primes up to %lu found in %.3f ms\n", PiN, N, duration);
    printf("Time for: Extend=%.2f ms, Delete=%.2f ms, Compress=%.2f ms, Other=%.2f ms",
      e_time*1E3, d_time*1E3, c_time*1E3, duration-(e_time+d_time+c_time+p_time)*1E3);
    if (printPrimes) printf(", Printing=%.2f ms", p_time*1E3);
    printf("\n");
}