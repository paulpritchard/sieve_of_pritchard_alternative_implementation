/* Sieve of Pritchard in C++, as described at https://en.wikipedia.org/wiki/Sieve_of_Pritchard
   arguments: N [-p]
      N: counts the primes up to N
     -p: (optional) print the primes up to N
   careful high-performance implementation using a simple array of integers and a bit array (using bytes for speed),
   including method for very fast counting ("The Count" remix)
   With -p: 2 <= N <= 1000000000
   (like the standard Sieve of Eratosthenes, this algorithm is not suitable for very large N due to memory requirements)
   Without -p: 2 <= N <= 100000000000
   (unlike the standard Sieve of Eratosthenes!) */

// compile with -DLOGGING to print array after each main step

#include <cstring>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctime>
#include <getopt.h>
#include <cmath>

double e_time = 0.0; // Extend
double d_time = 0.0; // Delete
double c_time = 0.0; // Compress
double p_time = 0.0; // Printing

void wCheck (uint32_t w[], uint32_t from, uint32_t to, bool d[]) { // for logging array contents
    for (uint32_t i=from; i <= to; i++) printf(" w[%d]=%d%s", i, w[i], (d[w[i]] ? "*" : ""));
    printf("\n");
}

void countsCheck (uint32_t nprime, uint32_t prime[], uint32_t nrWheelPrimes, uint64_t nrNoverpi[]) {
    for (uint32_t i=nrWheelPrimes; i < nprime; i++) printf(" nr[%d] for %d = %lu", i, prime[i], nrNoverpi[i]);
    printf("\n");
}

uint32_t find (uint32_t w[], uint32_t bound, uint32_t low, uint32_t high) {
    // returns max{i: w[i] <= bound} given w[low] <= bound < w[high]
    uint32_t i, j, mid;
    i = low; j = high;
    // invariant: w[i] <= bound < w[j]
    while (i+1 < j) if (w[mid = (i+j)/2] <= bound) i = mid; else j = mid;
    return(i);
}

void Extend (uint32_t* &w, uint32_t &w_end, uint64_t &length, uint64_t n, bool* &d) {
    // Rolls full wheel W up to n, and sets length=n
    int start_s = clock();
    // get bigger arrays:
    uint32_t w_end_new = (n/length)*(w_end+1)-1;
    uint32_t rem = n % length;
    uint32_t i, j, x;
    if (rem > 0) w_end_new += find(w, rem, 0, w_end+1)+1;
    uint32_t *w_next = new uint32_t[w_end_new+1];
    memcpy(w_next, w, (w_end+1)*sizeof(uint32_t));
    delete[] w; w = w_next;
    bool *d_next = new bool[n+1];
    memcpy(d_next, d, (length+1)*sizeof(bool));
    delete[] d; d = d_next;

    i = 0;
    for (j=w_end+1; j <= w_end_new; j++) {
        x = length + w[i++];
        // Append x to the ordered set W:
        w[j] = x; d[x] = false;
    }
    length = n; w_end = w_end_new;
    int stop_s=clock();
    e_time += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    #if LOGGING
    printf("After Extend(%lu):", n); wCheck(w, 0, w_end, d);
    #endif
}

void Delete (uint32_t w[], uint32_t p, uint32_t from, uint32_t to, bool d[]) {
    // Deletes multiples p*w[i] from W for from <= i <= to
    int start_s = clock();
    for (uint32_t i=from; i <= to; i++) d[p*w[i]] = true;
    int stop_s=clock();
    d_time += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    #if LOGGING
    printf("After Delete(%d, %d->%d):", p, from, to); wCheck(w, 0, to, d);
    #endif
}

void Compress(uint32_t w[], bool d[], uint32_t from, uint32_t &to, uint32_t w_end) {
    // Removes deleted values in w[from..to], updates to, and pads with zeros if to < w_end
    int start_s = clock();
    uint32_t to_init = to;
    uint32_t j = from;
    for (uint32_t i=from; i <= to_init; i++) if (!d[w[i]]) w[j++] = w[i];
    if (to_init < w_end) for (uint32_t k=j; k <= to_init; k++) w[k] = 0; // pad with zeros
    to = j-1; // return value
    int stop_s=clock();
    c_time += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    #if LOGGING
    printf("After Compress(%d,%d):", from, to); wCheck(w, 0, w_end, d);
    #endif
}

void Sift(uint32_t N, uint32_t* &prime, uint32_t &nr) {
    // stores the nr primes up to N in order in prime[] (with an extra courtesy element for use as sentinel!)
    uint32_t *w = new uint32_t[1];
    bool *d = new bool[2];
    /* representation invariant (for the main loop):
       if length < N (so W is a complete wheel), w[0..w_end] is the ordered set W;
       otherwise, w[0..w_end], omitting zeros and values w with d[w] true, is the ordered set W \union Pr - {p_i|i<=nrPrimes},
       and no zero or deleted values occur in w before a value between 1 and N/p' (where p' is the prime after p) */
    uint32_t p, p_index, p2, p_next, w_end_prev, limit;
    // W,k,length = {1},1,2:
    uint32_t w_end = 0; w[0] = 1; d[1] = false; uint64_t length = 2;
    uint32_t nr_appended = 0, nr_deleted = 0;
    // Pr = {2}:
    uint32_t *prime_deleted = new uint32_t[10];
    prime_deleted[0] = 2; uint32_t npdeleted = 1;
    p = 3; p_index = 1; p2 = p*p;
    uint32_t c_from = 1;
    uint32_t iNoverp = 0; // max{i: 0 < w[i] <= N/p (before deletion) }
    /* invariant: p = p_(k+1) and W = W_k inter {1,...,N} and length = min(P_k,N) and Pr = the first k primes
       (where p_i denotes the i'th prime, W_i denotes the i'th wheel, P_i denotes the product of the first i primes) */
    while (p2 <= N) { // p^2 <= N
        if (length < N) {
            w_end_prev = w_end;
            // Extend W,length to minimum of p*length,N:
            limit = p*length; if (limit > N) limit = N;
            Extend (w, w_end, length, limit, d);
            nr_appended += w_end-w_end_prev;
        }
        if (length < N) {
            prime_deleted[npdeleted++] = p;
            Delete(w, p, 0, w_end_prev, d); // starts from w[0] to delete p*w[0]=p*1=p
            nr_deleted += w_end_prev; // don't count prime deleted
            Compress(w, d, c_from, w_end, w_end); // c_from=1 since p=w[1] is deleted
        } else {
            if (iNoverp == 0) iNoverp = find(w, N/p, 1, w_end); // first time
            Delete(w, p, p_index, iNoverp, d); // starts from w[p_index]=p so p^2 is first deletion
            nr_deleted += iNoverp-p_index+1;
            p_next = w[++p_index];
            if (p_next == 0) break; // next p is after zeroed section so is too big
            while (w[iNoverp] > N/p_next) iNoverp--; // doubling then binary search could be used here
            if (p2 <= N/p_next) {
                while (w[c_from] < p2) c_from++; // doubling then binary search could be used here
                Compress(w, d, c_from, iNoverp, w_end);
            }
        }
        // p = next(W, 1):
        p = w[p_index]; p2 = p*p;
        // k++
    }
    if (length < N) {
        w_end_prev = w_end;
        Extend (w, w_end, length, N, d);
        nr_appended += w_end-w_end_prev;
    }
    // get all primes:
    int start_s = clock();
    uint32_t *pr = new uint32_t[1+nr_appended-nr_deleted+1]; // 1 extra
    uint32_t k;
    for (k=0; k < npdeleted; k++) pr[k] = prime_deleted[k];
    for (uint32_t i=1; i <= w_end; i++) {
        if (w[i] == 0 || d[w[i]]) continue;
        pr[k++] = w[i];
    }
    prime = pr; nr = k;
    int stop_s=clock();
    p_time += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    delete[] w; delete[] d; delete[] prime_deleted;
}

void getCounts (uint64_t N, uint32_t w[], uint32_t w_end, uint64_t length, uint32_t nrWheelPrimes, uint32_t nprime, uint32_t prime[], uint64_t nrNoverpi[]) {
    uint64_t limit, nr, rem;
    for (uint32_t i=nrWheelPrimes; i < nprime; i++) { // skip primes in wheel
        limit = N/prime[i];
        nr = (limit/length)*(w_end+1);
        rem = limit % length;
        if (rem > 0) nr += find(w, rem, 0, w_end+1)+1;
        nrNoverpi[i] = nr-i+nrWheelPrimes-1;
    }
}

void updateCounts (uint64_t N, uint32_t low, uint32_t w[], uint32_t w_end, uint32_t p, uint32_t p_index, uint32_t prime[], uint64_t nrNoverpi[]) {
    uint32_t imaxf;
    uint32_t i = low, high = w_end+1;
    while ((uint64_t)p*(uint64_t)p <= N/prime[i]) {
        imaxf = find(w, N/(p*prime[i]), 0, high);
        nrNoverpi[i] -= (imaxf-p_index+1);
        high = imaxf+1;
        i++;
    }
}

uint64_t Count(uint64_t N) {
    // returns the number of primes up to N
    uint32_t *w = new uint32_t[1];
    bool *d = new bool[3]; // since length=2
    /* representation invariant (for the main loop):
       if length < N (so W is a complete wheel), w[0..w_end] is the ordered set W;
       otherwise, w[0..w_end], omitting zeros and values w with d[w] true, is the ordered set (W \union Pr - {p_i|i<=nrPrimes})
       when restricted to the interval \inter {1,...,N/p},
       and no zero or deleted values occur in w before a value between 1 and N/p' (where p' is the prime after p),
       and prime[i] = p_(i+1) for 0 <= i < number of primes up to sqrt(N),
       and p = prime[ip], and nrNoverpi[i] = |W \inter {prime[i],...,N/prime[i]}| */
    uint32_t p, p_index, w_end_prev, rem;
    uint64_t p2, limit;
    uint32_t iNoverp = 0; // max {i: w[i]<=N/p}
    uint32_t iNoverp2 = 0; // max {i: w[i]<=N/p^2}
    uint32_t iNoverp3;
    // W,k,length = {1},1,2:
    uint32_t w_end = 0; w[0] = 1; d[1] = false; uint64_t length = 2;
    uint64_t nr_appended = 0, nr_deleted = 0; // notional
    // Pr = {2}:
    uint32_t nrWheelPrimes = 1;
    p = 3; p_index = 1; p2 = p*p;
    uint32_t c_from = 1;
    uint32_t sqrtN = std::sqrt(N);
    uint32_t npsqrtN;
    uint32_t *prime;
    uint64_t *nrNoverpi;
    if (N < 9) return((N+1)/2);
    Sift(sqrtN, prime, npsqrtN); prime[npsqrtN] = N+1; // sentinel
    nrNoverpi = new uint64_t[npsqrtN];
    uint32_t ip = 1;
    /* invariant for all main loops (notional): p = p_(k+1) and W = W_k inter {1,...,N} and length = min(P_k,N) and Pr = the first k primes
       (where p_i denotes the i'th prime, W_i denotes the i'th wheel, P_i denotes the product of the first i primes) */
    while (p*length <= N/p2) { // ensure last full wheel is not too big
        w_end_prev = w_end;
        Extend(w, w_end, length, p*length, d);
        nr_appended += w_end-w_end_prev;
        Delete(w, p, 0, w_end_prev, d); // starts from w[0] to delete p*w[0]=p*1=p
        nr_deleted += w_end_prev; // don't count prime deleted
        Compress(w, d, c_from, w_end, w_end); // c_from=1 since p=w[1] is deleted
        nrWheelPrimes++;
        // p = next(W, 1):
        p = prime[++ip];
        p2 = (uint64_t)p*(uint64_t)p;
        // k++
    }
    // extend to N notionally
    getCounts(N, w, w_end, length, nrWheelPrimes, npsqrtN, prime, nrNoverpi);
    nr_appended += (N/length-1)*(w_end+1);
    rem = N % length;
    if (rem > 0) nr_appended += find(w, rem, 0, w_end+1)+1;
    nr_deleted += nrNoverpi[ip];
    // use w up to N/p^2
    limit = N/p2;
    if (limit > length) {
        Extend(w, w_end, length, limit, d); 
    } else {
        while (w[w_end] > N/p2) w_end--; // doubling then binary search could be used here
    }
    length = N; // notional
    if (p2 <= N/p) updateCounts(N, ip+1, w, w_end, p, p_index, prime, nrNoverpi);
    if (p2 <= N/p2) {
        iNoverp3 = find(w, N/(p2*p), 0, w_end);
        Delete(w, p, 1, iNoverp3, d); // don't delete p
        while (w[c_from] < p2) c_from++; // doubling then binary search could be used here
        Compress(w, d, c_from, w_end, w_end);
    }
    p_index++;
    p = prime[++ip];
    p2 = (uint64_t)p*(uint64_t)p;
    while (p2 <= N/p2) { // p <= N^(1/4)
        // k++
        nr_deleted += nrNoverpi[ip];
        updateCounts(N, ip+1, w, w_end, p, p_index, prime, nrNoverpi);
        while (w[iNoverp3] > N/(p2*p)) iNoverp3--; // doubling then binary search could be used here
        Delete(w, p, p_index, iNoverp3, d); // don't delete p
        while (w[w_end] > N/p2) w_end--; // doubling then binary search could be used here
        while (w[c_from] < p2) c_from++; // doubling then binary search could be used here
        Compress(w, d, c_from, w_end, w_end);
        p_index++;
        p = prime[++ip];
        p2 = (uint64_t)p*(uint64_t)p;
    }
    while (p2 <= N/p) { // p <= N^(1/3)
        // k++
        nr_deleted += nrNoverpi[ip];
        updateCounts(N, ip+1, w, w_end, p, p_index, prime, nrNoverpi);
        p_index++;
        p = prime[++ip];
        p2 = (uint64_t)p*(uint64_t)p;
    }
    // finish sifting for p > N^(1/3)
    for (ip=ip; ip < npsqrtN; ip++) {
        // p = next(W, 1):
        p = prime[ip];
        // k++
        nr_deleted += nrNoverpi[ip];
    }
    delete[] w; delete[] d; delete[] prime; delete[] nrNoverpi;
    // printf("%lu primes up to %lu = 1 (for 2) + %lu (appended) - %lu (deleted).\n", (1+nr_appended-nr_deleted), N, nr_appended, nr_deleted);
    return(1+nr_appended-nr_deleted);
}

int main (int argc, char *argw[]) {
    bool error = false; bool printPrimes = false; uint64_t maxPrint = 1000000000; uint64_t maxCount = 1000000000000; 
    uint64_t N, nptoN;
    uint64_t max = maxCount;
    uint32_t *prime, np;
    if (argc == 3) {
        if (strcmp(argw[2], "-p") == 0) {
            printPrimes = true; max = maxPrint;
            argc--;
        } else {
            error = true;
        }
    }
    if (argc == 2) {
        N = atol(argw[1]);
        if (N < 2 || N > max) error = true;
    } else {
        error = true;
    }
    if (error) {
        printf("call with: %s N where 2 <= N <= %lu to count, or %s N -p where 2 <= N <= %lu to print\n", argw[0], maxCount, argw[0], maxPrint);
        exit(1);
    }
    int start_s = clock();
    if (printPrimes) {
        Sift(N, prime, np);
        for (uint32_t i=0; i < np; i++) printf("%d\n", prime[i]);
        nptoN = np;
    } else {
        nptoN = Count(N);
    }
    printf("%lu primes up to %lu\n", nptoN, N);
    int stop_s=clock();
    float duration = (stop_s-start_s)*1E3/double(CLOCKS_PER_SEC);
    printf("Total time %.2f ms\n", duration);
    printf("Time for: Extend=%.2f ms, Delete=%.2f ms, Compress=%.2f ms, Other=%.2f ms",
      e_time*1E3, d_time*1E3, c_time*1E3, duration-(e_time+d_time+c_time+p_time)*1E3);
    if (printPrimes) printf(", Printing=%.2f ms", p_time*1E3);
    printf("\n");
}