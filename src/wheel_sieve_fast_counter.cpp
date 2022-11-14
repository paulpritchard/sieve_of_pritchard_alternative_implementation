/* Sieve of Pritchard in C++, as described at https://en.wikipedia.org/wiki/Sieve_of_Pritchard */
/* arguments: N [-p] */
/*    N: counts the primes up to N */
/*   -p: (optional) print the primes up to N */
/* careful high-performance implementation using a simple array of integers and a bit array (using bytes for speed), */
/* including method for very fast counting ("The Count" remix) */
/* With -p: 2 <= N <= 1000000000 */
/* (like the standard Sieve of Eratosthenes, this algorithm is not suitable for very large N due to memory requirements) */
/* Without -p: 2 <= N <= 100000000000 */
/* (unlike the standard Sieve of Eratosthenes!) */

/* compile with -DLOGGING to print array after each main step */

#include <cstring>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctime>
#include <getopt.h>

double e_time = 0.0; /* Extend */
double d_time = 0.0; /* Delete */
double c_time = 0.0; /* Compress */
double p_time = 0.0; /* Printing */

void check(uint32_t w[], uint32_t from, uint32_t to, bool d[]) { /* for logging array contents */
    for (uint32_t i=from; i <= to; i++) printf(" w[%d]=%d%s", i, w[i], (d[w[i]] ? "*" : ""));
    printf("\n");
}

uint32_t find (uint32_t w[], uint32_t bound, uint32_t low, uint32_t high) {
    /* returns max{i: w[i] <= bound} given w[low] <= bound < w[high]*/
    uint32_t i, j, mid;
    i = low; j = high;
    /* invariant: w[i] <= bound < w[j] */
    while (i+1 < j) if (w[mid = (i+j)/2] <= bound) i = mid; else j = mid;
    return(i);
}

void Extend (uint32_t* &w, uint32_t &w_end, uint64_t &length, uint64_t n, bool* &d) {
    /* Rolls full wheel W up to n, and sets length=n */
    int start_s = clock();
    /* get bigger arrays: */
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
    /* */
    i = 0;
    for (j=w_end+1; j <= w_end_new; j++) {
        x = length + w[i++];
        /* Append x to the ordered set W: */
        w[j] = x; d[x] = false;
    }
    length = n; w_end = w_end_new;
    int stop_s=clock();
    e_time += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    #if LOGGING
    printf("After Extend(%lu):", n); check(w, 0, w_end, d);
    #endif
}

void Delete (uint32_t w[], uint32_t p, uint32_t from, uint32_t to, bool d[]) {
    /* Deletes multiples p*w[i] from W for from <= i <= to */
    int start_s = clock();
    for (uint32_t i=from; i <= to; i++) d[p*w[i]] = true;
    int stop_s=clock();
    d_time += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    #if LOGGING
    printf("After Delete(%d, %d->%d):", p, from, to); check(w, 0, to, d);
    #endif
}

void Compress(uint32_t w[], bool d[], uint32_t from, uint32_t &to, uint32_t w_end) {
    /* Removes deleted values in w[from..to], updates to, and pads with zeros if to < w_end */
    int start_s = clock();
    uint32_t to_init = to;
    uint32_t j = from;
    for (uint32_t i=from; i <= to_init; i++) if (!d[w[i]]) w[j++] = w[i];
    if (to_init < w_end) for (uint32_t k=j; k <= to_init; k++) w[k] = 0; /* pad with zeros */
    to = j-1; /* return value */
    int stop_s=clock();
    c_time += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    #if LOGGING
    printf("After Compress(%d,%d):", from, to); check(w, 0, w_end, d);
    #endif
}

void Sift(uint32_t N, uint32_t &nr) {
    /* prints the nr primes up to N */
    uint32_t *w = new uint32_t[1];
    bool *d = new bool[2];
    /* representation invariant (for the main loop): */
    /* if length < N (so W is a complete wheel), w[0..w_end] is the ordered set W; */
    /* otherwise, w[0..w_end], omitting zeros and values w with d[w] true, is the ordered set W \union Pr - {p_i|i<=nrPrimes}, */
    /* and no zero or deleted values occur in w before a value between 1 and N/p' (where p' is the prime after p) */
    uint32_t p, p_index, p2, p_next, w_end_prev;
    uint64_t limit;
    /* W,k,length = {1},1,2: */
    uint32_t w_end = 0; w[0] = 1; d[1] = false; uint64_t length = 2;
    uint32_t nr_appended = 0, nr_deleted = 0;
    /* Pr = {2}: */
    /* nrPrimes = 1; */ printf("%d\n", 2);
    p = 3; p_index = 1; p2 = p*p;
    uint32_t c_from = 1;
    uint32_t iNoverp = 0; /* max{i: 0 < w[i] <= N/p (before deletion) }*/
    /* invariant: p = p_(k+1) and W = W_k inter {1,...,N} and length = min(P_k,N) and Pr = the first k primes */
    /* (where p_i denotes the i'th prime, W_i denotes the i'th wheel, P_i denotes the product of the first i primes) */
    while (p2 <= N) { /* p^2 <= N */
        if (length < N) {
            w_end_prev = w_end;
            /* Extend W,length to minimum of p*length,N: */
            limit = p*length; if (limit > N) limit = N;
            Extend (w, w_end, length, limit, d);
            nr_appended += w_end-w_end_prev;
        }
        if (length < N) {
            /* nrPrimes++; */ printf("%d\n", p);
            Delete(w, p, 0, w_end_prev, d); /* starts from w[0] to delete p*w[0]=p*1=p */
            nr_deleted += w_end_prev; /* don't count prime deleted */
            Compress(w, d, c_from, w_end, w_end); /* c_from=1 since p=w[1] is deleted */
        } else {
            if (iNoverp == 0) iNoverp = find(w, N/p, 1, w_end); /* first time */
            Delete(w, p, p_index, iNoverp, d); /* starts from w[p_index]=p so p^2 is first deletion */
            nr_deleted += iNoverp-p_index+1;
            p_next = w[++p_index];
            if (p_next == 0) break; /* next p is after zeroed section so is too big */
            while (w[iNoverp] > N/p_next) iNoverp--; /* doubling then binary search could be used here */
            if (p2 <= N/p_next) {
                while (w[c_from] < p2) c_from++; /* doubling then binary search could be used here */
                Compress(w, d, c_from, iNoverp, w_end);
            }
        }
        /* p = next(W, 1): */
        p = w[p_index]; p2 = p*p;
        /* k++ */
    }
    if (length < N) {
        w_end_prev = w_end;
        Extend (w, w_end, length, N, d);
        nr_appended += w_end-w_end_prev;
    }
    /* print remaining primes: */
    int start_s = clock();
    for (uint32_t i=1; i <= w_end; i++) {
        if (w[i] == 0 || d[w[i]]) continue;
        printf("%d\n", w[i]);
    }
    int stop_s=clock();
    p_time += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    delete[] w; delete[] d;
    nr = 1+nr_appended-nr_deleted;
    printf("%d primes up to %d = 1 (for 2) + %d (appended) - %d (deleted).\n", nr, N, nr_appended, nr_deleted);
}

uint64_t Count(uint64_t N) {
    /* returns the number of primes up to N */
    if (N < 9) return((N+1)/2);
    uint32_t *w = new uint32_t[1];
    bool *d = new bool[3]; /* since length=2 */
    /* representation invariant (for the main loop): */
    /* if length < N (so W is a complete wheel), w[0..w_end] is the ordered set W; */
    /* otherwise, w[0..w_end], omitting zeros and values w with d[w] true, is the ordered set (W \union Pr - {p_i|i<=nrPrimes}) \inter {1,...,N/p}, */
    /* and no zero or deleted values occur in w before a value between 1 and N/p' (where p' is the prime after p) */
    uint32_t p, p_index, p_next, w_end_prev, rem;
    uint64_t p2, limit;
    uint32_t iNoverp = 0; /* max {i: w[i]<=N/p} */
    uint32_t iNoverp2 = 0; /* max {i: w[i]<=N/p^2} */
    /* W,k,length = {1},1,2: */
    uint32_t w_end = 0; w[0] = 1; d[1] = false; uint64_t length = 2;
    uint64_t nr_appended = 0, nr_deleted = 0; /* notional */
    uint64_t nr_compressed = 0;
    /* Pr = {2}: */
    /* nrPrimes = 1; */
    p = 3; p_index = 1; p2 = p*p;
    uint32_t c_from = 1;
    /* invariant: p = p_(k+1) and W = W_k inter {1,...,N} and length = min(P_k,N) and Pr = the first k primes */
    /* (where p_i denotes the i'th prime, W_i denotes the i'th wheel, P_i denotes the product of the first i primes) */
    while (p*length <= N/p) { /* ensure last full wheel is not too big */
        w_end_prev = w_end;
        Extend(w, w_end, length, p*length, d);
        nr_appended += w_end-w_end_prev;
        Delete(w, p, 0, w_end_prev, d); /* starts from w[0] to delete p*w[0]=p*1=p */
        nr_deleted += w_end_prev; /* don't count prime deleted */
        nr_compressed += w_end-c_from+1;
        Compress(w, d, c_from, w_end, w_end); /* c_from=1 since p=w[1] is deleted */
        /* p = next(W, 1): */
        p = w[p_index];
        p2 = (uint64_t)p*(uint64_t)p;
        /* k++ */
    }
    /* extend to N notionally */
    nr_appended += (N/length-1)*(w_end+1);
    rem = N % length;
    if (rem > 0) nr_appended += find(w, rem, 0, w_end+1)+1;
    limit = N/p;
    /* use w up to N/p */
    if (limit > length) {
        Extend(w, w_end, length, limit, d);
    } else { /* occurs e.g. when N=25,150 */
        while (w[w_end] > N/p) w_end--; /* doubling then binary search could be used here */
    }
    nr_deleted += w_end;
    length = N; /* notional */
    iNoverp2 = find(w, N/p2, 0, w_end);
    if (p2 <= N/p) {
        Delete(w, p, 1, iNoverp2, d);
        while (w[c_from] < p2) c_from++; /* doubling then binary search could be used here */
        nr_compressed += w_end-c_from+1;
        Compress(w, d, c_from, w_end, w_end);
    }
    iNoverp = w_end; /* notional */
    /* p = next(W, 1): */
    if (p_index == w_end) return(1+nr_appended-nr_deleted);
    p = w[++p_index];
    p2 = (uint64_t)p*(uint64_t)p;
    /* k++ */
    while (p2 <= N/p) {
        while (w[iNoverp] > N/p) iNoverp--; /* doubling then binary search could be used here */
        nr_deleted += iNoverp-p_index+1;
        while (w[iNoverp2] > N/p2) iNoverp2--; /* doubling then binary search could be used here */
        Delete(w, p, p_index, iNoverp2, d); /* don't delete p */
        while (w[c_from] < p2) c_from++; /* doubling then binary search could be used here */
        nr_compressed += iNoverp-c_from+1;
        Compress(w, d, c_from, iNoverp, w_end);
        /* p = next(W, 1): */
        p = w[++p_index];
        p2 = (uint64_t)p*(uint64_t)p;
        /* k++ */
    }
    while (p2 <= N) {
        while (w[iNoverp] > N/p) iNoverp--; /* doubling then binary search could be used here */
        nr_deleted += iNoverp-p_index+1;
        /* p = next(W, 1): */
        if (p_index > w_end) break; /* needed? */
        p = w[++p_index];
        if (p == 0) break; /* next p is after zeroed section so is too big */
        p2 = (uint64_t)p*(uint64_t)p;
        /* k++ */
    }
    delete[] w; delete[] d;
    printf("%lu compressed\n", nr_compressed);
    printf("%llu primes up to %llu = 1 (for 2) + %llu (appended) - %llu (deleted).\n", (1+nr_appended-nr_deleted), N, nr_appended, nr_deleted);
    return(1+nr_appended-nr_deleted);
}

int main (int argc, char *argw[]) {
    bool error = false; bool printPrimes = false; uint64_t maxPrint = 1000000000; uint64_t maxCount = 100000000000;
    uint32_t np;
    uint64_t N, nptoN;
    uint64_t max = maxCount;
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
        printf("call with: %s N where 2 <= N <= %llu to count, or %s N -p where 2 <= N <= %llu to print\n", argw[0], maxCount, argw[0], maxPrint);
        exit(1);
    }
    int start_s = clock();
    if (printPrimes) {
        Sift(N, np);
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
