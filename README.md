# An alternative implementation of the dynamic wheel sieve of Pritchard

This repository concerns the dynamic wheel sieve of Pritchard as described on its [Wikipedia page](https://en.wikipedia.org/wiki/Sieve_of_Pritchard).
It is an algorithm for finding the primes up to a given limit N that takes sublinear time O(N/log log N).
The repository contains several variations of an alternative implementation of the algorithm.

The original implementation is described in the paper
[Paul Pritchard, "A Sublinear Additive Sieve for Finding Prime Numbers", *Communications of the ACM*, vol. 24, no. 1, pp. 18–23](https://dl.acm.org/doi/10.1145/358527.358540).
A detailed code animation with a limit of N=150 is given in [this video](https://www.youtube.com/watch?v=GxgGMwLfTjE).

Since the dynamic wheel sieve was presented as a contribution to computational complexity (and as a beautiful algorithm),
the original implementation was designed to demonstrate the sublinear running time as simply and clearly as possible.

We provide here an alternative and much more efficient implementation.
it is both faster (by a constant factor) and requires less space than the original implementation.

Note that the pure bitmap implementation in the repository
[Sieve_of_Pritchard_Bitmap_Implementation](https://github.com/paulpritchard/Sieve_of_Pritchard_Bitmap_Implementation) is even faster again.
If it's pure speed you are after, proceed directly there. Do not pass GO.
However, the implementations is this repository are of independent interest because they naturally culminate in
algorithms for counting the primes much faster than they can be generated.

## implementations

The code for various implementations resides in the src directory.
These implementations are discussed below in order of increasing sophistication and better performance.

### wheel_sieve_simple.cpp

This is a self-contained C++ program.
It gives a simple but fast implementation of the dynamic wheel sieve,
to best illustrate the ideas behind the alternative implementations.
There are no changes to the abstract algorithm.

The sieve maintains an ordered set W which is the k'th wheel
(i.e. the ordered set of natural numbers up to the product of the first k primes that are not divisible by any of them)
restricted to the integers up to N.
To achieve a fast implementation, very simple data structures are used, which are processed with tight loops.
(This is why the sieve of Eratosthenes is so effective.
It is the Creedence Clearwater Revival of prime generators!)

The starting point is to represent W as simply as possible, with an array w[] containing the members in order;
i.e. w[i] is the i'th member (indexing from 0).
When the current wheel is extended by rolling it, the code simply iterates through the array w,
adding the length of the wheel to each member w[i] and appending the result.
The other step is to delete the composites formed by multiplying the values in the current wheel by the current prime p.
However, this presents problems, firstly because each deleted multiple cannot be found in w in O(1) time.
Accordingly, a bit array d[] (for "deleted") is introduced such that d[x] is initialized to false when a value x is appended to W,
and is set to true should x be deleted as a multiple x = p*w[i] of p.
Deletions are now fast, but the array is left containing deleted elements.

So if the new W will be extended in the next iteration, because its length < N, then the array w is compressed by eliminating the deleted values.
But once the length reaches N (which happens very quickly), it would be way too costly to compress w at the end of each iteration.
However, only the values in W up to N/p will be used as factors in the next lot of deletions.
So it suffices to compress only this initial section of w, then pad on the right with zeros.
When the remaining primes are gathered on completion, it is necessary to skip zero and deleted values in w.

Each low-level operation in the resulting algorithm can be associated with an abstract operation so that each of the latter gets O(1) operations.
So the resulting program still runs in O(N/log log N) time, and the implicit constant factor is quite small.
It should be more than competitive with a good implementation of the classic sieve of Eratosthenes (i.e., without wheels or segmentation).

### wheel_sieve.cpp

This is a self-contained C++ program, which dynamically allocates the growing array w[].
It is obtained from wheel_sieve_simple.cpp by minimizing the compressions performed.

First notice that once the length of W reaches N, it is not necessary to delete the current prime:
it is sufficient to begin deletions with the square of the current prime.
Now, after the Delete step, the values in w will have no deletions up to p^2.
So if p'*p^2 exceeds N (where p' is the next prime), no compression at all is needed.
Furthermore, compression that is needed can start with p^2.

Also, the code is engineered for high-performance.
Function Delete uses SIMD (short vector) processing.
Function Compress is optimized to avoid slow accesses of the bit array d.
It instead searches for successive deleted future factors, and uses memcpy to move the sections between them.
In order to permit fast doubling (then binary) searches,
gaps after compression are filled with an even number (rather than 0) that keeps w ordered.

The resulting code is still essentially a faithful implementation of the abstract algorithm.
A detailed high-level code animation (using 0 for padding) for N=2000 is shown in this video:

[<img src="https://user-images.githubusercontent.com/1209656/200105831-ca678d1f-eaab-4895-8dea-58f04001211f.jpg" width="50%">](https://www.youtube.com/watch?v=-q-CFXLRSY0 "Code animation for N=2000")

### wheel_sieve_fast_doubling.cpp

This is a self-contained C++ program.
It is obtained from wheel_sieve.cpp by using doubling then binary searches to compute key
indexes used for compression and deletion.
As shown by the "Other" timing category, this has almost no effect on performance,
so the increased code complexity is not justified.

### wheel_sieve_fast_counter.cpp

This is a self-contained C++ program, which behaves like wheel_sieve_fast.cpp when the printing option is chosen.
However, when only a count of the primes up to N is required, it runs much faster.
It does so by essentially simulating the fast implementation of the dynamic wheel sieve
with a cut-down version using much smaller arrays.

The number of primes is simply the difference between the total number of values appended to w,
and the total number of (non-prime) values deleted from w, plus 1 for the prime 2.
These two totals are increased by the functions Extend and Delete respectively.
This observation already makes it unnecessary to scan the final array to count the primes.

But it is possible to do much better: to compute the increments without doing all the work necessary to create the primes.

Computing the number of values appended is straightforward.
When the k'th wheel is extended to the next wheel with p_(k+1) times the circumference,
the number of values appended is p_(k+1)-1 times the current number of values P_k
(assuming the new length <= N, and it can otherwise be easily computed).

And when multiples of a prime p are deleted, the number of deletions is determined by the values in w up to N/p.

So the idea is to avoid creating a new wheel when the next length (p_(k+1) times the current length P_k) would exceed N/p_(k+1).
Instead, the current wheel is simply extended (if necessary) up to a new limit N/p_(k+1).
Then, the sieve is performed as usual with the limit N/p_(k+1),
which is sufficient to generate all the factors needed for deletions up to N (since they are at most N/p for each prime p >= p_(k+1)).

The main complication is that since compression is not performed up to N/p, it is not easy to count the factors up to N/p.
(If compression were performed up to N/p instead of N/p^2, the cost would be proportional to deletions needed for
the original limit N, for primes p up to N^(1/3), which is O(N/\log \log N), as for the full sieve).

So we need an alternative way of efficiently computing the number of factors up to N/p.
This is just the number of values when W_(k+1) is rolled up to N/p,
minus the number of deletions up to N/p,
minus the number of values before p (which are 1 and the primes from p_(k+1)).
The crucial part is the number of deletions.
But these are recorded in the bit-array d[].
It would be too slow to scan d each time to count them.
But suppose we split the interval 1 up to limit=N/p_{k+1} into sqrt(limit) intervals of sqrt(limit) values,
and maintain counts of the number of deletions in each interval.
Then the count of deletions up to N/p can be computed in O(sqrt(limit)) time.

It follows from the original complexity analysis that p_(k+1) is proportional to log N.
Thus even if we computed the count for each prime p up to sqrt(limit),
the total cost would be O(N/log^2 N).
(We can do much better, since for p > N^{1/3}, no more deletions are performed, since p^2 > N/p.)

So, when only counting is requested, the work performed is dominated by the cost of sifting up to limit,
which is O(N/(log N\*log log N)), a sppedup by a factor of log N.
This time-complexity is bettered only by the specialised algorithms for prime counting dating from Meissel's in 1870.
The space required for array w is O(N/(log N\*log log N)) words, and for d is O(N/log N) bits.

### wheel_sieve_fast_counter_two.cpp

[STOP PRESS: I think this can be improved by the counting techniqueused above. TBC ]

This is a self-contained C++ program, which behaves like wheel_sieve_fast.cpp when the printing option is chosen,
except that the code is modified to also return an array of the primes found.
However, when only a count of the primes up to N is required, it runs much faster.
Its performance is comparable to that of wheel_sieve_fast_counter.cpp, but may be improved by using the same counting techniques.

Additional information is maintained by this algorithm.
First, the primes up to sqrt(N) are found (by a call to the Sift function).
Second, an array is maintained that essentially records the number of elements in W up to N/p for each prime p up to sqrt(N).
Also, the code has been refactored to better reveal the different stages of the main loop over the primes.

When only counting is requested, the work performed is dominated by compression,
which is proportional to the number of deletions performed by the dynamic wheel sieve with a limit of N/log N,
limited to primes up to N^(1/3).
The space for the additional arrays is also (comparatively) insignificant.
Accordingly, the time-complexity is O(N/(log N\*log log N)) operations.
The space required for array w is O(N/((log N)^2\*log log N)) words,
and for d is O(N/((log N)^2) bytes (which again could be bits).

## optimizations

There are many opportunities for optimizing the code to speed it up (hopefully!) or reduce the memory required.
Here are some:

- Array w[] can store the gaps betwen successive elements rather than the elements themselves.
It can then be very efficiently extended using memcpy.
However, searches for values in w[] will have to be linear.

- Array d[] can be implemented as a true bit-array.
Furthermore, it can be compressed when extended by using the pre-extended wheel in w[].
This automates the schemes used for highly-optimized versions of the sieve of Eratosthenes,
which are typically done by precomputing constant arrays.

- There are opportunities to use multiprocessing.

Such optimizations inevitably reduce the readability of the code.
Given that the practical utility of the dynamic wheel sieve of Pritchard is limited by its memory requirement,
and that its performance for its range of applicability is already excellent,
there seems little reason to pursue them.