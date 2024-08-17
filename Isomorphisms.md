# Isomorphisms between triangulations



The program manages to gracefully handle millions of triangulations of $S^3$ and to check quite efficiently, given a triangulation, whether it is isomorphic to one of those.
How does it perform this task?



## Checking if two triangulations are isomorphic

The first thing to notice is that, given two triangulations, it is rather easy to check if they are isomorphic or not.
Here we are speaking, in particular, of (flag) triangulations of $S^3$, but the same holds in general for pseudomanifolds.

The idea is simple.
Assume that you have two triangulations $S$ and $T$ (in general, strongly connected pseudomanifolds of the same dimension), and you want to find an isomorphism from $S$ to $T$, if it exists.
Then a friend arrives and, somehow, tells you with confidence that you should map the vertices of one of the top-dimensional simplices of S to some specific vertices of $T$ (which, of course, span a simplex of the same dimension).
You then try to extend the partial map suggested by your friend to an isomorphism from $S$ to $T$.
You don't need much ingenuity to do that! Everything is forced: look at a top-dimensional simplex adjacent (sharing a facet) to the first one on which the map is already defined; it must go to the corresponding simplex in $T$, adjacent to the image of the first simplex.
So, you know how to extend the map to another vertex. And you continue in the same fashion, until you extend the map to an entire isomorphism or, if your friend was wrong, you find a "contradiction" at some point.
How fast you can actually perform this construction depends on how you are storing your simplicial complexes, but you can definitely arrange to do it in polynomial time in the sizes (number of simplices) of the complexes.

If you don't have a helpful friend telling you how to start defining your map, you can just try all possible "starting partial maps", and you get a quite efficient algorithm anyway.



## Keeping a set of pairwise-nonisomorphic triangulations

Suppose now that we are given a big quantity of triangulations, and our task is to subdivide them into isomorphism classes (or to extract a representative for each isomorphism class).
How can we do it?

Of course, one option is to perform an isomorphism-check for each pair of triangulations.
But this amount to performing a number of checks which is quadratic in the number of triangulations; could be feasible if you are dealing with thousands of them, not so much if it is _millions_.

We introduce two new ideas: labelled triangulations and encodings.

### Labelled triangulations and encodings

For us, a labelled triangulation is a triangulation together with an ordering of its vertices. In other words, vertices are numbered with integers, from $1$ to $N$, where $N$ is the number of vertices of the triangulation.
Note that a triangulation can be labelled in only finitely many ways (but they are a lot: $N!$).

Then, we want to "encode" labelled triangulations using some easy object; for instance, finite lists of integer numbers.
Since we are dealing with flag triangulations, which are uniquely determined by their 1-skeleton, we can use the following encoding scheme: the encoding is a finite list of integers, the first one being $N$, the number of vertices, and the remaining ones encoding the 1-simplices.
For instance, a 1-simplex with vertices labelled $u$ and $v$, with $u < v$, can be represented by the integer $(u-1)*N+(v-1)$.
We can list the numbers encoding the edges in increasing order, so that the encoding is uniquely dertermined from the labelled triangulation.

## Encoding an unlabelled triangulation

We wish to "encode" unlabelled triangulations, in such a way that the encoding is invariant under isomorphisms of triangulations.
We want to use encoding objects which are easy to compare, so that we can decide quickly whether two encodings are equal, or which is the "smallest".
Lists of integers are good: two of them can be compared lexicographically.

We start from this idea: consider all the labellings of a triangulation, computer their encodings, and keep the smallest one.
This realizes the wish expressed above, but computing this encoding takes too much time, because there are too many labellings.
The next idea is to restrict to a specific, much smaller set of "allowed" labellings, so that we can actually perform efficiently our "keep the smallest encoding" strategy.

One solution is very similar to what we did in the section where we discussed how to check whether two triangulations are isomorphic:
 - Start with a "partial labelling", defined only on the vertices of a top-dimensional simplex, assigning to them the smallest labels starting from $1$;
 - Continue with some fixed deterministic way to "visit" the rest of the vertices. For instance, take the smallest facet of one of the "already visited" top-dimensional simplices (such facets can be encoded by list of integers, since their vertices are already labelled, so they can be compared and ordered) which is adjacent to a "not yet visited" simplex, and assign the next label to the only new vertex in that adjacent "not yet visited" simplex.

The number of these "special encodings" is now much smaller than $N!$, providing us with an efficient way to encode triangulations in a isomorphism-invariant way.

Once we have such an encoding strategy, to accomplish our task we can:
 - Compute the encodings of all the triangulations;
 - Order the encodings, so that it is then immediate to spot equal ones or pick representatives.

The second step can be performed by executing a number of comparisons between encodings that grows as $L \cdot \log(L)$, where $L$ is the number of triangulations, which is significantly less than $L^2$ and is feasible even when $L$ is of the order of millions.

By using well-known data structures, one can keep a set of encodings (representing a set of pairwise-nonisomorphic triangulations) and, given another triangulation:
 - Compute the encoding of the triangulation;
 - Check efficiently if the encoding is already in the set, and possibly insert it in the set if it is new.

