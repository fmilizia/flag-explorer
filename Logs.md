# About log files

While it runs, the program prints to the standard output some lines describing what happened during the run, with information like the following:
 - A certain triangulation has been built;
 - A certain triangulation has been analyzed, that is, the program has checked the existence of simplicial maps from the triangulation to the two target triangulations (for details, read the relevant sections in the paper).

It is adviseable to redirect the standard output to a file, so that the information it is filled with can later be used to reconstruct the history of previous runs of the program.



## What information is kept during execution

During a run two things happen:
 - Triangulations are created;
 - Triangulations are analyzed.

The sequence of creations during a run should be thought of as steps of a random walk in the space of triangulations.
We call **event** each of the steps in this random walk; that is, an event corresponds to the "discovery" of a triangulation (which, however, could be isomorphic to a triangulation built in a previous event; stepping on already known triangulations counts as an event as well).

A **set of encountered triangulations** is kept in memory during execution.
When the constructed triangulation is new, not isomorphic to previously found ones (read [here](Isomorphisms.md) to learn how we manage to efficiently detect isomorphisms among millions of triangulations), it is inserted into the set.

The result of an **analysis** is a bitset containing 3 bits.
Each of these bits tells whether there is a simplicial map induced by a map from a LocalPicture of the triangulation to one of the 3 "target" LocalPictures (which come from the 2 target triangulations with 10 and 12 vertices; check the paper for details).

This is what is kept in memory:
 - The set of encountered (pairwise nonisomorphic) triangulations;
 - A list of events, i.e., a list of triangulations, but concretely this is implemented as a list of pointers to elements of the set of encountered triangulations;
 - Additional information attached to every element of the set of encountered triangulations, including the result of its analysis (if it has been performed), and an `ID` which is an identification number associated to the triangulation;
 - A "reverse id" map that keeps for every ID a pointer to the corresponding element in the set of encountered triangulations.

In the current implementation, as soon as a new triangulation is discovered, it is also analyzed (and an unused ID is associated to it); by changing the source code one could separate these two processes, or even parallelize them.



## Structure of a log file

A log file is a textfile composed by lines of two different types:
 - Log lines, with a precise format explained below;
 - Comment lines, which are lines starting with the character #, followed by any text.

While running, the program prints on the standard output only log lines, not comment lines.
Comment lines can be added by anyone who wants to add annotations to a log file without changing the "history" that the log file describes.



## Format of a log line

A log line begins with a character (an uppercase letter of the alphabet) telling what kind of information the line encodes; there are, in fact, 7 types of log lines, corresponding to letters from A to G.
The information following the first character depends on the type of the log line, as summarized in the table below, where:
 - `id` is a nonzero integer, the identification number of a triangulation;
 - `dom` is a bitset represented as three 0/1 characters;
 - `E` is a vector of integers that encodes a triangulations (read [here](Isomorphisms.md) for more about the encoding of triangulations).


| Format of the line | Meaning |
| ------------------ | ------- |
| A `id` `E`         | Associate the identification number `id` to the triangulation encoded by `E` |
| B `id`             | Event of finding the triangulation corresponding to the given id |
| C `dom` `id`       | Set equal to `dom` the bitset associated to the triangulation with the given id | 
| D `id` `E`         | Combination of A and B |
| E `dom` `id`       | Combination of B and C |
| F `dom` `id` `E`   | Combination of A and C |
| G `dom` `id` `E`   | Combination of A, B and C |

The log lines of types A, B and C, respectively, assign an ID to a triangulation, encode the event of finding a triangulation, and tell the result of an analysis performed on a triangulation.
The other types are combinations of these three basic types.
The current implementation only prints lines of type G and B, when "stepping" into new triangulations (and the triangulation is then immediately analyzed) and old triangulations, respectively. 



## Importing a log file

One or more log files (formatted as explained above) can be read at the start of a run by adjusting the `main` function in [flag-explorer.cpp](source/flag-explorer.cpp), invoking the function `ReadLogFile` as many times as needed.

Comment lines are ignored and log lines are interpreted as explained above.
When an ID is already associated to a triangulation, but is then associated (with a log line of type A, D, F or G) to a different triangulation, the "old" association is completely discarded, that is, the previously associated triangulation loses its ID.
This can happen when reading log files that where created by unrelated executions of the program.

