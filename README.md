# Flag explorer

This is an implementation of the algorithm described in the paper _Simplicial maps betweenn spheres and Davis' manifolds with positive simplicial volume_.
It is the C++ implementation which has been used to perform the experiments mentioned in the paper.



## What the code does

It builds, randomly, flag simplicial complexes homeomorphic to the 3-dimensional sphere, and then tries to find simplicial maps of nonzero degree going from these complexes to two specific triangulations of the 3-dimensional sphere.
These two specific triangulation have 10 and 12 vertices, respectively, and their 1-skeleta are described in the files contained in the [triangulations](triangulations/) directory.
More details are explained in the paper mentioned above.



## How to use the code

All the C++ code is located in the [source](source/) directory.
The file containing the `main` function is [flag-explorer.cpp](source/flag-explorer.cpp). All other files contain implementations of classes and functions used in the main file.

To compile, use, for instance, the following command:
```
g++ -Wall -Wextra -O2 -std=gnu++20 -o "flag-explorer" "source/flag-explorer.cpp"
```
This creates an executable file.
Notice that, at runtime, the program expects to find, under the working directory, the folder [localpictures](localpictures/), because it needs the files contained in it.
This behaviour can be changed by adjusting the relevant code in [flag-explorer.cpp](source/flag-explorer.cpp).

One piece of code that you might want to experiment with, trying some changes, is the function `Explore`, located in [flag-explorer.cpp](source/flag-explorer.cpp) just before the `main` function, which is at the end of the file.
By changing this function you can implement different strategies for the exploration of the space of triangulations.  
You can set the number of iterations performed in a run of the program in the `main` function, where the `Explore` function is called.

While running, the program produces a "history" of found triangulations, which is streamed on the standard output.
It is advisable to _redirect_ this output to a file, when executing the program.
The resulting "history file" can later be read at the beginning of a subsequent run; this behaviour can be activated by uncommenting the relevant line in the `main` function.

