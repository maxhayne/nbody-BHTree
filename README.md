## Build on Mac

> g++-10 -fopenmp -O2 -Wall Driver.cpp QuadTree.cpp

## Build on Linux

> g++ -fopenmp -O2 -Wall Driver.cpp QuadTree.cpp

## To run

> ./a.out

> python animate.py

Configure number of bodies, mass distribution, etc. in the driver file.

./a.out should produce a data file containing the location histories of all bodies. 'python animate.py' will produce an mp4 file of the simulation.

https://user-images.githubusercontent.com/32202629/166587947-fdbef902-22c5-4d1c-9975-8d8611701b46.mp4
