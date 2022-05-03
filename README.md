## Build on Mac

> g++-10 -fopenmp -O2 -Wall Driver.cpp QuadTree.cpp

## Build on Linux

> g++ -fopenmp -O2 -Wall Driver.cpp QuadTree.cpp

## To run

> ./a.out

> python animate.py

Configure number of bodies, mass distribution, etc. in the driver file.

./a.out should produce a data file containing the location histories of all bodies. 'python animate.py' will produce an mp4 file of the simulation.

### Three-body integration

https://user-images.githubusercontent.com/32202629/166587947-fdbef902-22c5-4d1c-9975-8d8611701b46.mp4

### Many-body *galaxy* integration

https://user-images.githubusercontent.com/32202629/166588398-c2e433ee-9a7c-4636-8b4a-df5b1b88e53d.mp4
