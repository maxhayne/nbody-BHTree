## To Compile

Apple's clang g++ doesn't have support for the '-fopenmp' option, so I've installed the GNU version through homebrew:

~~~console
brew install gcc
~~~

And then to compile:
~~~console
g++-11 -fopenmp -O3 -Wall Driver.cpp QuadTree.cpp
~~~
## on Linux
~~~console
g++ -fopenmp -O3 -Wall Driver.cpp QuadTree.cpp
~~~

## To run
~~~console
./a.out
python3 animate.py
~~~

Configure number of bodies, mass distribution, etc. in the driver file.

*./a.out* should produce a data file containing the location histories of all bodies. *python animate.py* will produce an mp4 file of the simulation.

### Three-body integration

https://user-images.githubusercontent.com/32202629/166587947-fdbef902-22c5-4d1c-9975-8d8611701b46.mp4

### Many-body *galaxy* integration

https://user-images.githubusercontent.com/32202629/166634755-c86ca09b-599d-41f7-8d2e-0f1c1e761f59.mp4
