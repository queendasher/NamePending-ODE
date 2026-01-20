# Introduction to Namepending-ODE


ASC-ODE is is a C++ library for solving ordinary differential equations (ODEs).
The equation is defined by the right hand side function.
ASC-ODE provides various time-steppers which may be used for odes with right hand sides
given by a function object.

The theory can be found here: https://jschoeberl.github.io/IntroSC/ODEs/ODEs.html

## Installation

Install Namepending-ODE via git-clone:

    git clone https://github.com/queendasher/NamePending-ODE.git



For the python applications, make sure you have the necessary dependencies installed, such as 
- Python 3.12 (or higher)
- pybind11
- pythreejs
- ...

They can be installed using apt on Linux:

```sh
sudo apt install python3.12
sudo apt install python3-pybind11
``` 

To build the library, navigate to the `src` folder and execute the following commands:

```sh
mkdir build
cd build
cmake .. 
make
```

After the build is complete, you can run the demonstration scripts. There should be an executable named `test_odes` in the `build` directory.

To run the demonstration scripts, execute the following command from the `build` directory:
```sh
./test_ode
```


This program will print the simulation results into a text file `output_test_ode.txt` depending on which folder the program has been called. 

The results can be quite different depending on the used time-stepping method. Available time-stepping methods are: 
- ExplicitEuler
- ImplicitEuler
- ImprovedEuler
- RungeKutta
- ExplicitRungeKutta
- ImplicitRungeKutta
- CrankNicolson


To evaluate the contents of the text file, a python script called `plotmassspring.py` can be called. This script loads the contents and generates two plots containing the simulation data.

```sh
python3 ./plotmassspring.py
```



   
