# NaBBODES
*Not a Black Box ODE Solver*

A package for ODE solvers that anyone should be able to take, hack, and use for specific problems they might be interested in solving. The plan is to be able to use the skeleton of RK methods, and be able to pass as arguments, everything from just the Butcher tableau to specific step control, and optimizing every step for a chosen RK method. 


---
This repository exists, in order to package the ODE solvers from [*ASAP*](https://dkaramit.github.io/ASAP/).

*ASAP* is a repository in which we try to have simple code in order to show how everything works, having in mind that someone who is not an expert in C++ or python can understand and translate each algorithm to a familiar language. In *NaBBODES* we don't really care about that. The code will still be simple (to avoid maintenance cost), but here we are going to use features of C++ and python, to make the code feel more professional and be faster.   

The plan for the moment is to start with the C++ ODE solvers from *ASAP*, and put them together. Later, we'll try to optimize them, and write some kind of wrapper for python (ctypes or pybind11).


That's it for now.

Enjoy,
Dimitris
