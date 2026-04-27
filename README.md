

# NaBBODES
*Not a Black Box ODE Solver*

![](https://img.shields.io/badge/language-C++-black.svg)  ![](https://tokei.rs/b1/github/dkaramit/NaBBODES)
  
![](https://img.shields.io/github/repo-size/dkaramit/NaBBODES?color=blue)

ODE solvers that anyone should be able to take, hack, and use for specific problems they might be interested in solving. The plan is to be able to use the skeleton of RK methods, and be able to pass as arguments everything from just the Butcher tableau to specific step control, and optimizing every step for a chosen RK method. Examples on how to use the code are given, although documentation should be included.

There is some rudimentary documentation, that will be expanded in the future:

1. Read a bit about the use of RKF and Rosenbock in [Usage.md](Usage.md).

2. Dive a bit deeper if you are interested in hacking the code in [Hacking.md](Hacking.md).

3. If you are interested in seeing some basic maths (and conventions) used in RKF and Rosenbock methods go to [BackgroundMaths](BackgroundMaths), run the `make` (to create the Basics.pdf from the Basics.tex file), and read! 

If you want to make suggestion, contribute, or you can't make the code work, contact me at <dkaramit@yahoo.com> or <dimitrios.karamitros@pd.infn.it>. 

---
NaBBODES started as part of [*ASAP*](https://dkaramit.github.io/ASAP/), but it evolved to an independent thing. Here we will be able to further improve it without keeping things as simple as possible.  That said, the code should still be kept somewhat simple to allow the user to understand and modify it.

That's it for now.

Enjoy,

Dimitris
