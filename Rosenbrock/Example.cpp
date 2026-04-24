// This is how you run Rosenbrock, using a functor
#include<iostream>
#include<fstream>
#include<cmath>
#include"Ros.hpp"
#include "Jacobian/Jacobian.hpp"//this is treated as user input, since one may have an analytic form.
#include "METHOD.hpp"

using std::cout;
using std::endl;
/*--------------------------------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------------------------------*/
#ifndef LONG
#define LONG 
#endif


#define LD LONG double


#ifndef METHOD
#define METHOD ROS34PW2
#endif


// this is how the diffeq should look like
#define n_eqs 1 //number of equations
using Array =  std::array<LD, n_eqs>;//define an array type of length n_eqs
//-------------------------------------------------------------------------//


// you can use a function, but with a class you can also hold data that can be useful.
class diffeq{
    public:
    LD c;
    diffeq(LD c):c(c){};

    void operator()(Array &lhs, Array &y  , LD t){
        lhs[0]=t*c;
    }

};



using SOLVER = Ros<n_eqs, METHOD<LD> ,Jacobian<n_eqs,LD> , LD>;

int main(int argc, const char** argv) {
    
    Array y0 = {2};
    diffeq dydt(2);


    SOLVER System(dydt,y0, 1e4,
        {
            .initial_step_size=1e-2,
            .minimum_step_size=1e-8,
            .maximum_step_size=1e3,
            .maximum_No_steps=1000000,
            .absolute_tolerance=1e-8,
            .relative_tolerance=1e-8,
            .beta=0.5,
            .fac_max=1.01,
            .fac_min=0.9
        }
    );
    
    System.solve();

    int step=0;
    for (auto _t: System.time){
        printf("%e ",(double)_t);
        for( unsigned int eq = 0; eq < n_eqs; eq++){ printf("%e ", (double)System.solution[eq][step]);}
        for( unsigned int eq = 0; eq < n_eqs; eq++){ printf("%e " ,(double)System.error[eq][step]);}
        printf("\n");
        step++;
    }

    
    return 0;
 }