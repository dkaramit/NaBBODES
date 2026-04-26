// This is how you run RK. 
#include<iostream>
#include<stdio.h>
#include<cmath>
#include"RKF.hpp"
#include"METHOD.hpp"


#ifndef LONG
#define LONG 
#endif


#define LD LONG double


#ifndef METHOD
#define METHOD RKF45
#endif


// this is how the diffeq should look like
#define n_eqs 3 //number of equations
using Array =  std::array<LD,n_eqs>;//define an array type of length n_eqs
//-------------------------------------------------------------------------//

using std::pow;

// you can use a function, but with a class you can also hold data that can be useful.
class diffeq{
    public:
    diffeq()=default;
    ~diffeq()=default;

    void operator()(Array &lhs, const Array &y, const LD& t){
        lhs[0]=-20*y[0]*pow(t,2) ;
        lhs[1]=5*y[0]*pow(t,2)+2*(-pow( y[1],2  )+pow( y[2],2 ) )*pow(t,1);
        lhs[2]=15*y[0]*pow(t,2)+2*(pow( y[1],2  )-pow( y[2],2 ) )*pow(t,1);
    }

};

// choose step controller (if you don't choose, it will use PI by default)
// using SOLVER = RKF::Solver<n_eqs,METHOD<LD>,LD, RKF::step_controlers::PI>;
using SOLVER = RKF::Solver<n_eqs,METHOD<LD>,LD, RKF::step_controlers::simple>;

int main(int argc, const char** argv) {
    Array y0 = {8,12,4};

    diffeq dydt;

    SOLVER System(dydt,y0,1e1,
        {.initial_step_size=1e-2,
        .minimum_step_size=1e-5,
        .maximum_step_size=1e1,
        .maximum_No_steps=1000000,
        .absolute_tolerance=1e-15,
        .relative_tolerance=1e-15,
        .beta=0.95,
        .fac_max=1.1,
        .fac_min=0.8
        }
    );
    
    System.solve();

    int step=0;
    for (auto _t: System.get_t()){
        printf("%e ",(double)_t);
        for( unsigned int eq = 0; eq < n_eqs; eq++){ printf("%e ", (double)System.get_solution(eq,step));}
        for( unsigned int eq = 0; eq < n_eqs; eq++){ printf("%e " ,(double)System.get_error(eq,step));}
        printf("\n");
        step++;
    }




    return 0;
 }