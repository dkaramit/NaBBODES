// This is how you run RK. 
#include<iostream>
#include<stdio.h>
#include<cmath>
#include"RKF.hpp"


#ifndef LONG
#define LONG 
#endif


#define LD LONG double


#ifndef METHOD
#define METHOD RKF45
#endif



using std::pow;

// you can use a function, but with a class you can also hold data that can be useful.
class diffeq{
    public:
    diffeq()=default;
    ~diffeq()=default;

    void operator()(std::vector<LD> &lhs, const std::vector<LD> &y, const LD& t){
        lhs[0]=-20*y[0]*pow(t,2) ;
        lhs[1]=5*y[0]*pow(t,2)+2*(-pow( y[1],2  )+pow( y[2],2 ) )*pow(t,1);
        lhs[2]=15*y[0]*pow(t,2)+2*(pow( y[1],2  )-pow( y[2],2 ) )*pow(t,1);
    }

};

// using SOLVER = rkf::Solver<LD,rkf::METHOD<LD>, rkf::step_controllers::PI>;
using SOLVER = rkf::Solver<LD, rkf::METHOD<LD>, rkf::step_controllers::simple>;

// this also works with DormandPrince as default method and PI as step_controller.
// using SOLVER = rkf::Solver<LD>;

int main(int argc, const char** argv) {
    std::vector<LD> y0 = {8,12,4};

    diffeq dydt;

    SOLVER System(dydt,y0,1e1,
        {.initial_step_size=1e-2,
        .minimum_step_size=1e-5,
        .maximum_step_size=1e2,
        .maximum_No_steps=1000000,
        .absolute_tolerance=1e-8,
        .relative_tolerance=1e-8,
        .beta=0.95,
        .fac_max=1.5,
        .fac_min=0.5
        }
    );
    
    System.solve();

    unsigned int n_eqs=y0.size();
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