// This is how you run RK. 
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
#define n_eqs 3 //number of equations
using Array =  std::array<LD, n_eqs>;//define an array type of length n_eqs
//-------------------------------------------------------------------------//

using std::pow;


// you can use a function, but with a class you can also hold data that can be useful.
class diffeq{
    public:
    diffeq(){};
    ~diffeq(){};

    void operator()(Array &lhs, const Array &y  , const LD& t){
        lhs[0]=-2*y[0]*pow(t,2) ;
        lhs[1]=2*y[0]*pow(t,2)+2*(-pow( y[1],2  )+pow( y[2],2 ) )*pow(t,1);
        lhs[2]=4*y[0]*pow(t,2)+2*(pow( y[1],2  )-pow( y[2],2 ) )*pow(t,1);
    }

};



using SOLVER = Ros<n_eqs, METHOD<LD>, LD>;

int main(int argc, const char** argv) {
    
    Array y0 = {8,12,4};
    diffeq dydt;
    
    // we use the default Jacobian. The final parameter in the constructor is the Jacobian's value for h
    SOLVER System(dydt,y0, 1e4,
        {
        .initial_step_size = 1e-2,
        .minimum_step_size = 1e-8,
        .maximum_step_size = 1e4,
        .maximum_No_steps = 1000000,
        .absolute_tolerance = 1e-8,
        .relative_tolerance = 1e-8,
        .beta = 0.95,
        .fac_max = 1.1,
        .fac_min = 0.7
        },
        1e-4
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