// This is how you run RK. 
#include<iostream>
#include<fstream>
#include<cmath>
#include"Rosenbrock.hpp"

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

using std::pow;

// you can use a function, but with a class you can also hold data that can be useful.
class diffeq{
    public:
    diffeq(){};
    ~diffeq(){};

    void operator()(std::vector<LD> &lhs, const std::vector<LD> &y  , const LD& t){
        lhs[0]=-2*y[0]*pow(t,2) ;
        lhs[1]=2*y[0]*pow(t,2)+2*(-pow( y[1],2  )+pow( y[2],2 ) )*pow(t,1);
        lhs[2]=4*y[0]*pow(t,2)+2*(pow( y[1],2  )-pow( y[2],2 ) )*pow(t,1);
    }

};



// choose step controller (if you don't choose, it will use PI by default)
using SOLVER = rosenbrock::Solver<LD, rosenbrock::METHOD<LD>, rosenbrock::step_controllers::PI>;
// using SOLVER = rosenbrock::Solver<LD, rosenbrock::METHOD<LD>, rosenbrock::step_controllers::simple>;

int main(int argc, const char** argv) {
    
    std::vector<LD> y0 = {8,12,4};
    diffeq dydt;
    
    unsigned int n_eqs=y0.size();

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