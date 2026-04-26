// This is how you run Rosenbrock, using a functor
#include<iostream>
#include<fstream>
#include<cmath>
#include"Ros.hpp"
// #include "Jacobian/Jacobian.hpp"//this is treated as user input, since one may have an analytic form.
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
    diffeq(const LD& c):c(c){};

    void operator()(Array &lhs, const Array &y  , const LD& t)const{
        lhs[0]=t*c;
    }

};

using SOLVER = Ros<n_eqs, METHOD<LD>, LD>;

int main(int argc, const char** argv){
    
    Array y0 = {2};
    diffeq dydt(2);

    // this is the numerical Jacobian
    // Jacobian<n_eqs,LD> Jac(dydt,1e-10);
    
    //we know the analytical jacobian here (type deduction is powerful :P)
    SOLVER::Jacobian Jac=[&dydt](auto& J, auto& dfdt,auto& y, auto& t){ dfdt[0]=dydt.c; J.fill({0}); };

    SOLVER System(dydt,y0, 1e4, Jac,
        {
            .initial_step_size=1e-2,
            .minimum_step_size=1e-8,
            .maximum_step_size=1e3,
            .absolute_tolerance=1e-8,
            .relative_tolerance=1e-8,
            .beta=0.5,
            .fac_max=1.2,
            .fac_min=0.6
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

    //this is how you can change parameters
    System.reset(y0,1e4,{.initial_step_size=20,.absolute_tolerance=1e-2,.relative_tolerance=1e-2});
    System.solve();
    
    step=0;
    step=System.get_current_step()-1;
    printf("%e ",(double)System.get_t().at(step));
    for( unsigned int eq = 0; eq < n_eqs; eq++){ printf("%e ", (double)System.get_solution(eq,step));}
    for( unsigned int eq = 0; eq < n_eqs; eq++){ printf("%e " ,(double)System.get_error(eq,step));}
    printf("\n");
    std::cout<<System.get_current_step();
    printf("\n");

    std::cout<<System.get_parameters().initial_step_size.value()<<std::endl;
    std::cout<<System.get_parameters().minimum_step_size.value()<<std::endl;
    

    
    return 0;
 }