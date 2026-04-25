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
#define n_eqs 1 //number of equations
using Array =  std::array<LD,n_eqs>;//define an array type of length n_eqs
//-------------------------------------------------------------------------//

/*std::function supports any callable object, so we cna simply use this without overloading 
operator=, or defining copy constructor.*/
class diffeq{
    public:
    LD c;
    diffeq(double c):c(c){}; 

    void operator()(Array &lhs, const Array &y, const LD& t){
        lhs[0]=t*c;
    }

};

using SOLVER = RKF<n_eqs,METHOD<LD>,LD>;

int main(int argc, const char** argv) {
    Array y0={0};

    diffeq dydt(2);

    SOLVER System(dydt,y0,1e1, {.initial_step_size=1e-5,.minimum_step_size=1e-5,.maximum_step_size=1e1});
    
    System.solve();

    int step=0;
    for (auto _t: System.get_t()){
        printf("%e ",(double)_t);
        for( unsigned int eq = 0; eq < n_eqs; eq++){ printf("%e ", (double)System.get_solution(eq,step));}
        for( unsigned int eq = 0; eq < n_eqs; eq++){ printf("%e " ,(double)System.get_error(eq,step));}
        printf("\n");
        step++;
    }

    System.reset(y0,1e1,{.initial_step_size=20,.absolute_tolerance=1e-2,.relative_tolerance=1e-2});
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