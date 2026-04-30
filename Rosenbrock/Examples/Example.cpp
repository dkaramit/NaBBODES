// This is how you run Rosenbrock, using a functor
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


// you can use a function, but with a class you can also hold data that can be useful.
class diffeq{
    public:
    LD c;
    diffeq(const LD& c):c(c){};

    void operator()(std::vector<LD> &lhs, const std::vector<LD> &y  , const LD& t)const{
        lhs[0]=t*c;
    }

};

// choose step controller (if you don't choose, it will use PI by default)
// using SOLVER = rosenbrock::Solver<LD, rosenbrock::METHOD<LD>, rosenbrock::step_controllers::PI>;
// using SOLVER = rosenbrock::Solver<LD, rosenbrock::METHOD<LD>, rosenbrock::step_controllers::simple>;
using SOLVER = rosenbrock::Solver<LD>;

int main(int argc, const char** argv){
    
    std::vector<LD> y0 = {2};
    diffeq dydt(2);

    unsigned int n_eqs=y0.size();


    //we know the analytical jacobian here (type deduction is powerful :P)
    SOLVER::Jacobian_type Jac=[&dydt](auto& J, auto& dfdt,auto& y, auto& t){ dfdt[0]=dydt.c; J.assign(1,std::vector<LD>(1,0)); };

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


    std::cout<<System.get_parameters().initial_step_size.value()<<std::endl;
    std::cout<<System.get_parameters().minimum_step_size.value()<<std::endl;
    std::cout<<System.get_parameters().maximum_No_steps.value()<<std::endl;

    //this is how you can change parameters
    System.reset(y0,1e4,{.initial_step_size=20,.maximum_No_steps=500,.absolute_tolerance=1e-2,.relative_tolerance=1e-2});
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
    std::cout<<System.get_parameters().maximum_No_steps.value()<<std::endl;
    

    
    return 0;
 }