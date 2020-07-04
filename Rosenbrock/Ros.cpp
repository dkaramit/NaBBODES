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


#define initial_step_size 1e-2
#define minimum_step_size 1e-8
#define maximum_step_size 1e4
#define maximum_No_steps 1000000
#define absolute_tolerance 1e-8
#define relative_tolerance 1e-8
#define beta 0.95
#define fac_max 1.1
#define fac_min 0.7

#define N_out 500
// this is how the diffeq should look like
#define n_eqs 3 //number of equations
using Array =  LD[n_eqs];//define an array type of length n_eqs
//-------------------------------------------------------------------------//

using std::pow;


// you can use a function, but with a class you can also hold data that can be useful.
class diffeq{
    public:
    diffeq(){};
    ~diffeq(){};

    void operator()(Array &lhs, Array &y  , LD t){
        lhs[0]=-2*y[0]*pow(t,2) ;
        lhs[1]=2*y[0]*pow(t,2)+2*(-pow( y[1],2  )+pow( y[2],2 ) )*pow(t,1);
        lhs[2]=4*y[0]*pow(t,2)+2*(pow( y[1],2  )-pow( y[2],2 ) )*pow(t,1);
    }

};



using SOLVER = Ros<diffeq,n_eqs, METHOD<LD> ,Jacobian<diffeq,n_eqs,LD>, N_out , LD>;

int main(int argc, const char** argv) {
    
    Array y0 = {8,12,4};
    diffeq dydt;

    Jacobian<diffeq,n_eqs,LD> jac(dydt);

    SOLVER System(dydt,y0, 1e4,
    initial_step_size,  minimum_step_size,  maximum_step_size, maximum_No_steps, 
    absolute_tolerance, relative_tolerance, beta, fac_max, fac_min);
    
    System.solve(true);

    std::cout<<System.time.size()<<"\n";
    // std::cout<<System.Deltas.size()<<"\n"; // this it the same as time_full
    std::cout<<System.time_full.size()<<"\n";

    // this prints only N_out even if one runs System.solve(true)
    int step=0;
    for (auto _t: System.time){
        printf("%e ",(double)_t);
        for( int eq = 0; eq < n_eqs; eq++){ printf("%e ", (double)System.solution[eq][step]);}
        for( int eq = 0; eq < n_eqs; eq++){ printf("%e " ,(double)System.error[eq][step]);}
        printf("%e\n" ,(double)System.hist[step]) ; 
        step++;
    }
    // print the deltas. Should be <~ 1
    for(auto _del : System.Deltas) {  printf("%e \n",(double)_del ) ; }


    step=0;
    // print the full solution
    for(auto _t: System.time_full ) {  
        printf("%e ",(double)_t) ; 
        for( int eq = 0; eq < n_eqs; eq++){ printf("%e ", (double)System.solution_full[eq][step]);    }
        for( int eq = 0; eq < n_eqs-1; eq++){ printf("%e ", (double)System.error_full[eq][step]);}
        printf("%e\n", (double)System.error_full[n_eqs-1][step]);
        ++step;
    }
    
    return 0;
 }