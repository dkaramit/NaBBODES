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

#define LD LONG double



// this is how the diffeq should look like
#define n_eqs 3 //number of equations
typedef LD Array[n_eqs];//define an array type of length n_eqs
typedef void (*diffeq)(Array &lhs, Array &y  , LD t );
//-------------------------------------------------------------------------//

void sys( Array &lhs, Array &y  , LD t )
        {
            //lhs is an array that gets the return value (the left hand side of the equation)
            //y is an array with values of y
            //t is the value of the variable t
            
            lhs[0]=-20*y[0]*pow(t,2) ;
            lhs[1]=5*y[0]*pow(t,2)+2*(-pow( y[1],2  )+pow( y[2],2 ) )*pow(t,1);
            lhs[2]=15*y[0]*pow(t,2)+2*(pow( y[1],2  )-pow( y[2],2 ) )*pow(t,1);

        };


#define initial_step_size 1e-5
#define minimum_step_size 1e-15
#define maximum_step_size 1e-3
#define maximum_No_steps 1000000
#define absolute_tolerance 1e-10
#define relative_tolerance 1e-10
#define beta 0.9
#define fac_max 100

#define N_out 500

//you can also define then using the -D flag ( as -DMETHOD=ROS34PW2 for example)
// #define METHOD ROS3w //2nd order
// #define METHOD ROS34PW2 //3rd order
// #define METHOD RODASPR2  //4th order

int main(int argc, const char** argv) {
    
    Array lhs;
    Array y0;
    y0[0]=5;
    y0[1]=10;
    y0[2]=0;
    Jacobian<diffeq,n_eqs,LD> jac(sys);

    // this is not looing good.. Just pass the number of equations as argument to get rid of Array and Matrix...
    Ros<diffeq,n_eqs, METHOD<LD> ,Jacobian<diffeq,n_eqs,LD>, N_out , LD > System(sys,y0, 
     initial_step_size,  minimum_step_size,  maximum_step_size, maximum_No_steps, 
     absolute_tolerance, relative_tolerance, beta, fac_max);
    // System.next_step();
System.solve(true);
    // 

    std::cout<<N_out<<"\n";
    std::cout<<System.Deltas.size()<<"\n";
    std::cout<<System.current_step<<"\n";
    

    for (int i = 0; i < N_out; i++){
        printf("%e ",(double)System.time[i]);

        for( int eq = 0; eq < n_eqs; eq++){ printf("%e ", (double)System.solution[eq][i]);    }
        for( int eq = 0; eq < n_eqs; eq++){ printf("%e " ,(double)System.error[eq][i]) ; }
        printf("%e\n" ,(double)System.hist[i]) ; 
    }
    

    for(int i=0; i< System.Deltas.size() ; ++i) {  printf("%e \n",(double)System.Deltas[i]) ; }

    for(int i=0; i< System.time_full.size() ; ++i) {  
        printf("%e ",(double)System.time_full[i]) ; 
        for( int eq = 0; eq < n_eqs-1; eq++){ printf("%e ", (double)System.solution_full[eq][i]);    }
        printf("%e\n", (double)System.solution_full[n_eqs-1][i]);
        }
    return 0;
 }