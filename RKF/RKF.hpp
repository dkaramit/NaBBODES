#ifndef RKF_headers
#define RKF_headers
#include <vector>




/*
diffeq is a class of the system of  equations to be solved 
N_eqs is ten number of equations to be solved
RKF_method is the method (the DormandPrince seems to be the standard here)
N_out number of output points (to be taken in intervals of approximately  1/(N_out-1) )
*/

#define _RKF_template_ template<class diffeq, int N_eqs, class RKF_method, int N_out, class LD> 

#include "RKF_class.hpp"

#define _RKF_Func_ void  RKF<diffeq, N_eqs, RKF_method, N_out, LD>
#define _RKF_Cosnt_ RKF<diffeq, N_eqs, RKF_method, N_out, LD>


#include "RKF_costructor.hpp"
#include "RKF_calc_k.hpp"
#include "RKF_sums.hpp"
#include "RKF_step_control.hpp"
#include "RKF_steps.hpp"


#endif