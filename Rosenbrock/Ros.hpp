#ifndef Ros_headers
#define Ros_headers
#include <vector>
#include "LU/LU.hpp"




/*
diffeq is a class of the system of  equations to be solved 
N_eqs is ten number of equations to be solved
RK_method is the method (I trust ROS34PW2)
N_out number of output points (to be taken in intervals of approximately  1/(N_out-1) )
*/

#define _Ros_template_ template<class diffeq, int N_eqs, class RK_method, class jacobian, int N_out, class LD> 

#include "Ros_class.hpp"

#define _Ros_Func_ void  Ros<diffeq, N_eqs,  RK_method,jacobian, N_out, LD>
#define _Ros_Cosnt_ Ros<diffeq, N_eqs, RK_method,  jacobian, N_out, LD>



#include "Ros_costructor.hpp"
#include "Jacobian.hpp"
#include "Ros_LU.hpp"
#include "Ros_calc_k.hpp"
#include "Ros_sums.hpp"
// #include "Ros_step_control-simple.hpp"
#include "Ros_step_control-PI.hpp"

#include "Ros_steps.hpp"


#endif