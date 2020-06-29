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


#include "Ros_class.hpp"
#include "Ros_costructor.hpp"
#include "Jacobian.hpp"
#include "Ros_LU.hpp"
#include "Ros_calc_k.hpp"
#include "Ros_sums.hpp"
// #include "Ros_step_control-simple.hpp"
#include "Ros_step_control-PI.hpp"
#include "Ros_steps.hpp"

#endif