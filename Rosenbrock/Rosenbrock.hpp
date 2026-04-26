#ifndef Ros_headers
#define Ros_headers

#include "LU/LU.hpp"
//You need this for LU decomposition. 
// I have not seen a reason for this to be changed.
template<class LD> constexpr LD _tiny=1e-25;

#include"Jacobian/Jacobian.hpp"

#include"METHOD.hpp"

#include "Ros_class.hpp"
#include "Ros_reset.hpp"
#include "Ros_LU.hpp"
#include "Ros_calc_k.hpp"
#include "Ros_sums.hpp"
#include "Ros_step_control_simple.hpp"
#include "Ros_step_control_PI.hpp"
#include "Ros_steps.hpp"



#endif