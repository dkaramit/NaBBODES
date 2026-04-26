#ifndef RKF_step_control_simple
#define RKF_step_control_simple
#include "RKF_class.hpp"


namespace rkf{
    
/*-----------------------Begin: step_control---------------------------------*/

template<unsigned int N_eqs, class RK_method, class LD, step_controllers step_controller>
void Solver<N_eqs, RK_method, LD, step_controller>::step_control_simple(){
    LD Delta=0.;
    LD _sc;
    LD fac=beta;

    for (unsigned int eq = 0; eq < N_eqs; eq++){
        _sc=std::max(std::abs( ynext[eq] ), std::abs( yprev[eq] ));
        _sc=abs_tol+rel_tol*_sc;
        Delta+= (abs_delta[eq]/_sc)*(abs_delta[eq]/_sc);  
    ;}

    Delta=std::sqrt(1./N_eqs*Delta);
    if(Delta<1) { 
        delta_acc=Delta;
        h_stop=true; 
    }

    //step size cotrol from "Solving Ordinary Differential Equations I"
    fac*=std::pow( Delta , -1./( static_cast<LD>(RK_method::p + 1)));
    if(fac> fac_max){fac = fac_max;}
    if(fac< fac_min){fac = fac_min;}
    
    h_trial= h_trial*fac ;

    if (h_trial>hmax ){ h_trial=hmax; h_stop=true;}
    if (h_trial<hmin ){ h_trial=hmin; h_stop=true;}
}
/*-----------------------End: step_control---------------------------------*/

}

#endif
