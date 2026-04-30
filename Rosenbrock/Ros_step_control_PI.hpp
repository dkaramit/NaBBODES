#ifndef Ros_step_control_PI
#define Ros_step_control_PI
#include "Ros_class.hpp"


namespace rosenbrock{

// Keep in mind that here delta_acc=Deltas.back(), while 
// delta_rej is the previous Delta (not the accepted one).

/*-----------------------Begin: step_control---------------------------------*/
template<class LD, class RK_method, step_controllers step_controller> 
void Solver<LD, RK_method, step_controller>::step_control_PI(){
    LD Delta=0.;
    LD _sc;
    LD fac=beta;
    // regulated_delta = (1-gam*h*J)^{-1}*abs_delta
    std::vector<LD> regulated_delta(N_eqs);
    //PI step size cotrol from "Solving Ordinary Differential Equations II"
    //We rescale the error by multiplying it with (1 + h \gamma J)^-1
    dot<LD>( _inv, abs_delta , regulated_delta);

    for (unsigned int eq = 0; eq < N_eqs; eq++){
        _sc=std::max(std::abs( ynext[eq] ), std::abs( yprev[eq] ));
        _sc=abs_tol+rel_tol*_sc;
        Delta+= (regulated_delta[eq]/_sc)*(regulated_delta[eq]/_sc);
    }
    Delta=std::sqrt(1./N_eqs*Delta);
    if(Delta==0){Delta=abs_tol+rel_tol;}

    // I use the step control from 
    // https://www.sciencedirect.com/science/article/pii/S147466701751767X
    if(Delta<=1) { 
        if(delta_rej<=1){fac*=h_trial/h_old; }
        fac*=std::pow(Delta, -0.65/( static_cast<LD>(RK_method::p + 1)) );   
        fac*=std::pow( delta_acc/Delta, 0.3/ ( static_cast<LD>(RK_method::p + 1) ) );   
        // fac*=std::pow(delta_acc/Delta/Delta, 1/((LD) RK_method::p+1.) ) ;
        h_stop=true ;
        delta_acc=Delta;
        
    }else{
        fac*=std::pow( Delta , -1./(static_cast<LD>(RK_method::p +1.) ) );
        delta_rej=Delta;
    }
    
    //this is an alternative. Not very good for some reason. 
    // fac*=h/h_old*std::pow(delta_acc/Delta/Delta, 1/((LD) RK_method::p+1.) ) ;
    
    if(fac> fac_max){fac = fac_max;}
    if(fac< fac_min){fac = fac_min;}
    h_trial=h_trial*fac ;

    if (h_trial>hmax ){ h_trial=hmax; h_stop=true;}
    if (h_trial<hmin ){ h_trial=hmin; h_stop=true;}
}
/*-----------------------End: step_control---------------------------------*/

}

#endif