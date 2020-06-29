#ifndef Ros_step_control
#define Ros_step_control
#include "Ros_class.hpp"

#define max(a,b)  (a <= b) ? b : a


/*-----------------------Begin: step_control---------------------------------*/
Ros_Template
void Ros_Namespace::step_control(){
    LD Delta=0.;
    LD _sc;
    LD fac=beta;
    
    for (int eq = 0; eq < N_eqs; eq++){
        _sc=max(std::abs( ynext[eq] ), std::abs( ynext_star[eq] ));
        _sc=abs_tol+rel_tol*_sc;
        Delta+= std::pow((abs_delta[eq]/_sc),2.);
        
    ;}
    Delta=std::sqrt(1./N_eqs*Delta);
    if(Delta<1) { h_stop=true ; }
    //step size cotrol from "Solving Ordinary Differential Equations I"
    fac*=std::pow( 1/Delta , 1./( (LD)method.p ) );
    if(fac > fac_max){fac=fac_max;}
    
    h0= beta*h0*fac ;

    if (h0>hmax ){ h0=hmax; h_stop=true ; }
    if (h0<hmin ){ 
        h0=hmin; h_stop=true;
        // std::cout<<"#minimum stepsize reached. Try increasing it to improve accuracy!\n"; 
    }


    if(h_stop){Deltas.push_back( Delta);}
    if (tn+h0>tmax ){ h0=tmax-tn;  }
    
    

}
/*-----------------------End: step_control---------------------------------*/


#endif