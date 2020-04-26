#ifndef Ros_step_control
#define Ros_step_control
#include "Ros_class.hpp"

#define max(a,b)  (a <= b) ? b : a
#define min(a,b)  (a >= b) ? b : a


/*-----------------------Begin: step_control---------------------------------*/
_Ros_template_
_Ros_Func_::step_control(){
    
    //calculate the absolute value of delta
    
    LD Delta=0.;
    LD _sc;
    LD fac;
    

    dot<N_eqs,LD>( _inv,abs_delta , regulated_delta);
    for (int eq = 0; eq < N_eqs; eq++){
        _sc=max(fabs( ynext[eq] ), fabs( ynext_star[eq] ));
        _sc=abs_tol+rel_tol*_sc;
        Delta+= pow((regulated_delta[eq]/_sc),2.);
    }

    Delta=pow(1./N_eqs*Delta,0.5);
    if(Delta==0){Delta=abs_tol;}
    
    
    //PI step size cotrol from "Solving Ordinary Differential Equations II"
    fac=  h0/h1*pow(Deltas[current_step-1]/Delta/Delta, 1/((LD) method.p+1.) ) ;
    
    //this step controller depends on the next to the previous accepted step. 
    // So, update it in order to use it for the next step.
    if(Delta<1) { h_stop=true ; Deltas.push_back( Delta);  h1=h0; } 
    
    if(fac> fac_max){fac = fac_max;}

    h0= beta*h0*fac ;

    if (h0>hmax ){ h0=hmax;  }
    if (h0<hmin ){ h1=h0; 
        h0=hmin; h_stop=true ;  Deltas.push_back( Delta) ; 
        // std::cout<<"#minimum stepsize reached. Try increasing it to improve accuracy!\n"; 
    }

    if (tn+h0>1. ){ h0=1-tn;  }
    
    

}
/*-----------------------End: step_control---------------------------------*/


#endif