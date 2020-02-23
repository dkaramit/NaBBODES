#ifndef RKF_step_control
#define RKF_step_control
#include "RKF_class.hpp"

#define max(a,b)  (a <= b) ? b : a
#define min(a,b)  (a >= b) ? b : a


/*-----------------------Begin: step_control---------------------------------*/
_RKF_template_
_RKF_Func_::step_control(){
    
    //calculate the absolute value of delta
    


    LD Delta=0.;
    LD _sc;
    LD fac;
    
    for (int eq = 0; eq < N_eqs; eq++){
        _sc=max(fabs( ynext[eq] ), fabs( ynext_star[eq] ));
        _sc=abs_tol+rel_tol*_sc;
        Delta+= pow((abs_delta[eq]/_sc),2.);
        
    ;}
    Delta=pow(1./N_eqs*Delta,0.5);


    if(Delta<1) { h_stop=true ; Deltas.push_back( Delta) ;}
    //step size cotrol from "Solving Ordinary Differential Equations I"
    fac=min(pow( 1/Delta , 1./(method.p+1) ) , fac_max );
    
    h0= beta*h0*fac ;

    if (h0>hmax ){ h0=hmax;  }
    if (h0<hmin ){ 
        h0=hmin; h_stop=true ;  Deltas.push_back( Delta) ;
        std::cout<<"#minimum stepsize reached. Try increasing it to improve accuracy!\n"; 
    }

    if (tn+h0>1. ){ h0=1-tn;  }
    
    

}
/*-----------------------End: step_control---------------------------------*/


#endif