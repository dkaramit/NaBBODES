#ifndef RKF_step_control
#define RKF_step_control
#include "RKF_class.hpp"

#define max(a,b)  (a <= b) ? b : a
#define min(a,b)  (a >= b) ? b : a


/*-----------------------Begin: step_control---------------------------------*/

RKF_Template
void RKF_Namespace::step_control(){
    
    //calculate the absolute value of delta
    


    LD Delta=0.;
    LD _sc;
    LD fac;
    
    for (int eq = 0; eq < N_eqs; eq++){
        _sc=max(fabs( ynext[eq] ), fabs( ynext_star[eq] ));
        _sc=abs_tol+rel_tol*_sc;
        Delta+= pow((abs_delta[eq]/_sc),2.);
        
    ;}
    Delta=std::sqrt(1./N_eqs*Delta);
    if(Delta==0){Delta=abs_tol;}

    //step size cotrol from "Solving Ordinary Differential Equations I"
    fac=min(std::pow( 1/Delta , 1./(method.p) ) , fac_max );

    h0= beta*h0*fac ;

    // std::cout<<h0<<"\t"<<fac<<"\t"<<Delta<<"\t"<<tn<<std::endl;

    if(Delta<1) { h_stop=true;}
    if (h0>hmax ){ h0=hmax;  }
    if (h0<hmin ){ 
        h0=hmin; h_stop=true ;
        // std::cout<<"#minimum stepsize reached. Try increasing it to improve accuracy!\n"; 
    }

    if(h_stop){ Deltas.push_back( Delta );}
    if (tn+h0>tmax ){ h0=tmax-tn;  }
    
    

}
/*-----------------------End: step_control---------------------------------*/


#endif