#ifndef RKF_step_control
#define RKF_step_control
#include "RKF_class.hpp"

#define max(a,b)  (a <= b) ? b : a


/*-----------------------Begin: step_control---------------------------------*/

RKF_Template
void RKF_Namespace::step_control(){
    LD Delta=0.;
    LD _sc=0;
    LD fac=beta;
    
    for (int eq = 0; eq < N_eqs; eq++){
        _sc=max(std::abs( ynext[eq] ), std::abs( ynext_star[eq] ));
        _sc=abs_tol+rel_tol*_sc;
        Delta+= std::pow((abs_delta[eq]/_sc),2.);
        
    ;}
    Delta=std::sqrt(1./N_eqs*Delta);
    if(Delta==0){Delta=abs_tol*rel_tol;}

    // I use the step control from 
    // https://www.sciencedirect.com/science/article/pii/S147466701751767X
    if(Delta<=1) { 
        if(h_stop==false){fac*=h0/h1;}
        fac*=std::pow(Delta, -0.65/( (LD) method.p + 1.) );   
        fac*=std::pow( Deltas[current_step-1]/Delta, 0.3/ ( (LD) method.p + 1. ) );   
        h_stop=true ;
    }else{
        fac*=std::pow( Delta , -1./((LD) method.p + 1.) );
    }
    
    if(fac>fac_max){fac=fac_max;}

    h0= h0*fac ;
    // std::cout<<h0<<"\t"<<fac<<"\t"<<Delta<<"\t"<<tn<<std::endl;

    
    if (h0>hmax ){ h0=hmax; h_stop=true;}
    
    if (h0<hmin ){ 
        h0=hmin; h_stop=true ;
        // std::cout<<"#minimum stepsize reached. Try increasing it to improve accuracy!\n"; 
    }

    if(h_stop){ Deltas.push_back( Delta ); }
    if (tn+h0>tmax ){ h0=tmax-tn;  }
}
/*-----------------------End: step_control---------------------------------*/


#endif