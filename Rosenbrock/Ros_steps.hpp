#ifndef Ros_steps
#define Ros_steps
#include "Ros_class.hpp"



/*---------------------------------------------------Begin: Get next step-------------------------------------------------------------------------------*/
template<unsigned int N_eqs, class RK_method, class jacobian, class LD> 
void Ros<N_eqs, RK_method,  jacobian, LD>::next_step(){
    if (tn + h_trial > tmax){h_trial = tmax - tn;} //make sure that h_trial does not push tn>tmax. 
    //set h_stop=false, to start looking for stepsize
    h_stop=false;
    
    //for the PI controller: h has to be the laste *accepted* stepsize. Remember that h_trial is the last updaed h, once the last one is accepted. So, the last accepted step size is h, not h_trial. 
    h_old=h; 
    
    //for the PI controller: initialize  delta_rej to delta_acc.
    delta_rej=delta_acc;
    
    // Jacobian does not depend on h, so we can calculate it once at the beginning of the step control loop.
    Jac(J,dfdt,yprev,tn);
    
    
    //calculate ynext and ynext_star until h_stop=true 
    while (true){
        
        //calculate the LU decomposition of (1-\gamma*h*J) and find its inverse. It depends on h, so it needs to be updated as we try to find suitable h. 
        LU();

        // calculate \vec{k}:
        calc_k();
        
        // now we can calulate \vec{y}_{n+1}
        // for this we need sum_{i}^{s}b_{i}\vec{k}_i *h. That is, call sum_bk
        sum_bk();
        
        // having bk, we now have \vec{y}_{n+1} \vec{y}^{\star}_{n+1}. 
        for (unsigned int eq = 0; eq < N_eqs; eq++){   
            ynext[eq] =  yprev[eq] + bk[eq];
            ynext_star[eq] =  yprev[eq] + bstark[eq];       
            // calculate the error
            abs_delta[eq]= ynext[eq] - ynext_star[eq] ;
        }
        
        h=h_trial;//note that h is the trial step size, while h is the accepted one (the one used to compute y_next).

        // call step_control to see if the error is acceptable
        step_control();
            
        if(h_stop){break;}
             
    }
}
/*---------------------------------------------------End: Get next step-------------------------------------------------------------------------------*/



/*---------------------------------------------------Begin: solve-------------------------------------------------------------------------------*/

template<unsigned int N_eqs, class RK_method, class jacobian, class LD> 
void Ros<N_eqs, RK_method,  jacobian, LD>::solve(){
    unsigned int current_step=0;
    while (true){
        //increase current_step
        current_step++;
        // check if it's done
        if(tn>=tmax  or current_step == max_N){break;}
        // run next step
        next_step();
        
        // set previous y to last one
        for (unsigned int eq = 0; eq < N_eqs; eq++){yprev[eq]=ynext[eq];}
        // increase time by the accepted step size (not the trial one, which is h)
        tn+=h;

        // store everything
        time.push_back(tn);
        for (unsigned int eq = 0; eq < N_eqs; eq++){ 
            solution[eq].push_back( ynext[eq] ); 
            error[eq].push_back(ynext[eq] - ynext_star[eq]);
        }
    }
} 
/*---------------------------------------------------End: solve-------------------------------------------------------------------------------*/




#endif
