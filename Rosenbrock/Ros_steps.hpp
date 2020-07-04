#ifndef Ros_steps
#define Ros_steps
#include "Ros_class.hpp"



/*---------------------------------------------------Begin: Get next step-------------------------------------------------------------------------------*/
Ros_Template
void Ros_Namespace::next_step(){
    //set h_stop=false, to start looking for stepsize
    h_stop=false;
    h_old=h;//for the PI controller
    delta_acc=Deltas.back();//for the PI controller
    delta_rej=delta_acc;
    //calculate the LU decomposition of (1-\gamma*h*J) and find its inverse before you enter the loop. 
    LU();
    //calculate ynext and ynext_star until h_stop=true 
    while (true) {
        // calculate \vec{k}:
        calc_k();
        // now we can calulate \vec{y}_{n+1}
        // for this we need sum_{i}^{s}b_{i}\vec{k}_i *h. That is, call sum_bk
        sum_bk();
        // having bk, we now have \vec{y}_{n+1} \vec{y}^{\star}_{n+1}. 
        for (int eq = 0; eq < N_eqs; eq++){   
            ynext[eq] =  tmp_sol[eq] + bk[eq];
            ynext_star[eq] =  tmp_sol[eq] + bstark[eq];       
            // calculate the error
            abs_delta[eq]= ynext[eq] - ynext_star[eq] ;
        }
        // call step_control to see if the error is acceptable
        step_control();
        if(h_stop){break; }
    }
}
/*---------------------------------------------------End: Get next step-------------------------------------------------------------------------------*/



/*---------------------------------------------------Begin: solve-------------------------------------------------------------------------------*/

Ros_Template
void Ros_Namespace::solve(bool _full_){
    if( _full_ ){
        //the initial values are not set in the contructor because one may not want the full solution
        for (int eq = 0; eq < N_eqs; eq++){ 
            solution_full[eq].push_back( tmp_sol[eq] ); 
            error_full[eq].push_back(0);
            }
        time_full.push_back(0);//tn=0 here
    }

    int tmp_step=1;
    // Use this to count how many steps you take between entries in time and solution. 
    // This basically makes a histogram os number of steps, which is stored in the vector hist.
    int _hist_steps=0; 
    while (true ){
        //increase current_step
        current_step++;
        
        if( tn>=tmax  or current_step == max_N  ) {break ;}
        
        next_step();
        
        //set tmp_sol = ynext, since ynext changes until h_stop=false in the next step
        for (int eq = 0; eq < N_eqs; eq++){tmp_sol[eq]=ynext[eq];}
        
        tn+= h;
        
        _hist_steps++;
        if (  tn >=( (LD) tmp_step)/( (LD)N_out-1.)*tmax   ){
            time.push_back(tn);
            hist.push_back(_hist_steps);
            for (int eq = 0; eq < N_eqs; eq++){
                solution[eq].push_back(ynext[eq]);
                error[eq].push_back(ynext[eq] - ynext_star[eq]);
            }                
            tmp_step++;
            _hist_steps=0;
        }

        if( _full_ ){
            for (int eq = 0; eq < N_eqs; eq++){ 
                solution_full[eq].push_back( ynext[eq] ); 
                error_full[eq].push_back(ynext[eq] - ynext_star[eq]);
            }
            time_full.push_back(tn);
        }
    }
} 
/*---------------------------------------------------End: solve-------------------------------------------------------------------------------*/




#endif
