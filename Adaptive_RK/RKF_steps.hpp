#ifndef RKF_steps
#define RKF_steps
#include "RKF_class.hpp"



/*---------------------------------------------------Begin: Get next step-------------------------------------------------------------------------------*/
template<class diffeq, int N_eqs, class RKF_method>
void RKF<diffeq, N_eqs, RKF_method>::next_step(){

        //set h_stop=false, to start looking for stepsize
        h_stop=false;

        //calculate ynext and ynext_star until h_stop=true 
        while (true) 
        {
            // calculate \vec{k}:
            calc_k();
        
            // now we can calulate \vec{y}_{n+1}
            // for this we need sum_{i}^{s}b_{i}\vec{k}_i *h. That is, call sum_bk
            sum_bk();
            
            // having bk, we now have \vec{y}_{n+1} \vec{y}^{\star}_{n+1}. 
            for (int eq = 0; eq < N_eqs; eq++)
            {   
                ynext[eq] =  solution[eq][current_step-1] + bk[eq];
                //  std::cout<<bk[eq] << "  ";
                
                ynext_star[eq] =  solution[eq][current_step-1] + bstark[eq];       
                // std::cout<<ynext_star[eq]<< "  "<<ynext[eq]<< "  "<<ynext_star[eq]-ynext[eq]<<"\n";
                
                // calculate the error^2
                abs_delta[eq]= ynext[eq] - ynext_star[eq] ;
                
            }

            // call step_control to see if the error is acceptable
            step_control();
            if(h_stop){break;}

        }
            
}
/*---------------------------------------------------End: Get next step-------------------------------------------------------------------------------*/



/*---------------------------------------------------Begin: solve-------------------------------------------------------------------------------*/

template<class diffeq, int N_eqs, class RKF_method>
void RKF<diffeq, N_eqs, RKF_method>::solve(){

        while (true )
        {
            //increase current_step
            current_step++;

            if( tn>=1.  or current_step == max_N  ) {break ;}
            
            next_step();

            
            for (int eq = 0; eq < N_eqs; eq++){solution[eq][current_step]=ynext[eq];}
            tn+= h0;
            // std::cout<< tn <<"  "<< current_step <<"\n";
                    
            steps[current_step] = tn;

            
        }

    } 
/*---------------------------------------------------End: solve-------------------------------------------------------------------------------*/




#endif