#ifndef RKF_steps
#define RKF_steps
#include "RKF_class.hpp"



/*---------------------------------------------------Begin: Get next step-------------------------------------------------------------------------------*/
_RKF_template_
_RKF_Func_::next_step(){

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
                ynext[eq] =  tmp_sol[eq] + bk[eq];
                //  std::cout<<bk[eq] << "  ";
                
                ynext_star[eq] =  tmp_sol[eq] + bstark[eq];       
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
_RKF_template_
_RKF_Func_::solve(){
        int tmp_step=1;
        
        int _hist_steps=0; 
        while (true )
        {
            //increase current_step
            current_step++;

            if( tn>=1.  or current_step == max_N  ) {break ;}
            
            next_step();
            _hist_steps++;

            
            for (int eq = 0; eq < N_eqs; eq++){tmp_sol[eq]=ynext[eq];}
            tn+= h0;

            if (  tn >=( (double ) tmp_step)/( (double )N_out-1.)   ){
                
                time[tmp_step]=tn;
                hist[tmp_step]=_hist_steps;

                
                for (int eq = 0; eq < N_eqs; eq++){
                    solution[eq][tmp_step] =  tmp_sol[eq] + bk[eq];
                    error[eq][tmp_step]= bk[eq] - bstark[eq];

                    }                

                tmp_step++;
                _hist_steps=0;

            }

                    

            
        }

    } 
/*---------------------------------------------------End: solve-------------------------------------------------------------------------------*/




#endif