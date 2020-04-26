#ifndef Ros_steps
#define Ros_steps
#include "Ros_class.hpp"



/*---------------------------------------------------Begin: Get next step-------------------------------------------------------------------------------*/
_Ros_template_
_Ros_Func_::next_step(){

        //set h_stop=false, to start looking for stepsize
        h_stop=false;

        //calculate the LU decomposition of (1-\gamma*h*J) and find its inverse before you enter the loop. 
        LU();

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
            
            if(h_stop){break; }

        }
            
    

        // std::cout<<current_step<<std::endl;
}
/*---------------------------------------------------End: Get next step-------------------------------------------------------------------------------*/



/*---------------------------------------------------Begin: solve-------------------------------------------------------------------------------*/

_Ros_template_
_Ros_Func_::solve(bool _full_){

        if( _full_ )
        {
            for (int eq = 0; eq < N_eqs; eq++){ solution_full[eq].push_back( tmp_sol[eq] ); }
            time_full.push_back(tn);
        }
    

        int tmp_step=1;
        
        // Use this to count how many steps you take between entries in time and solution. 
        // This basically makes a histogram os number of steps, which is stored in the vector hist.
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

            if (  tn >=( (LD ) tmp_step)/( (LD )N_out-1.)   ){
                
                time.push_back(tn);
                hist.push_back(_hist_steps);

                
                for (int eq = 0; eq < N_eqs; eq++)
                {
                    solution[eq].push_back(ynext[eq]);
                    error[eq].push_back(ynext[eq] - ynext_star[eq]);
                }                

                tmp_step++;
                _hist_steps=0;

            }

            if( _full_ )
            {
                for (int eq = 0; eq < N_eqs; eq++){ solution_full[eq].push_back( ynext[eq] ); }
                time_full.push_back(tn);
            }

            
        }


    } 
/*---------------------------------------------------End: solve-------------------------------------------------------------------------------*/




#endif