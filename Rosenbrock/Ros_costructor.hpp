#ifndef Ros_constructor
#define Ros_constructor
#include "Ros_class.hpp"



//The constructor. Remember that N has default value
_Ros_template_
_Ros_Cosnt_:: Ros(diffeq dydt, LD (&init_cond)[N_eqs] , 
        LD initial_step_size, LD minimum_step_size, LD maximum_step_size,int maximum_No_steps, 
        LD absolute_tolerance,LD relative_tolerance,LD beta,LD fac_max){
        // Initialize inputs
        this->dydt=dydt;
        this->Jac=jacobian(dydt);

        this->h0=initial_step_size;
        this->h1=initial_step_size;

        this->hmin=minimum_step_size;
        this->hmax=maximum_step_size;
        this->max_N=maximum_No_steps;
        this->abs_tol=absolute_tolerance;
        this->rel_tol=relative_tolerance;
        this->beta=beta;
        this->fac_max=fac_max;


        // ---------------------------------------------------------------------------------- //
       //define tmp_sol[N_eqs]. It is also good to initialize ynext.
        for(int i = 0; i < N_eqs ;++i) {
                this->tmp_sol[i]=init_cond[i];
                this->ynext[i]=init_cond[i];
                (this->solution)[i].push_back( init_cond[i]);
                (this->error)[i].push_back(0);
                
        }
        (this->time).push_back(0);
        (this->hist).push_back(0);

        // ---------------------------------------------------------------------------------- //
        
        // define k[N_eqs][method.s]. Initialize also k=0.
        this->k=new LD*[N_eqs];
        for(int i = 0; i < N_eqs ;++i) {
                this->k[i] = new LD[ this->method.s];
                for(int j =0 ; j<(this->method.s)-1; j++ ){this->k[i][j] =0;  }
                } 
        

        // calculate sums over gamma for all stages 
        this->sum_gamma=new LD[this->method.s];
        for(int stage = 0; stage < this->method.s; stage++){this->sum_gamma[stage]=0;
          for(int j =0 ; j<(this->method.s)-1; j++ ) {  this->sum_gamma[stage]+=this->method.g[stage][j]; }
          }


        //initialize tn, current_step, and End
        this->tn=0;
        this->current_step=0;

        //initialize this to 1 for the PI step control
        this->Deltas.push_back(1);
        
        };

//The destructor
_Ros_template_
_Ros_Cosnt_::~Ros(){
        // std::cout << "I'm done" << std::endl;
        delete[] this->k;

    };



#endif