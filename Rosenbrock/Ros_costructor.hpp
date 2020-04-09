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
        this->hmin=minimum_step_size;
        this->hmax=maximum_step_size;
        this->max_N=maximum_No_steps;
        this->abs_tol=absolute_tolerance;
        this->rel_tol=relative_tolerance;
        this->beta=beta;
        this->fac_max=fac_max;

        // ---------------------------------------------------------------------------------- //
       //define tmp_sol[N_eqs]. It is also good to initialize ynext.
        this->tmp_sol = new LD[N_eqs]; 
        for(int i = 0; i < N_eqs ;++i) {
                this->tmp_sol[i]=init_cond[i];
                this->ynext[i]=init_cond[i];
        }

        this->hist = new int[N_out];//make a list in which you'll put the steps it took between time[i] and time[i+1] in order to make a histogram
        this->time = new LD[N_out];//make a list in which you'll put the steps (these will be approximately at intervals of 1/(N_out-1))
        this->time[0]=0;
        this->hist[0]=0;


        this->solution = new LD*[N_eqs];
        this->error = new LD*[N_eqs];
        for(int i = 0; i < N_eqs ;++i) {
                this->solution[i] = new LD[ N_out ];
                this->error[i] = new LD[ N_out ];
                this->error[i][0]=0;
                this->solution[i][0]=init_cond[i];
            } 
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
        delete[] this->solution;
        delete this->time;
        delete this->error;
        delete this->hist;

        delete this->tmp_sol;
        delete[] this->k;

    };



#endif