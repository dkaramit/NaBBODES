#ifndef RKF_constructor
#define RKF_constructor
#include "RKF_class.hpp"



//The constructor. Remember that N has default value
_RKF_template_
_RKF_Cosnt_:: RKF(diffeq dydt, LD (&init_cond)[N_eqs] ,
        LD initial_step_size, LD minimum_step_size, LD maximum_step_size,int maximum_No_steps, 
        LD absolute_tolerance,LD relative_tolerance,LD beta,LD fac_max){
        // Initialize inputs
        this->dydt=dydt;
        this->h0=initial_step_size;
        this->hmin=minimum_step_size;
        this->hmax=maximum_step_size;
        this->max_N=maximum_No_steps;
        this->abs_tol=absolute_tolerance;
        this->rel_tol=relative_tolerance;
        this->beta=beta;
        this->fac_max=fac_max;

        // ---------------------------------------------------------------------------------- //
        // later, I'll make steps and tmp_sol std::vector
        
        //define tmp_sol[N_eqs]
        this->tmp_sol = new LD[N_eqs]; 
        for(int i = 0; i < N_eqs ;++i) {this->tmp_sol[i]=init_cond[i];}

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
        
        // define k[N_eqs][method.s]
        this->k=new LD*[N_eqs];
        for(int i = 0; i < N_eqs ;++i) {
                this->k[i] = new LD[ this->method.s];
                } 
        
        
        //initialize tn, current_step, and End
        this->tn=0;
        this->current_step=0;
        
        };

//The destructor
_RKF_template_
_RKF_Cosnt_::~RKF(){
        // std::cout << "I'm done" << std::endl;
        delete[] this->solution;
        delete this->time;
        delete this->error;
        delete this->hist;

        delete this->tmp_sol;
        delete[] this->k;
    };



#endif