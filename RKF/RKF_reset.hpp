#ifndef RKF_constructor
#define RKF_constructor
#include "RKF_class.hpp"


template<unsigned int N_eqs, class RK_method, class LD>
void RKF<N_eqs, RK_method, LD>::set_parameters(const parameters<LD>& opt){

    this->h_trial=opt.initial_step_size;
    this->h_acc=this->h_trial;
    
    this->hmin=opt.minimum_step_size;
    this->hmax=opt.maximum_step_size;
    this->max_N=opt.maximum_No_steps;
    this->abs_tol=opt.absolute_tolerance;
    this->rel_tol=opt.relative_tolerance;
    this->beta=opt.beta;
    this->fac_max=opt.fac_max;
    this->fac_min=opt.fac_min;
    
}


template<unsigned int N_eqs, class RK_method, class LD>
void RKF<N_eqs, RK_method, LD>::reset(const std::array<LD,N_eqs>& init_cond, const LD& tmax, const parameters<LD>& opt){
    
    set_parameters(opt);
    
    this->tmax=tmax;
    // ---------------------------------------------------------------------------------- //
    this->time.clear();
    this->time.push_back(0);
    //define yprev[N_eqs]. It is also good to initialize ynext.
    for(unsigned int i = 0; i < N_eqs ;++i) {
        this->yprev[i]=init_cond[i];
        this->ynext[i]=init_cond[i];

        (this->solution)[i].clear();
        (this->solution)[i].push_back(init_cond[i]);
        
        (this->error)[i].clear();
        (this->error)[i].push_back(0);
    }

    // ---------------------------------------------------------------------------------- //

    // Initialize k=0 for definiteness.
    for(unsigned int i = 0; i < N_eqs ;++i){for(unsigned int j=0; j<RK_method::s; j++ ){ this->k[i][j]=0;}} 

    //initialize tn, current_step, and End
    this->tn=0;
    //initialize delta_acc
    delta_acc=1.;
}


#endif