#ifndef Ros_constructor
#define Ros_constructor
#include "Ros_class.hpp"

template<unsigned int N_eqs, class RK_method, class jacobian, class LD> 
void Ros<N_eqs, RK_method,  jacobian, LD>::set_parameters(const parameters<LD>& opt){

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

template<unsigned int N_eqs, class RK_method, class jacobian, class LD> 
void Ros<N_eqs, RK_method,  jacobian, LD>::reset(const std::array<LD,N_eqs>& init_cond, LD tmax, const parameters<LD>& opt){
    
    this->tmax=tmax;
    set_parameters(opt);
    
    // ---------------------------------------------------------------------------------- //
    //define yprev[N_eqs]. It is also good to initialize ynext.
    (this->time).clear();
    (this->time).push_back(0);
    for(unsigned int i = 0; i < N_eqs ;++i) {
        this->yprev[i]=init_cond[i];
        this->ynext[i]=init_cond[i];
        
        (this->solution)[i].clear();
        (this->solution)[i].push_back( init_cond[i]);

        (this->error)[i].clear();
        (this->error)[i].push_back(0);   
    }

    // ---------------------------------------------------------------------------------- //

    //Initialize also k=0.
    for(unsigned int i = 0; i < N_eqs ;++i){for(unsigned int j =0 ; j<RK_method::s; j++ ){this->k[i][j] =0;}} 


    // calculate sums over gamma for all stages 
    for(unsigned int stage = 0; stage < RK_method::s; stage++){
        this->sum_gamma[stage]=0;
        for(unsigned int j =0 ; j<RK_method::s; j++ ){ this->sum_gamma[stage]+=RK_method::g[stage][j];}
    }

    //initialize tn
    this->tn=0;
    //initialize delta_acc
    delta_acc=1.;
}

#endif