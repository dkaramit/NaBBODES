#ifndef RKF_constructor
#define RKF_constructor
#include "RKF_class.hpp"

// change the parameters of the solver
template<unsigned int N_eqs, class RK_method, class LD>
void RKF<N_eqs, RK_method, LD>::set_parameters(const parameters<LD>& opt){

    // change the parameters that exist in opt (inclusing an update to the corresponding parameter in params). 
    // If a parameter does not have a value. set it to the existing params.
    if(opt.initial_step_size.has_value()) {params.initial_step_size=h_trial=opt.initial_step_size.value();}else{h_trial=params.initial_step_size.value();}
    if(opt.minimum_step_size.has_value()) {params.minimum_step_size.value()=hmin=opt.minimum_step_size.value();}else{hmin=params.minimum_step_size.value();}
    if(opt.maximum_step_size.has_value()) {params.maximum_step_size.value()=hmax=opt.maximum_step_size.value();}else{hmax=params.maximum_step_size.value();}
    if(opt.maximum_No_steps.has_value())  {params.maximum_No_steps=max_N=opt.maximum_No_steps.value();}else{max_N=params.maximum_No_steps.value();}
    if(opt.absolute_tolerance.has_value()){params.absolute_tolerance=abs_tol=opt.absolute_tolerance.value();}else{abs_tol=params.absolute_tolerance.value();}
    if(opt.relative_tolerance.has_value()){params.relative_tolerance=rel_tol=opt.relative_tolerance.value();}else{rel_tol=params.relative_tolerance.value();}
    if(opt.beta.has_value())              {params.beta=beta=opt.beta.value();}else{beta=params.beta.value();}
    if(opt.fac_max.has_value())           {params.fac_max=fac_max=opt.fac_max.value();}else{fac_max=params.fac_max.value();}
    if(opt.fac_min.has_value())           {params.fac_min=fac_min=opt.fac_min.value();}else{fac_min=params.fac_min.value();}
    
    h_acc=h_trial;
}

// reset the solver (set clear arrays, set t=0, set y to its initial condition, etc.)
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