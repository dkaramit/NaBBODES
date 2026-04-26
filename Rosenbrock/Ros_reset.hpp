#ifndef Ros_constructor
#define Ros_constructor
#include "Ros_class.hpp"

#define parameter_check(opt_name)  \
    if(!params.opt_name.has_value()){params.opt_name=default_parameters<LD>.opt_name.value();}

#define option_check(opt_name,class_name)  \
    if(opt.opt_name.has_value()) {params.opt_name=class_name=opt.opt_name.value();}else{class_name=params.opt_name.value();}

namespace Rosenbrock{

template<unsigned int N_eqs, class RK_method, class LD, step_controlers step_controler> 
void Solver<N_eqs, RK_method, LD, step_controler>::set_parameters(const parameters<LD>& opt){

    // if some parameter in opt does not have a value, use the corresponding parameter from default.default_parameters 
    parameter_check(initial_step_size)
    parameter_check(minimum_step_size)
    parameter_check(maximum_step_size)
    parameter_check(maximum_No_steps)
    parameter_check(absolute_tolerance)
    parameter_check(relative_tolerance)
    parameter_check(beta)
    parameter_check(fac_max)
    parameter_check(fac_min)

    // change the parameters that exist in opt (inclusing an update to the corresponding parameter in params). 
    // If a parameter does not have a value. set it to the existing params.
    option_check(initial_step_size,h_trial)
    option_check(minimum_step_size,hmin)
    option_check(maximum_step_size,hmax)
    option_check(maximum_No_steps,max_N)
    option_check(absolute_tolerance,abs_tol)
    option_check(relative_tolerance,rel_tol)
    option_check(beta,beta)
    option_check(fac_max,fac_max)
    option_check(fac_min,fac_min)
    
    h_acc=h_trial;
}

template<unsigned int N_eqs, class RK_method, class LD, step_controlers step_controler> 
void Solver<N_eqs, RK_method, LD, step_controler>::reset(const std::array<LD,N_eqs>& init_cond, LD tmax, const parameters<LD>& opt){
    
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

}

#undef parameter_check
#undef option_check

#endif