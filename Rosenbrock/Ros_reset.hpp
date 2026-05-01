#ifndef Ros_constructor
#define Ros_constructor
#include "Ros_class.hpp"
#include<vector>

#define parameter_check(opt_name)  \
    if(!params.opt_name.has_value()){params.opt_name=default_parameters<LD>.opt_name.value();}

#define option_check(opt_name,class_name)  \
    if(opt.opt_name.has_value()) {params.opt_name=class_name=opt.opt_name.value();}else{class_name=params.opt_name.value();}

namespace rosenbrock{

template<class LD, class RK_method, step_controllers step_controller> 
void Solver<LD, RK_method, step_controller>::set_parameters(const parameters<LD>& opt){

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

template<class LD, class RK_method, step_controllers step_controller> 
void Solver<LD, RK_method, step_controller>::reset(const std::vector<LD>& init_cond, LD tmax, const parameters<LD>& opt){
    k.assign(N_eqs,std::vector<LD>(RK_method::s,0));
    ak.assign(N_eqs,0);
    gk.assign(N_eqs,0);
    Jk.assign(N_eqs,0);
    bk.assign(N_eqs,0);
    bstark.assign(N_eqs,0);
    abs_delta.assign(N_eqs,0);
    sum_gamma.assign(RK_method::s,0);
    dfdt.assign(N_eqs,0);
    _inv.assign(N_eqs,std::vector<LD>(N_eqs,0));
    L.assign(N_eqs,std::vector<LD>(N_eqs,0));
    U.assign(N_eqs,std::vector<LD>(N_eqs,0));
    P.assign(N_eqs,0);
    lu_sol.assign(N_eqs,0);
    J.assign(N_eqs,std::vector<LD>(N_eqs,0));
    yprev.assign(N_eqs,0);
    ynext.assign(N_eqs,0);
    ynext_star.assign(N_eqs,0);
    solution.assign(N_eqs,std::vector<LD>(0,0));
    error.assign(N_eqs,std::vector<LD>(0,0));

    if(N_eqs!=init_cond.size()){
        throw std::logic_error(
                "Size mismatch: number of equation used = " + std::to_string(N_eqs) +
                ", the new initial conditions have size = " + std::to_string(init_cond.size())
            );
    }

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