#ifndef Ros_class
#define Ros_class

#include<array>
#include<vector>
#include<functional>
#include<optional>

#include"Jacobian/Jacobian.hpp"

namespace rosenbrock{
    
// struct for the parameters of the RKF algorithm. 
// It is useful because we can now pass them named parameters!
template<class LD>
struct parameters{
    std::optional<LD> initial_step_size{}; 
    std::optional<LD> minimum_step_size{}; 
    std::optional<LD> maximum_step_size{};
    std::optional<unsigned int> maximum_No_steps{}; 
    std::optional<LD> absolute_tolerance{};
    std::optional<LD> relative_tolerance{};
    std::optional<LD> beta{};
    std::optional<LD> fac_max{}; 
    std::optional<LD> fac_min{};
};

//these are the default set of parameters
template<class LD>
inline constexpr parameters<LD> default_parameters {
    .initial_step_size=1e-5, 
    .minimum_step_size=1e-9, 
    .maximum_step_size=1e-3,
    .maximum_No_steps=1000000, 
    .absolute_tolerance=1e-6,
    .relative_tolerance=1e-6,
    .beta=0.85,
    .fac_max=3, 
    .fac_min=0.3
};

// these will be passed as template arguments to chose step controller
enum class  step_controllers{
    simple,
    PI
};

/*
LD it the numeric type to use.
RK_method is the method (the RODAS5 is the default)
step_controller chooses the step controller from the enum step_controllers (PI is the default)
*/

//This is a general implementation of embedded Rosenbrock solver of
// a system of differential equations in the interval [0,tmax].
template<class LD, class RK_method=RODAS5<LD>, step_controllers step_controller=step_controllers::PI> 
//Note that you can use template to pass the method
class Solver{
    public:
        // maybe it is useful to know the type of the equations
        using system_type=std::function<void(std::vector<LD> &lhs, const  std::vector<LD> &y, const LD &t)>;
        using Jacobian_type=std::function<void(std::vector<std::vector<LD>> &J, std::vector<LD> &dfdt, const std::vector<LD> &y, const LD& t)>;

    private:
        system_type dydt;
        Jacobian_type Jac;
        
        parameters<LD> params; //use this to get and change parameters if needed


        LD hmin, hmax, abs_tol, rel_tol, beta, fac_max, fac_min;
        unsigned int max_N;
        unsigned int N_eqs;
        LD h_old,h_trial,h_acc,delta_acc,delta_rej;//these will be initialized at the beginning of next_step
        bool h_stop;//h_stop becomes true when suitable stepsize is found.    
        
        LD tmax, tn;
        
        
        //these are here to hold the k's, sum_i b_i*k_i, sum_i b_i^{\star}*k_i, and sum_j a_{ij}*k_j 
        std::vector<std::vector<LD>> k;
        std::vector<LD> ak,gk,Jk, bk,bstark;
        // abs_delta=abs(ynext-ynext_star)
        std::vector<LD> abs_delta;
        
        /*--These are specific to Rosenbrock methods*/
        // need this to store the sum over \gammas (see the contructor)
        std::vector<LD> sum_gamma;
        // the t component of the jacobian
        std::vector<LD> dfdt; 
        //define the coefficient. This will become (I-\gamma*h*J). _inv is its inverse
        std::vector<std::vector<LD>> _inv;
        // There are for the LUP-decomposition of (I-\gamma*h*J) 
        std::vector<std::vector<LD>> L;
        std::vector<std::vector<LD>> U;
        std::vector<int> P;
        //lu_sol will capture the sulution of (I-\gamma*h*J)* k = rhs (i.e. k = (I-\gamma*h*J)^{-1} rhs)
        std::vector<LD> lu_sol;
        std::vector<std::vector<LD>> J;//this is here to hold values of the Jacobian
        
        std::vector<LD> yprev;// previously accepted step. maybe the name is not good.
        std::vector<LD> ynext;//this is here to hold the prediction
        std::vector<LD> ynext_star;//this is here to hold the second prediction
        
        std::vector<LD> time;
        std::vector<std::vector<LD>> solution;
        std::vector<std::vector<LD>> error;
        
        void LU();
        void calc_k();
        void calc_Jk();
        
        void sum_ak(const unsigned int& stage); // calculate sum_j a_{ij}*k_j and passit to this->ak
        void sum_gk(const unsigned int& stage); // calculate sum_j g_{ij}*k_j and passit to this->gk
        void sum_bk();// calculate sum_i b_i*k_i and passit to this->bk 
        
        void step_control_PI();//adjust stepsize until error is acceptable
        void step_control_simple();//adjust stepsize until error is acceptable
        
        void step_control(){
            if constexpr (step_controller==step_controllers::PI){
                step_control_PI();
            }else{
                step_control_simple();
            }
        };//adjust stepsize until error is acceptable

    public:

        // Notice that if you use the default Jacobian, you have the option to change its default value for h.
        Solver(const system_type& dydt, const std::vector<LD> &init_cond, LD tmax, const parameters<LD>& opt=default_parameters<LD>, const LD& Jacobian_h=1e-8)
            :dydt(dydt),Jac(Jacobian<LD>(dydt,Jacobian_h)), params(opt), N_eqs(init_cond.size()){reset(init_cond,tmax,opt);}
            
            
        Solver(const system_type& dydt, const std::vector<LD> &init_cond, LD tmax, Jacobian_type Jac, const parameters<LD>& opt=default_parameters<LD>)
            :dydt(dydt),Jac(Jac), params(opt), N_eqs(init_cond.size()){reset(init_cond,tmax,opt);}

        
        ~Solver()=default;

        const std::vector<LD>& get_t() const { return time; }
        auto get_t(const unsigned int& step) const { return time.at(step); }
        //access the array of solution[eq]
        const std::vector<LD>& get_solution(const unsigned int& eq) const { return solution.at(eq); }
        //access the element solution[eq][step]
        auto get_solution(const unsigned int& eq, const unsigned int& step) const { return solution.at(eq).at(step); }
        
        //access the array of error[eq]
        const std::vector<LD>& get_error(const unsigned int& eq) const { return error.at(eq); }
        //access the element error[eq][step]
        auto get_error(const unsigned int& eq, const unsigned int& step) const { return error.at(eq).at(step); }
        
        
        void next_step();
        void solve();

        void set_parameters(const parameters<LD>& opt=default_parameters<LD>);
        void reset(const std::vector<LD>& init_cond, LD tmax, const parameters<LD>& opt=default_parameters<LD>);

        //generally helpful, but not very important
        auto get_current_step() const {return time.size();}
        auto get_current_step_size() const {return h_acc;}
        const parameters<LD>& get_parameters() const {return params;}
};

}

#endif
