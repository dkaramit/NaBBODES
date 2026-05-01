#ifndef RKF_class
#define RKF_class

#include<vector>
#include<functional>
#include<optional>

#include"METHOD.hpp"

namespace rkf{
    
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

//This is a general implementation of explicit embedded RK solver of
// a system of differential equations in the interval [0,tmax].

/*
LD it the numeric type to use.
RK_method is the method (the DormandPrince is the default)
step_controller chooses the step controller from the enum step_controllers (PI is the default)
*/

//This is a general implementation of explicit embedded RK solver of
// a system of differential equations in the interval [0,tmax].
template<class LD, class RK_method=DormandPrince<LD>, step_controllers step_controller=step_controllers::PI>
class Solver{
    public:
    // maybe it is useful to know the type of the equation
    // system_type is a collable of the system of  equations to be solved 
    using system_type=std::function<void(std::vector<LD> &lhs, const  std::vector<LD> &y, const LD &t)>;

    private:
        parameters<LD> params; //use this to get and change parameters if needed

        //these are not constant because reset can update them
        LD hmin, hmax, abs_tol, rel_tol, beta, fac_max, fac_min;

        LD h_old,h_trial,h_acc,delta_acc, delta_rej;//these will be initialized at the beginning of next_step
        unsigned int max_N;
        bool h_stop;//h_stop becomes true when suitable stepsize is found.    
        
        //Inputs. The initial condition is given as a Array (the type is users choice as long as it can be called with [])
        const system_type dydt;
        // integrate up to t=tmax
        LD tmax;
        
        LD tn;

        unsigned int N_eqs;

        std::vector<std::vector<LD>> k;
        
        //these are here to hold the k's, sum_i b_i*k_i, sum_i b_i^{\star}*k_i, and sum_j a_{ij}*k_j 
        std::vector<LD> ak, bk, bstark;
        // abs_delta=abs(ynext-ynext_star)
        std::vector<LD> abs_delta;
        
        std::vector<LD> yprev;
        std::vector<LD> ynext;//this is here to hold the prediction
        std::vector<LD> ynext_star;//this is here to hold the second prediction
        
        std::vector<LD> time;
        std::vector<std::vector<LD>> solution;
        std::vector<std::vector<LD>> error;

        void calc_k(); //calculate the values of k
        void sum_ak(const unsigned int& stage); // calculate sum_j a_{ij}*k_j and pass it to this->ak
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
        
        Solver(const system_type&  dydt, const std::vector<LD>& init_cond, const LD& tmax, const parameters<LD>& opt=default_parameters<LD>): 
            dydt(dydt),params(opt),N_eqs(init_cond.size()){reset(init_cond,tmax,opt); };
        
        ~Solver()=default;

        void next_step();
        void solve();

        
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
        
        
        void set_parameters(const parameters<LD>& opt=default_parameters<LD>);
        void reset(const std::vector<LD>& init_cond, const LD& tmax, const parameters<LD>& opt=default_parameters<LD>);
        
        //generally helpful, but not very important
        auto get_current_step() const {return time.size();}
        auto get_current_step_size() const {return h_acc;}
        const parameters<LD>& get_parameters() const {return params;}
    };

}

#endif
