#ifndef Ros_class
#define Ros_class

#include<array>
#include<vector>
#include<functional>
#include<optional>


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

/*
N_eqs is ten number of equations to be solved RK_method 
is the method
*/

//This is a general implementation of explicit embedded RK solver of
// a system of differential equations in the interval [0,tmax].
template<unsigned int N_eqs, class RK_method, class jacobian, class LD> 
//Note that you can use template to pass the method
class Ros{
    private:
        using diffeq=std::function<void(std::array<LD, N_eqs> &lhs, const  std::array<LD, N_eqs> &y, const LD &t)>;
        diffeq dydt;
        jacobian Jac;

        parameters<LD> params; //use this to get and change parameters if needed


        LD hmin, hmax, abs_tol, rel_tol, beta, fac_max, fac_min;
        unsigned int max_N;
        LD h_old,h_trial,h_acc,delta_acc,delta_rej;//these will be initialized at the beginning of next_step
        bool h_stop;//h_stop becomes true when suitable stepsize is found.    
        
        LD tmax, tn;
        std::array<LD, N_eqs> yprev;// previously accepted step. maybe the name is not good.
    
    
        std::vector<LD> time;
        std::array<std::vector<LD>, N_eqs> solution;
        std::array<std::vector<LD>, N_eqs> error;
        
        
        //these are here to hold the k's, sum_i b_i*k_i, sum_i b_i^{\star}*k_i, and sum_j a_{ij}*k_j 
        std::array<std::array<LD,RK_method::s>,N_eqs> k;
        std::array<LD,N_eqs> ak,gk,Jk, bk,bstark;
        // need this to store the sum over \gammas (see the contructor)
        std::array<LD,RK_method::s> sum_gamma;
        // abs_delta=abs(ynext-ynext_star)
        std::array<LD, N_eqs> abs_delta;
        
        std::array<LD, N_eqs> ynext;//this is here to hold the prediction
        std::array<LD, N_eqs> ynext_star;//this is here to hold the second prediction
    
        
        /*--These are specific to Rosenbrock methods*/
        std::array<LD, N_eqs> dfdt; 
        //define the coefficient. This will become (I-\gamma*h*J). _inv is its inverse
        std::array<std::array<LD, N_eqs>, N_eqs> _inv;
        // There are for the LUP-decomposition of (I-\gamma*h*J) 
        std::array<std::array<LD, N_eqs>, N_eqs> L;
        std::array<std::array<LD, N_eqs>, N_eqs> U;
        std::array<int,N_eqs> P;
        //lu_sol will capture the sulution of (I-\gamma*h*J)* k = rhs (i.e. k = (I-\gamma*h*J)^{-1} rhs)
        std::array<LD, N_eqs> lu_sol;
        std::array<std::array<LD, N_eqs>, N_eqs> J;//this is here to hold values of the Jacobian

        void LU();
        void calc_k();
        void calc_Jk();
        
        void sum_ak(const unsigned int& stage); // calculate sum_j a_{ij}*k_j and passit to this->ak
        void sum_gk(const unsigned int& stage); // calculate sum_j a_{ij}*k_j and passit to this->ak
        void sum_bk();// calculate sum_i b_i*k_i and passit to this->bk 
        
        void step_control();//adjust stepsize until error is acceptable
    public:

        Ros(const diffeq& dydt, const std::array<LD, N_eqs> &init_cond, LD tmax, const parameters<LD>& opt=default_parameters<LD>)
                    : dydt(dydt),Jac(jacobian(dydt)),params(opt) {
                // if some parameter in opt does not have a value, use the corresponding parameter from default.default_parameters 
                if(!params.initial_step_size.has_value()){params.initial_step_size=default_parameters<LD>.initial_step_size.value();}
                if(!params.minimum_step_size.has_value()){params.minimum_step_size=default_parameters<LD>.minimum_step_size.value();}
                if(!params.maximum_step_size.has_value()){params.maximum_step_size=default_parameters<LD>.maximum_step_size.value();}
                if(!params.maximum_No_steps.has_value()){params.maximum_No_steps=default_parameters<LD>.maximum_No_steps.value();}
                if(!params.absolute_tolerance.has_value()){params.absolute_tolerance=default_parameters<LD>.absolute_tolerance.value();}
                if(!params.relative_tolerance.has_value()){params.relative_tolerance=default_parameters<LD>.relative_tolerance.value();}
                if(!params.beta.has_value()){params.beta=default_parameters<LD>.beta.value();}
                if(!params.fac_max.has_value()){params.fac_max=default_parameters<LD>.fac_max.value();}
                if(!params.fac_min.has_value()){params.fac_min=default_parameters<LD>.fac_min.value();}

                reset(init_cond,tmax,opt); 
            };

        
        ~Ros()=default;

        const std::vector<LD>& get_t() const { return time; }
        const std::vector<LD>& get_t(const unsigned int& step) const { return time.at(step); }
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
        void reset(const std::array<LD,N_eqs>& init_cond, LD tmax, const parameters<LD>& opt=default_parameters<LD>);

        //generally helpful, but not very important
        auto get_current_step() const {return time.size();}
        auto get_current_step_size() const {return h_acc;}
        const parameters<LD>& get_parameters() const {return params;}
};



#endif
