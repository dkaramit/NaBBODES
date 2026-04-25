#ifndef RKF_class
#define RKF_class

#include<array>
#include<vector>
#include<functional>

// struct for the parameters of the RKF algorithm. 
// It is useful because we can now pass them named parameters!
template<class LD>
struct parameters{
    const LD initial_step_size=1e-5; 
    const LD minimum_step_size=1e-11; 
    const LD maximum_step_size=1e-3;
    const unsigned int maximum_No_steps=1000000; 
    const LD absolute_tolerance=1e-8;
    const LD relative_tolerance=1e-8;
    const LD beta=0.85;
    const LD fac_max=3; 
    const LD fac_min=0.3;
};

//This is a general implementation of explicit embedded RK solver of
// a system of differential equations in the interval [0,tmax].

/*
diffeq is a class of the system of  equations to be solved 
N_eqs is ten number of equations to be solved
RKF_method is the method (the DormandPrince seems to be the standard here)
*/

template<unsigned int N_eqs, class RK_method, class LD>
class RKF{
    private:
        using diffeq=std::function<void(std::array<LD, N_eqs> &lhs, std::array<LD, N_eqs> &y, LD t)>;

        //these are not constant because reset can update them
        LD hmin, hmax, abs_tol, rel_tol, beta, fac_max, fac_min;

        LD h_old,h_trial,h_acc,delta_acc, delta_rej;//these will be initialized at the beginning of next_step
        unsigned int max_N;
        bool h_stop;//h_stop becomes true when suitable stepsize is found.    
        
        //Inputs. The initial condition is given as a Array (the type is users choice as long as it can be called with [])
        const diffeq dydt;
        // integrate up to t=tmax
        LD tmax;
        
        LD tn;

        std::array<std::array<LD,RK_method::s>,N_eqs> k;
        
        //these are here to hold the k's, sum_i b_i*k_i, sum_i b_i^{\star}*k_i, and sum_j a_{ij}*k_j 
        std::array<LD,N_eqs> ak, bk, bstark;
        // abs_delta=abs(ynext-ynext_star)
        std::array<LD,N_eqs> abs_delta;
        
        std::array<LD, N_eqs> yprev;
        std::array<LD,N_eqs> ynext;//this is here to hold the prediction
        std::array<LD,N_eqs> ynext_star;//this is here to hold the second prediction
        
        std::vector<LD> time;
        std::array<std::vector<LD>, N_eqs> solution;
        std::array<std::vector<LD>, N_eqs> error;

        void calc_k(); //calculate the values of k
        void sum_ak(const unsigned int& stage); // calculate sum_j a_{ij}*k_j and pass it to this->ak
        void sum_bk();// calculate sum_i b_i*k_i and passit to this->bk 
        void step_control();//adjust stepsize until error is acceptable

        public:
        
        RKF(const diffeq&  dydt, const std::array<LD,N_eqs>& init_cond, const LD& tmax, 
            const parameters<LD>& opt=parameters<LD>{}): dydt(dydt) { reset(init_cond,tmax,opt); };
        
        ~RKF()=default;

        void next_step();
        void solve();

        
        const std::vector<LD>& get_t() const { return time; }
        //access the array of solution[eq]
        const std::vector<LD>& get_solution(const unsigned int& eq) const { return solution[eq]; }
        //access the element solution[eq][step]
        const LD get_solution(const unsigned int& eq, const unsigned int& step) const { return solution[eq][step]; }
        
        //access the array of error[eq]
        const std::vector<LD>& get_error(const unsigned int& eq) const { return error[eq]; }
        //access the element error[eq][step]
        const LD get_error(const unsigned int& eq, const unsigned int& step) const { return error[eq][step]; }

        void set_parameters(const parameters<LD>& opt=parameters<LD>{});
        void reset(const std::array<LD,N_eqs>& init_cond, const LD& tmax, const parameters<LD>& opt=parameters<LD>{});
};





#endif
