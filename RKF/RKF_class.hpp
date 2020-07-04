#ifndef RKF_class
#define RKF_class

//This is a general implementation of explicit embedded RK solver of
// a system of differential equations in the interval [0,tmax].

/*
diffeq is a class of the system of  equations to be solved 
N_eqs is ten number of equations to be solved
RKF_method is the method (the DormandPrince seems to be the standard here)
N_out number of output points (to be taken in intervals of approximately  1/(N_out-1) )
*/

#define RKF_Template template<class diffeq, int N_eqs, class RKF_method, int N_out, class LD> 
#define RKF_Namespace RKF<diffeq, N_eqs, RKF_method, N_out, LD>


RKF_Template
class RKF
{
    /*      */
private://There is no reason to make things private (if you break it it's not my fault)... 
    
public:
    //Inputs. The initial condition is given as a Array (the type is users choice as long as it can be called with [])
    diffeq dydt;
    RKF_method method;
    LD tmax, h, hmin, hmax, abs_tol, rel_tol, beta, fac_max, fac_min;
    int max_N;
    LD h_old,delta_acc, delta_rej;//these will be initialized at the beginning of next_step
    
    //things that we'll need
    int current_step;
    bool h_stop;//h_stop becomes true when suitable stepsize is found.    
    LD tn;
    LD tmp_sol[N_eqs];
    
    
    std::vector<LD> solution[N_eqs];
    std::vector<LD> error[N_eqs];
    std::vector<LD> time;
    std::vector<int> hist;

    std::vector<LD> Deltas;//this will hold the Dy/scale (this is what we try to send to 1 by adjusting the stepsize )
    

    std::vector<LD> time_full;
    std::vector<LD> solution_full[N_eqs];
    std::vector<LD> error_full[N_eqs];
    

   
    //these are here to hold the k's, sum_i b_i*k_i, sum_i b_i^{\star}*k_i, and sum_j a_{ij}*k_j 
    LD **k;
    LD ak[N_eqs], bk[N_eqs],bstark[N_eqs];
    LD abs_delta[N_eqs];

    LD ynext[N_eqs];//this is here to hold the prediction
    LD ynext_star[N_eqs];//this is here to hold the second prediction

    LD yn[N_eqs];//thi i here to hold accepted ynext (redundant, I'll remove it later)
    LD fyn[N_eqs];//this is here to get dydt in each step
    
    RKF(diffeq  dydt, LD (&init_cond)[N_eqs], LD tmax,
LD initial_step_size=1e-5, LD minimum_step_size=1e-11, LD maximum_step_size=1e-3,int maximum_No_steps=1000000, 
        LD absolute_tolerance=1e-8,LD relative_tolerance=1e-8,LD beta=0.85,LD fac_max=3, LD fac_min=0.3);
    
    ~RKF();

    /*-------------------it would be nice to have a way to define these sums more generaly-----------------*/
    void next_step();

    void calc_k();

    void sum_ak(int stage); // calculate sum_j a_{ij}*k_j and passit to this->ak
    void sum_bk();// calculate sum_i b_i*k_i and passit to this->bk 
    void sum_bstark();// calculate sum_i b^{\star}_i*k_i and passit to this->bk
    

    
    void step_control();//adjust stepsize until error is acceptable
    
    // the default value of _full_ is false. So if you use solve() the full output will not be saved. To get the full output call solve(true)
    void solve(bool _full_=false);


};



#endif
