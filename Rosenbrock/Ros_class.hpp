#ifndef Ros_class
#define Ros_class

//This is a general implementation of explicit embedded RK solver of
// a system of differential equations in the interval [0,1].


_Ros_template_//Note that you can use template to pass the method
class Ros
{
    /*      */
private://There is no reason to make things private (if you break it it's not my fault)... 
    
public:
    //Inputs. The initial condition is given as a Array (the type is users choice as long as it can be called with [])
    diffeq dydt;
    RK_method method;
    jacobian Jac;
    


    // h0 is the current stepsize (used to update the stepsize)
    // h1 the previous accepted stepsize.
    LD h0,h1, hmin, hmax, abs_tol, rel_tol, beta, fac_max;
    int max_N;
    
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


   
    //these are here to hold the k's, sum_i b_i*k_i, sum_i b_i^{\star}*k_i, and sum_j a_{ij}*k_j 
    LD **k;
    LD ak[N_eqs],gk[N_eqs],Jk[N_eqs], bk[N_eqs],bstark[N_eqs];
    
    // abs_delta=fabs(ynext-ynext_star)
    LD abs_delta[N_eqs];
    // regulated_delta = (1-gam*h*J)^{-1}*abs_delta
    LD regulated_delta[N_eqs];

    LD ynext[N_eqs];//this is here to hold the prediction
    LD ynext_star[N_eqs];//this is here to hold the second prediction

    
    LD yn[N_eqs];//this is here to hold values I might need
    LD dfdt[N_eqs];//this is here to hold values I might need
    
    /*--These are specific to Rosenbrock methods*/

    //define the coefficient. This will become (I-\gamma*h*J). _inv is its inverse
    LD _coeff[N_eqs][N_eqs], _inv[N_eqs][N_eqs] ;
    // There are for the LUP-decomposition of (I-\gamma*h*J) 
    LD L[N_eqs][N_eqs];
    LD U[N_eqs][N_eqs];
    int P[N_eqs];

    // fyn will be passed to dydt and get the result. rhs is the rhs side of the equation to be solved by LU.  
    LD fyn[N_eqs],rhs[N_eqs];

    // to make it more clear, we are going to separate the rhs in three different parts
    LD rhs1[N_eqs],rhs2[N_eqs];

    //lu_sol will capture the sulution of (I-\gamma*h*J)* k = rhs (it is basically k)
    LD lu_sol[N_eqs];

    // need this to sore the sum over \gammas
    LD *sum_gamma;
       
    LD J[N_eqs][N_eqs];//this is here to hold values I might need
    /*----------------------------------------------------------------------------------------------------*/
    Ros(diffeq dydt, LD (&init_cond)[N_eqs], 
        LD initial_step_size=1e-5, LD minimum_step_size=1e-11, LD maximum_step_size=1e-3,int maximum_No_steps=1000000, 
        LD absolute_tolerance=1e-15,LD relative_tolerance=1e-15,LD beta=0.85,LD fac_max=3);
    
    ~Ros();

    /*-------------------it would be nice to have a way to define these sums more generaly-----------------*/
    void next_step();
    void LU();

    void calc_k();
    void calc_Jk();

    void sum_ak(int stage); // calculate sum_j a_{ij}*k_j and passit to this->ak
    void sum_gk(int stage); // calculate sum_j a_{ij}*k_j and passit to this->ak
    void sum_bk();// calculate sum_i b_i*k_i and passit to this->bk 
    void sum_bstark();// calculate sum_i b^{\star}_i*k_i and passit to this->bk
    

    
    void step_control();//adjust stepsize until error is acceptable
    void solve(bool _full_=false);


};



#endif