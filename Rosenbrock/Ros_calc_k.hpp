#ifndef Ros_calc_k
#define Ros_calc_k
#include "Ros_class.hpp"
#include<vector>


namespace rosenbrock{

/*-----------------------Begin: calc_Jk---------------------------------*/
template<class LD, class RK_method, step_controllers step_controller> 
void Solver<LD, RK_method, step_controller>::calc_Jk(){
    /*
    Calculate product of matrix with vector. 
    We use it to get J*gk.
    */
    for(size_t i=0; i<N_eqs ; i++){
        Jk[i]=0;
        for(size_t j=0; j<N_eqs ; j++){ Jk[i]+=J[i][j]*gk[j]; }
    }         
}
/*-----------------------End: calc_Jk---------------------------------*/

/*-----------------------Begin: calc_k---------------------------------*/
template<class LD, class RK_method, step_controllers step_controller> 
void Solver<LD, RK_method, step_controller>::calc_k(){
    // since LU decomposition is not updated as you try to find suitable step, put it ouside the step control loop (should be faster)!
    // LU();
    // calculate k for the other stages
    std::vector<LD> yn(N_eqs,0);
    // rhs is the rhs side of the equation to be solved by LU.  
    std::vector<LD> rhs(N_eqs,0);
    // to make it more clear, we are going to separate the rhs in three different parts
    std::vector<LD> rhs1(N_eqs,0),rhs2(N_eqs,0);

    for(size_t stage = 0; stage < RK_method::s; stage++){
        sum_ak(stage);
        sum_gk(stage);
        // setup the argument for dydt (it is evaluated at y_n+\sum a*k )
        for(size_t eq=0; eq<N_eqs ; eq++){yn[eq]=yprev[eq]+ak[eq];}
        /*--- Get the rhs terms ---*/
        // first term
        dydt(rhs1, yn, tn+RK_method::c[stage]*h_trial);
        // second term
        for(size_t eq=0; eq<N_eqs ; eq++) { rhs2[eq] =  (RK_method::gamma+sum_gamma[stage])*dfdt[eq];}
        // third term
        calc_Jk();
        // then the rhs becomes
        for(size_t eq=0; eq<N_eqs ; eq++) { rhs[eq] = h_trial*rhs1[eq] + h_trial*h_trial*rhs2[eq] + h_trial*Jk[eq];}
        /*------------------------*/


        // Solve  L*U*P*lu_sol=rhs.
        Solve_LU<LD>(L,U,P,rhs, lu_sol);
        for(size_t eq = 0; eq < N_eqs; eq++ ){ k[eq][stage]=lu_sol[eq]; }
    }
}
/*-----------------------End: calc_k---------------------------------*/

}

#endif