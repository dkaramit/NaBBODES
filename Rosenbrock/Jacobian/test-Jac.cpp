#include<iostream>
#include<vector>
#include"Jacobian.hpp"
#define _class //run it with class with overloaded operator()
// #define _function //run it with function pointer

#define LD double

using namespace rosenbrock; 




#ifdef _function
typedef void (*diffeq)(std::vector<LD> &lhs, const std::vector<LD> &y  , const LD& t); // define function pointer
void sys(std::vector<LD> &lhs, const std::vector<LD> &y  , const LD& t)
{

    lhs[0]=y[1]*y[0]*t;
    lhs[1]=y[0]*y[1]+t;

}
#endif


#ifdef _class
class Cdiffeq{
    public:
    Cdiffeq(){};
    ~Cdiffeq(){};
    void operator()(std::vector<LD> &lhs, const std::vector<LD> &y  , const LD& t)
    {
        lhs[0]=y[1]*y[0]*t;
        lhs[1]=y[0]*y[1]+t;
    }


};
#endif


using std::cout;
using std::endl;

unsigned int n_eqs = 2;

int main(){
    #ifdef _class
        Cdiffeq dydt;
        Jacobian<LD> jac(dydt);
    #endif

    #ifdef _function
        Jacobian<LD> jac(sys);
    #endif

    std::vector<LD> dfdt(n_eqs,0);
    // Matrix J={{0,0},{0,0}};
    std::vector<std::vector<LD>> J(n_eqs,std::vector<LD>(n_eqs,0));

    std::vector<LD> y={1,4.2};
    LD t=0.3;

    jac(J,dfdt,y,t);

    cout<<"dfdt=[";
    for(int i=0; i<n_eqs ; i++){  cout<<dfdt[i]; if(i!=n_eqs-1){cout<<",";}   }
    cout<<"]"<<endl;


    cout<<"J=["<<endl;
    
    for (int i = 0; i < n_eqs; i++){
         cout<<"[";
        for (int j = 0; j < n_eqs; j++) {cout<<J[i][j]; if(j!=n_eqs-1){cout<<",";} }
          cout<<"]";if(i!=n_eqs-1){cout<<",";}cout<<endl;
        
        
    }  
    cout<<"]"<<endl;
    





    return 0;
}