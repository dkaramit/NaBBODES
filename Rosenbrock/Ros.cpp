// This is how you run RK. 
#include<iostream>
#include<fstream>
#include<cmath>
#include"Ros.hpp"
#include "Jacobian/Jacobian.hpp"//this is treated as user input, since one may have an analytic form.
using std::cout;
using std::endl;
/*--------------------------------------------------------------------------------------------------------*/
class ROS3w{
    public:
        const int s=3;
        const int p=2;
        double c[3]={0,2/3.,4/3.};
        double b[3]={0.25,0.25,0.5 };
        double bstar[3]={ 0.746704703274 , 0.1144064078371 , 0.1388888888888 };
        double gamma=0.4358665215084;
        
        double a[3][3] , g[3][3];

        ~ROS3w(){};
        ROS3w(){
            
            for (int i = 0; i < s; i++){ for (int j = 0; j < s; j++){a[i][j]=0; g[i][j]=0; } }
            
            a[1][0]=2/3.;
            a[2][0]=2/3.;
            a[2][1]=2/3.;
            
            g[1][0]=0.3635068368900;
            g[2][0]=-0.8996866791992;
            g[2][1]=-0.1537997822626;
        };
};
/*--------------------------------------------------------------------------------------------------------*/
class ROS34PW2{
    public:
        const int s=4;
        const int p=3;
        double c[4];// c[i]=sum_{j} a[i][j]
        double b[4]={2.4212380706095346e-1 , -1.2232505839045147 , 1.5452602553351020 , 4.3586652150845900e-1};
        double bstar[4]={ 3.7810903145819369e-1 , -9.6042292212423178e-2 , 0.5 , 2.1793326075422950e-1};
        double gamma=4.3586652150845900e-1; 
        
        double a[4][4] , g[4][4];

        ~ROS34PW2(){};
        ROS34PW2(){
            
            for (int i = 0; i < s; i++){c[i]=0; for (int j = 0; j < s; j++){  a[i][j]=0; g[i][j]=0; } }
            
            a[1][0]=8.7173304301691801e-1;
            
            a[2][0]=8.4457060015369423e-1;
            a[2][1]=-1.1299064236484185e-1;
           
            a[3][2]=1.;

            for (int i = 0; i < s; i++){for (int j = 0; j < s; j++){  c[i]+=a[i][j]; } }
           

            g[1][0]=-8.7173304301691801e-1;
            
            g[2][0]=-9.0338057013044082e-1;
            g[2][1]=5.4180672388095326e-2;
            
            g[3][0]=2.4212380706095346e-1;
            g[3][1]=-1.2232505839045147;
            g[3][2]=5.4526025533510214e-1;
           
        };
};
/*--------------------------------------------------------------------------------------------------------*/
class RODASPR2{
    public:
        const int s=6;
        const int p=4;
        double c[6];// c[i]=sum_{j} a[i][j]
        double b[6]={5.1944159827788361e-1,3.9955867781540699e-2,-4.7356407303732290e-1,9.4907420451383284e-1,-3.4740759753593431e-1,3.125e-1  };
        double bstar[6]={-1.7746585073632790e-1,-5.8241418952602364e-1,6.8180612588238165e-1,7.6557391437996980e-1,3.125e-1,0 };
        double gamma=3.125e-1; 
        
        double a[6][6] , g[6][6];

        ~RODASPR2(){};
        RODASPR2(){
            
            for (int i = 0; i < s; i++){c[i]=0; for (int j = 0; j < s; j++){  a[i][j]=0; g[i][j]=0; } }
            a[1][0]=9.375e-1;
            
            a[2][0]= -4.7145892646261345e-2; 
            a[2][1]=5.4531286650471122e-1;
            
            
            a[3][0]= 4.6915543899742240e-1;    
            a[3][1] = 4.4490537602383673e-1; 
            a[3][2] =-2.2498239334061121e-1 ;
            
            a[4][0]= 1.0950372887345903;     
            a[4][1] = 6.3223023457294381e-1;  
            a[4][2] = -8.9232966090485821e-01; 
            a[4][3]= 1.6506213759732410e-1;

            a[5][0]= -1.7746585073632790e-1;    
            a[5][1] =-5.8241418952602364e-1; 
            a[5][2] = 6.8180612588238165e-1;
            a[5][3]= 7.6557391437996980e-1;    
            a[5][4] = 3.125e-1;

            for (int i = 0; i < s; i++){for (int j = 0; j < s; j++){  c[i]+=a[i][j]; } }
           

            g[1][0]=-9.3750000000000000e-01;
            
            g[2][0]=-9.7580572085994507e-02;
            g[2][1]=-5.8666328499964138e-01;
            
            g[3][0]=-4.9407065013256957e-01;
            g[3][1]=-5.6819726428975503e-01;
            g[3][2]=5.0318949274167679e-01;
            
            g[4][0]=-1.2725031394709183;
            g[4][1]=-1.2146444240989676;
            g[4][2]=1.5741357867872399;
            g[4][3]=6.0051177678264578e-01;
            
            g[5][0]=6.9690744901421153e-01;
            g[5][1]=6.2237005730756434e-01;
            g[5][2]=-1.1553701989197045;
            g[5][3]=1.8350029013386296e-01;
            g[5][4]=-6.5990759753593431e-01;
           
           
        };
};
/*--------------------------------------------------------------------------------------------------------*/




// this is how the diffeq should look like
#define n_eqs 3 //number of equations
typedef double Array[n_eqs];//define an array type of length n_eqs
typedef void (*diffeq)(Array &lhs, Array &y  , double t );
//-------------------------------------------------------------------------//

void sys( Array &lhs, Array &y  , double t )
        {
            //lhs is an array that gets the return value (the left hand side of the equation)
            //y is an array with values of y
            //t is the value of the variable t
            
            lhs[0]=-20*y[0]*pow(t,3.);
            lhs[1]=5*y[0]*pow(t,2)+2*(-pow( y[1],2  )+pow( y[2],2 ) )*t;
            lhs[2]=15*y[0]*pow(t,2)+2*(pow( y[1],2  )-pow( y[2],2 ) )*t;

        };


#define initial_step_size 1e-3
#define minimum_step_size 1e-15
#define maximum_step_size 1e-2
#define maximum_No_steps 1000000
#define absolute_tolerance 1e-12
#define relative_tolerance 1e-12
#define beta 0.85
#define fac_max 3

// #define METHOD ROS3w //2nd order
// #define METHOD ROS34PW2 //3rd order
#define METHOD RODASPR2  //4th order

int main(int argc, const char** argv) {
    
    Array lhs;
    Array y0;
    y0[0]=5;
    y0[1]=10;
    y0[2]=0;
    Jacobian<diffeq,n_eqs> jac(sys);

    // this is not looing good.. Just pass the number of equations as argument to get rid of Array and Matrix...
    Ros<diffeq,n_eqs, METHOD ,Jacobian<diffeq,n_eqs> > System(sys,y0, 
     initial_step_size,  minimum_step_size,  maximum_step_size, maximum_No_steps, 
     absolute_tolerance, relative_tolerance, beta, fac_max);
    // System.next_step();
    System.solve();
    // 
    

    #if 0
    for (int i = 0; i < System.max_N; i++)
    {
       System.next_step();

       std::cout<<System.tn<<"\n";

       if(System.tn==1){break;} 
    }
    #endif

    
    #if 0
    for (int i = 0; i < System.current_step; i++)
    {

        std::cout<<i<<"  "<<System.steps[i]<<"\t\t\t\t";
        for(int eq =0; eq<n_eqs; eq++){ std::cout<< System.solution[eq][i]<< "\t\t\t"; }
        std::cout<<"\n";
        
        if(System.steps[i]==1){break;}
    }
    #endif


    std::ofstream f1,f2,f3,t,err;
    f1.open ("y1.dat");
    f2.open ("y2.dat");
    f3.open ("y3.dat");
    t.open ("t.dat");
    err.open ("err.dat");
    

   
    for (int i = 0; i < System.current_step; i++)
    {
        f1 << System.solution[0][i] ;
        f1 << "\n";
        f2 << System.solution[1][i] ;
        f2 << "\n";
        f3 << System.solution[2][i] ;
        f3 << "\n";
        t << System.steps[i] ;
        t << "\n";
        err << System.err[i] ;
        err << "\n";
        
        // std::cout<<System.steps[i]<<"\t"<< System.solution[0][i] << "\t"<< System.solution[1][i] << "\t"<< System.solution[2][i] <<"\n";
        if(System.steps[i]==1){break;}
    }

    f1.close();
    f2.close();
    f3.close();
    t.close();

    // std::cout<<System.current_step<<std::endl;

    return 0;
 }