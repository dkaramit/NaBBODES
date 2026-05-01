#ifndef Jac_head
#define Jac_head
#include<functional>
#include<vector>
#include<cmath>

namespace rosenbrock{

// This is the default Jacobian class.

template<class LD>
class Jacobian{
    private:
        using diffeq = std::function<void(std::vector<LD> &lhs, const std::vector<LD> &y, const LD& t)>;
        
        LD h;
        diffeq dydt;
        
        public:
        
        Jacobian(const diffeq& dydt, const LD& h=1e-8):  dydt(dydt), h(h) {};
        ~Jacobian()=default;
        
        void update_step_size(const LD& h){this->h=h;}

        void operator()(std::vector<std::vector<LD>> &J, std::vector<LD> &dfdt, const std::vector<LD> &y, const LD& t)const{
            
            unsigned int N_eqs=y.size();


            std::vector<LD> y0(N_eqs),y1(N_eqs),dydt0(N_eqs),dydt1(N_eqs);

            // you can use something like this to scale the stepsize according to the scale of t
            LD a=h+h*t;
            // take the time derivative
            dydt(dydt0,y,t-a);
            dydt(dydt1,y,t+a);
            for (unsigned int i = 0; i < N_eqs; i++){ 
                dfdt[i]=(dydt1[i]-dydt0[i])/(2*a); 
                // since you run this loop, initialize y0 and y1
                y0[i]=y[i]; y1[i]=y[i];
            }

            // take the derivatives over y
            for (unsigned int j = 0; j < N_eqs; j++){
                // you can use something like this to scale the stepsize according to the scale of y[j]
                a=h+h*std::abs(y0[j]);
                
                y0[j]=y0[j]-a;
                y1[j]=y1[j]+a;
                
                dydt(dydt0,y0,t);
                dydt(dydt1,y1,t);

                for (unsigned int i = 0; i < N_eqs; i++){J[i][j]=(dydt1[i]-dydt0[i])/(2*a);}
                // if(isnan(J[i][j])){J[i][j]=0;}

                y0[j]=y[j];
                y1[j]=y[j];
            }
    };





};

}

#endif
