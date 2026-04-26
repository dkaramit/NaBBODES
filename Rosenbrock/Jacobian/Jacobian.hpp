#ifndef Jac_head
#define Jac_head
#include<functional>
#include<array>
// This is an example Jacobian class.


template<int N_eqs, class LD>
class Jacobian{
    private:
        using diffeq = std::function<void(std::array<LD, N_eqs> &lhs, const std::array<LD, N_eqs> &y, const LD& t)>;
        
        LD h;
        diffeq dydt;
        
        public:
        
        Jacobian(const diffeq& dydt, LD h=1e-10):  dydt(dydt), h(h) {};
        ~Jacobian()=default;
        
        void operator()(std::array<std::array<LD, N_eqs>, N_eqs> &J, std::array<LD, N_eqs> &dfdt, const std::array<LD, N_eqs> &y, const LD& t)const{
            
            std::array<LD, N_eqs> y0,y1,dydt0,dydt1;

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


#endif
