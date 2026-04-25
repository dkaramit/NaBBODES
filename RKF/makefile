#set PATH variable to working directory


METHOD=DormandPrince #higher rel_tol and abs_tol
# METHOD=RKF45  #low rel_tol and abs_tol, but the result seems better than DormandPrince
# METHOD=CashKarp #similar to RKF45


# LONG= #leave empty to use doubles
LONG=long #set to long , in order to use long doubles



Path=$(CURDIR)
CC=g++

FLG=  -std=c++20  -I "$(Path)" -lm  -DMETHOD=$(METHOD) -DLONG=$(LONG)  -O3 -Wpedantic


all: RKF.run Example.run

RKF.run: METHOD.hpp  RKF_calc_k.hpp  RKF_class.hpp  RKF_reset.hpp\
	RKF.hpp  RKF_step_control-simple.hpp RKF_step_control-PI.hpp\
	RKF_steps.hpp  RKF_sums.hpp  RKF.cpp makefile
	$(CC) -o "$(Path)/RKF.run" "$(Path)/RKF.cpp" $(FLG) 

Example.run: METHOD.hpp  RKF_calc_k.hpp  RKF_class.hpp  RKF_reset.hpp\
	RKF.hpp  RKF_step_control-simple.hpp RKF_step_control-PI.hpp\
	RKF_steps.hpp  RKF_sums.hpp  Example.cpp makefile
	$(CC) -o "$(Path)/Example.run" "$(Path)/Example.cpp" $(FLG) 


clean:
	rm -f *.run
	