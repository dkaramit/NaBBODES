#set PATH variable to working directory


METHOD=DormandPrince


# LONG= #leave empty to use doubles
LONG=long #set to long , in order to use long doubles



Path=$(PWD)
CC=g++
FLG=  -std=c++17  -I "$(Path)" -lm  -DMETHOD=$(METHOD) -DLONG=$(LONG)  -O3


all: RKF

RKF: METHOD.hpp  RKF_calc_k.hpp  RKF_class.hpp  RKF_costructor.hpp    RKF.hpp  RKF_step_control.hpp  RKF_steps.hpp  RKF_sums.hpp  RKF.cpp
	$(CC) -o "$(Path)/RKF" "$(Path)/RKF.cpp" $(FLG) 


clean:
	@[ -f "$(Path)/RKF" ] && rm "$(Path)/RKF" || true
	