#set PATH variable to working directory


# METHOD=ROS34PW2
METHOD=RODASPR2
# METHOD=ROS3w


# LONG=  #leave empty to use doubles
LONG=long #set to long , in order to use long doubles

 

Path=$(PWD)


CC=g++
FLG=  -std=c++17  -I "$(Path)" -lm  -DMETHOD=$(METHOD) -DLONG=$(LONG)  -O3


all: TwoState.run

TwoState.run: makefile TwoState.cpp 
	$(CC) -o "$(Path)/TwoState.run" "$(Path)/TwoState.cpp" $(FLG) 


clean:
	@[ -f "$(Path)/TwoState.run" ] && rm "$(Path)/TwoState.run" || true
	 