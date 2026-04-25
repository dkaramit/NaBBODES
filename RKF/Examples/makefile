METHOD=DormandPrince #higher rel_tol and abs_tol
# METHOD=RKF45  #low rel_tol and abs_tol, but the result seems better than DormandPrince
# METHOD=CashKarp #similar to RKF45


# LONG= #leave empty to use doubles
LONG=long #set to long , in order to use long doubles


rootDir="../"
IncludeDir=$(rootDir)
IncludeHPP = $(shell find $(rootDir) -type f -name '*.hpp')

Path=$(CURDIR)
CC=g++

FLG=  -std=c++20  -I "$(rootDir)" -lm  -DMETHOD=$(METHOD) -DLONG=$(LONG)  -O3 -Wpedantic


all: RKF.run Example.run

RKF.run: $(IncludeHPP) RKF.cpp makefile
	$(CC) -o "$(Path)/RKF.run" "$(Path)/RKF.cpp" $(FLG) 

Example.run: $(IncludeHPP) Example.cpp makefile
	$(CC) -o "$(Path)/Example.run" "$(Path)/Example.cpp" $(FLG) 


clean:
	rm -f *.run
	