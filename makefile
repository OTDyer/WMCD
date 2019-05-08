SHELL := /bin/bash

Nruns = 100
Batch = 1
Index = 1
Continue = 0
# set Continue=0 unless you are continuing from a previous run
# if you don't you risk appending data where you don't want to
Param = 0.0
# set Param to 0.0 if using parameters as set by input.c
# if want to change one from outside, adjust the line under the 
# call to input() as appropriate and set Param to desired value


# Commands for compiling code

wmcd:
	gcc -Wall -lm WMCD_main.c WMCD_maths.c WMCD_phys.c WMCD_transforms.c WMCD_input.c WMCD_output.c pcg_basic.c -o WMCD.out


wmcdicc:
	icc -O3 WMCD_main.c WMCD_maths.c WMCD_phys.c WMCD_transforms.c WMCD_input.c WMCD_output.c pcg_basic.c -o WMCD.out


# Commands for running code

run:
	./WMCD.out $(Batch) $(Index) $(Continue) $(Param)


multirun:
	for ((n=1; n<= $(Nruns); n++)); do printf "Run $$n of $(Nruns)\n"; ./WMCD.out $$n $(Nruns) 0 0.0; done



