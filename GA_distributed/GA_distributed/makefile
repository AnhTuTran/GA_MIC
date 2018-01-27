all: cpu-GA mic-GA
	
cpu-GA: Distributed_GA.cpp
	mpiicpc -o $@ $< -qopenmp -D numMachine=500 -D numTask=500 -D sizePop=200 -D generation=500 -D rateMutan=0.005

mic-GA: Distributed_GA.cpp
	mpiicpc -o $@ $< -mmic -qopenmp -D numMachine=500 -D numTask=500 -D sizePop=200 -D generation=500 -D rateMutan=0.005
	scp $@ mic0:~/
	scp $@ mic1:~/

run:
	mpiexec -host node6 -np 1 ./cpu-GA : -host mic0 -np 1 ~/mic-GA : -host mic1 -np 1 ~/mic-GA

help:
	export I_MPI_FABRICS=shm:tcp
	export I_MPI_MIC=1
