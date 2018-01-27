all: cpu-GA mic-GA
	
cpu-GA: Distributed_GA.cpp
	mpiicpc -o $@ $<

mic-GA: Distributed_GA.cpp
	mpiicpc -o $@ $< -mmic
	scp $@ mic0:~/
	scp $< mic0:~/

run:
	mpiexec -host node6 -np 1 ./cpu-GA : -host mic0 -np 1 ~/mic-GA : -host mic0 -np 1 ~/mic-GA

help:
	export I_MPI_FABRICS=shm:tcp
	export I_MPI_MIC=1
