all: cpu-app mic-app
	
cpu-app: exp1.cpp
	mpiicpc -o $@ $<

mic-app: exp1.cpp
	mpiicpc -o $@ $< -mmic
	scp $@ mic0:~/
	scp $< mic0:~/

run:
	mpiexec -host node6 -np 1 ./cpu-app : -host mic0 -np 1 ~/mic-app

help:
	export I_MPI_FABRICS=shm:tcp
	export I_MPI_MIC=1
