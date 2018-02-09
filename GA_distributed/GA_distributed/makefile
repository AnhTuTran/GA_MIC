all: cpu-GA mic-GA
	
cpu-GA: Distributed_GA.cpp
	mpiicpc -o $@ $< -qopenmp -D numMachine=500 -D numTask=500 -D sizePop=200 -D generation=500 -D rateMutan=0.005

mic-GA: Distributed_GA.cpp
	mpiicpc -o $@ $< -mmic -qopenmp -D numMachine=500 -D numTask=500 -D sizePop=200 -D generation=500 -D rateMutan=0.005
	#scp $@ mic0:~/
	#scp $@ mic1:~/

run:
	mpiexec -host node6 -np 1 ./cpu-GA : -host mic0 -np 1 ~/mic-GA : -host mic1 -np 1 ~/mic-GA

help:
	export I_MPI_FABRICS=shm:tcp
	export I_MPI_MIC=1

run_test:
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=200
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=200
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=400
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=400
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=600
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=600
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=800
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=800
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=1000
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=1000
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=1200
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=1200
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=1400
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=1400
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=1600
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=1600
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=1800
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=1800
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=2000
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=2000
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA

run_test_big:
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=2200
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=2200
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=2400
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=2400
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=2600
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=2600
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=2800
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=2800
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=3000
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=500 -D numTask=500 -D generation=500 -D rateMutan=0.005 -D sizePop=3000
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA


run_test_1000:
	#mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=200
	#mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=200
	#scp mic-GA mic0:~/GA_distributed
	#scp mic-GA mic1:~/GA_distributed
	#mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	#
	#mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=400
	#mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=400
	#scp mic-GA mic0:~/GA_distributed
	#scp mic-GA mic1:~/GA_distributed
	#mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	#
	#mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=600
	#mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=600
	#scp mic-GA mic0:~/GA_distributed
	#scp mic-GA mic1:~/GA_distributed
	#mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	#
	#mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=800
	#mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=800
	#scp mic-GA mic0:~/GA_distributed
	#scp mic-GA mic1:~/GA_distributed
	#mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	#
	#mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=1000
	#mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=1000
	#scp mic-GA mic0:~/GA_distributed
	#scp mic-GA mic1:~/GA_distributed
	#mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	#
	#mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=1200
	#mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=1200
	#scp mic-GA mic0:~/GA_distributed
	#scp mic-GA mic1:~/GA_distributed
	#mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=1400
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=1400
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=1600
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=1600
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=1800
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=1800
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	
	mpiicpc -o cpu-GA Distributed_GA.cpp -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=2000
	mpiicpc -o mic-GA Distributed_GA.cpp -mmic -qopenmp -D numMachine=1000 -D numTask=1000 -D generation=500 -D rateMutan=0.005 -D sizePop=2000
	scp mic-GA mic0:~/GA_distributed
	scp mic-GA mic1:~/GA_distributed
	mpiexec -host node25 -np 1 ./cpu-GA : -host mic0 -np 1 ~/GA_distributed/mic-GA : -host mic1 -np 1 ~/GA_distributed/mic-GA
	
	
