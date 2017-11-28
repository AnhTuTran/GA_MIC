all: GA_MIC

GA_MIC: GA_baseline.cpp
	icpc -o $@ $< -qopenmp -D numMachine=500 -D numTask=500 -D sizePop=500 -D generation=500 -D rateMutan=0.005 -qopt-report=5 -qopt-report-file=$@.optrpt

run:
	./GA_MIC

run_test:
	icpc -o GA_MIC GA_baseline.cpp -qopenmp -D numMachine=500 -D numTask=500 -D sizePop=200 -D generation=500 -D rateMutan=0.005
	./GA_MIC
	icpc -o GA_MIC GA_baseline.cpp -qopenmp -D numMachine=500 -D numTask=500 -D sizePop=400 -D generation=500 -D rateMutan=0.005
	./GA_MIC
	icpc -o GA_MIC GA_baseline.cpp -qopenmp -D numMachine=500 -D numTask=500 -D sizePop=600 -D generation=500 -D rateMutan=0.005
	./GA_MIC
	icpc -o GA_MIC GA_baseline.cpp -qopenmp -D numMachine=500 -D numTask=500 -D sizePop=800 -D generation=500 -D rateMutan=0.005
	./GA_MIC
	
	icpc -o GA_MIC GA_baseline.cpp -qopenmp -D numMachine=500 -D numTask=500 -D sizePop=1000 -D generation=500 -D rateMutan=0.005
	./GA_MIC
	icpc -o GA_MIC GA_baseline.cpp -qopenmp -D numMachine=500 -D numTask=500 -D sizePop=1200 -D generation=500 -D rateMutan=0.005
	./GA_MIC
	icpc -o GA_MIC GA_baseline.cpp -qopenmp -D numMachine=500 -D numTask=500 -D sizePop=1400 -D generation=500 -D rateMutan=0.005
	./GA_MIC
	icpc -o GA_MIC GA_baseline.cpp -qopenmp -D numMachine=500 -D numTask=500 -D sizePop=1600 -D generation=500 -D rateMutan=0.005
	./GA_MIC
	icpc -o GA_MIC GA_baseline.cpp -qopenmp -D numMachine=500 -D numTask=500 -D sizePop=1800 -D generation=500 -D rateMutan=0.005
	./GA_MIC
	icpc -o GA_MIC GA_baseline.cpp -qopenmp -D numMachine=500 -D numTask=500 -D sizePop=2000 -D generation=500 -D rateMutan=0.005
	./GA_MIC

run_test_thread:
	icpc -o GA_MIC GA_baseline.cpp -qopenmp -D numMachine=500 -D numTask=500 -D sizePop=1000 -D generation=500 -D rateMutan=0.005
	export KMP_PLACE_THREADS=60C,1t
	./GA_MIC
	export KMP_PLACE_THREADS=60C,2t
	./GA_MIC
	export KMP_PLACE_THREADS=60C,3t
	./GA_MIC
	export KMP_PLACE_THREADS=60C,4t
	./GA_MIC
env:
	export KMP_PLACE_THREADS=60C,1t
