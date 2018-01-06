all: GA_baseline

GA_baseline: GAtuantu.cpp
	icpc -o $@ $< -D numMachine=500 -D numTask=500 -D sizePop=500 -D generation=500 -D rateMutan=0.005 -qopt-report-file=$@.optrpt -qopenmp -mmic

run:
	./GA_baseline

run_test_gcc:
	g++ -o GA_baseline GAtuantu.cpp -D numMachine=500 -D numTask=500 -D sizePop=200 -D generation=500 -D rateMutan=0.005
	./GA_baseline
	g++ -o GA_baseline GAtuantu.cpp -D numMachine=500 -D numTask=500 -D sizePop=400 -D generation=500 -D rateMutan=0.005
	./GA_baseline
	g++ -o GA_baseline GAtuantu.cpp -D numMachine=500 -D numTask=500 -D sizePop=600 -D generation=500 -D rateMutan=0.005
	./GA_baseline
	g++ -o GA_baseline GAtuantu.cpp -D numMachine=500 -D numTask=500 -D sizePop=800 -D generation=500 -D rateMutan=0.005
	./GA_baseline
	
	g++ -o GA_baseline GAtuantu.cpp -D numMachine=500 -D numTask=500 -D sizePop=1000 -D generation=500 -D rateMutan=0.005
	./GA_baseline
	g++ -o GA_baseline GAtuantu.cpp -D numMachine=500 -D numTask=500 -D sizePop=1200 -D generation=500 -D rateMutan=0.005
	./GA_baseline
	g++ -o GA_baseline GAtuantu.cpp -D numMachine=500 -D numTask=500 -D sizePop=1400 -D generation=500 -D rateMutan=0.005
	./GA_baseline
	g++ -o GA_baseline GAtuantu.cpp -D numMachine=500 -D numTask=500 -D sizePop=1600 -D generation=500 -D rateMutan=0.005
	./GA_baseline
	g++ -o GA_baseline GAtuantu.cpp -D numMachine=500 -D numTask=500 -D sizePop=1800 -D generation=500 -D rateMutan=0.005
	./GA_baseline
	g++ -o GA_baseline GAtuantu.cpp -D numMachine=500 -D numTask=500 -D sizePop=2000 -D generation=500 -D rateMutan=0.005
	./GA_baseline
