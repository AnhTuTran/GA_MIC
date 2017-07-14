all: GA_MIC

GA_MIC: GA_baseline.cpp
	icpc -o $@ $< -qopenmp -D numMachine=500 -D numTask=500 -D sizePop=512 -D generation=100 -D rateMutan=0.005 -qopt-report=5 -qopt-report-file=$@.optrpt

run:
	./GA_MIC
