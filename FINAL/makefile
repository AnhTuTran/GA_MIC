all: GA_baseline

GA_baseline: GAtuantu.cpp
	icpc -o $@ $< -D numMachine=500 -D numTask=500 -D sizePop=512 -D generation=100 -D rateMutan=0.005 -qopt-report-file=$@.optrpt

run:
	./GA_baseline
