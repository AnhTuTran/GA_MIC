all: GA_baseline

GA_baseline: GAtuantu.cpp
	icpc -xhost -o $@ $< -D numMachine=500 -D numTask=500 -D sizePop=512 -D generation=100

run:
	./GA_baseline
