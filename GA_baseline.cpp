/**
GA tuan tu tren CPU
	Version 1.0
	Developer : Le Thanh Tan
				Chiem Thach Phat
*/
#pragma offload_attribute(push, target(mic))
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <string>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <omp.h>
#pragma offload_attribute(pop) 

using namespace std;

// define
#ifndef numTask
#define numTask				10		// number of Tasks, task 0, 1, ...
#endif

#ifndef numMachine
#define numMachine			10		// number of Machines, machine 0, 1, ...,		criteria: <=1000
#endif

#ifndef sizePop
#define sizePop				1024 	// size of population
#endif

#define edgeSizePop			64		// edge of population

#ifndef generation
#define	generation			1000	// number of generations
#endif

#define timeStopGa  		3600	// thoi gian toi da chay giai thuat ga, tinh bang giay
#define maxCost				1000	// max cost excuting task, use in random create
#define minCost				13		// min cost //

#ifndef rateMutan
#define rateMutan			0.005	// rate of mutation is 0.5%
#endif

#define blockSize			16		// threads per block
#define threadsPerBlock  	256
#define locateIndex			8		// index of locate, use to recombine
#define rateRe				8		// so hang xom cua 1 thread

#define TaskCore			0
#define TaskMips			1
#define TaskMemory			2
#define TaskStorage			3
#define TaskNetwork			4
#define TaskStartTime 		5
#define TaskDuration 		6

#define MachineMipsPerCore	0
#define MachineCore			1
#define MachinePower		2
#define MachineBasePower	3
#define MachineMemory		4
#define MachineStorage		5
#define MachineNetwork		6

#define MAX_FIELD 			8
#define ID					0
#define CORES 				1
#define MIPS 				2
#define MEMORY 				3
#define STORAGE 			4
#define NETWORK 			5
#define TASK_START_TIME 	6
#define DURATION 			7
#define MAX_POWER 			6
#define BASE_POWER 			7

#define COMMENT		";"
#define outfile		"nst.txt"

#define ALLOC 	alloc_if(1) free_if(0)
#define FREE 	alloc_if(0) free_if(1)
#define REUSE 	alloc_if(0) free_if(0)

//-----------------------------sort field------------------------------
#pragma offload_attribute(push, target(mic))
class keySortLess { 
public: 
	int key;
	int follow;
	bool operator<(const keySortLess &other) const { 
		return key < other.key;
	}
};
class keySortGreater { 
public: 
	float key;
	int   follow;
	bool operator<(const  keySortGreater & other) const { 
		return key > other.key;
	}
};
void sort_by_key_less(int *keyarray, int numsize, int *valuearray) {
	keySortLess *dataSort = new keySortLess[numsize];

    //#pragma omp parallel for simd
	for (int i = 0; i < numsize; i++) {
		dataSort[i].key    = keyarray[i];
		dataSort[i].follow = valuearray[i];
	}
	stable_sort(dataSort, dataSort + numsize);

	//#pragma omp parallel for simd
    for (int i = 0; i < numsize; i++ ) {
		keyarray[i]   = dataSort[i].key;
		valuearray[i] = dataSort[i].follow;
	}
	delete [] dataSort;
}
void sort_by_key_greater(float *keyarray, int numsize, int *valuearray) {
	keySortGreater *dataSort = new keySortGreater[numsize];
	
    //#pragma omp parallel for simd
    for (int i = 0; i < numsize; i++) {
		dataSort[i].key 	= keyarray[i];
		dataSort[i].follow 	= valuearray[i];
	}
	stable_sort( dataSort, dataSort + numsize);
	
    //#pragma omp parallel for simd
    for (int i = 0; i < numsize; i++) {
		keyarray[i]  	= dataSort[i].key ;
		valuearray[i] 	= dataSort[i].follow ;
	}
	delete [] dataSort;
}
#pragma offload_attribute(pop) 
//---------------------------------------------------------------------

class Tasks {
public :
	int * cost;			// time need to complete task
	int * start_time;	// 
	int * ncore;		// number of required cpu cores for this task
	int * nmips;		// number of required MIPS/core of this task
	int * duration;		// deadline of task 
	int   nTask;		// real number of tasks
	int * memory;
	int * storage;
	int * network;

	void host_init_Task();
public :
	Tasks ( ) {
		cost 		= NULL;
		start_time 	= NULL;
		ncore 		= NULL;
		nmips 		= NULL;
		duration 	= NULL;
		memory 		= NULL;
		storage 	= NULL;
		network 	= NULL;
		nTask 		= numTask;
		host_init_Task();
		cout<<"Tasks database was created\n";
	}
	~Tasks ( ) {
		delete [] cost;
		delete [] start_time;
		delete [] ncore;
		delete [] nmips;
		delete [] duration;
		delete [] memory;
		delete [] storage;
		delete [] network;
		cout<<"Tasks database was destroyed\n";
	}
	int getCores 	(int i) { return ncore[i];		}
	int getMips 	(int i) { return nmips[i];		}
	int getMemory 	(int i) { return memory[i];		}
	int getStorage 	(int i) { return storage[i];	}
	int getNetwork 	(int i) { return network[i];	}
	int getNumTask 	(	  ) { return nTask; 		}
	int getStartTime(int i) { return start_time[i];	}
	int getDuration (int i) { return duration[i];	}
	
	void showTask() {
		cout<<"Cost\tStart time\tCores\tMips\tDura\tMem\tSto\tNet\n";
		for (int i = 0; i < nTask; i++) {
			cout << i << "\t" << cost[i] << "\t" << start_time[i] << "\t" 
				 << ncore[i] << "\t" << nmips[i] << "\t" << duration[i] 
				 << "\t" << memory[i] << "\t" << storage[i] << "\t" << network[i] << "\n";
		}
	}
	void createTask(int, int, int, int, int, int, int, int);
};
class Machines {
public :
	int * mipsPerCore;	// maximum of MIPS in each cpu core of this machine
	int * cores;		// maximum of cpu cores in this machine
	int   nMachine;		// real number of machines
	int * maxPower;		// maximum power for each machine
	int * basePower;	// machine base power
	int * memory;		// machine memory
	int * storage;		// machine storage
	int * network;		// machine network bandwitdth
//	long * allocatedMips; // number of allocated MIPS for assignment of each task	
	
	void host_init_Machine();
public :
	Machines () {
		mipsPerCore = NULL;
		cores 		= NULL;
		maxPower 	= NULL;
		basePower 	= NULL;
		memory 		= NULL;
		storage 	= NULL;
		network 	= NULL;
		nMachine 	= numMachine;
		host_init_Machine();
		cout << "Machine database was created\n";
	}
	~Machines () {
		delete [] mipsPerCore;
		delete [] cores;
		delete [] maxPower;
		delete [] basePower;
		delete [] memory;
		delete [] storage;
		delete [] network;
		cout << "Machine database was destroyed\n";
	}
	
	void createMachine(int, int, int, int, int, int, int, int);
	void showMachine() {
		cout<<"Mip/C\tCore\tMax\tpwr\tBs\tpwr\tMem\tSto\tNet\tMPP\n";
		for (int i = 0; i < nMachine; i++) {
			cout << i << "\t" << mipsPerCore[i] << "\t" << cores[i] << "\t" 
				 << maxPower[i] << "\t" << basePower[i] << "\t" << memory[i] 
				 << "\t" << storage[i] << "\t" << network[i] << "\t" 
				 << (float) mipsPerCore[i]*cores[i] / maxPower[i] <<"\n";
		}
	}
	int getNumMachine 	(	  ) { return nMachine;		}
	int getCores 		(int i) { return cores[i];		}
	int getMipsPerCore 	(int i) { return mipsPerCore[i];}
	int getMemory 		(int i) { return memory[i];		}
	int getStorage 		(int i) { return storage[i];	}
	int getNetwork 		(int i) { return network[i];	}
	int getMaxPower 	(int i) { return maxPower[i];	}
	int getBasePower 	(int i) { return basePower[i];	}

};

class Population {
private:
	int 	 * mang_NST;		// set of NSTs
	float 	 * mang_fitness;	// do tuong thich
	int 	 * mang_makespan;	// thoi gian hoan thanh solution
	int 	   size_pop;		// so luong NST
	int 	 * rankings;		// rankings[0] = 5 meaning  NST 5 is the best solution.
	Tasks  	 * dataT;			// du lieu cac Tasks
	Machines * dataM;			// du lieu cac Machines

	void host_init_NST();
	void parseValue(string, bool);
	void extractField(float*, bool);

    __attribute__ ((target (mic))) void cross_over(int*, int*);    
    __attribute__ ((target (mic))) void evalNST(int *, float*, int*, int*, int*, int*, int*, int*, int*, int*);
	__attribute__ ((target (mic))) void selection(int*, float*, int*, float*, int*);
	__attribute__ ((target (mic))) void mutan(int*);
	
public:
	Population() {
		mang_NST 		= NULL;
		mang_fitness 	= NULL;
		mang_makespan 	= NULL;
		rankings 		= NULL;
		dataT 			= new Tasks();
		dataM 			= new Machines();
		size_pop 		= sizePop;
		host_init_NST();
		cout<<"The population was created\n";
	}
	~Population() {
		delete [] mang_NST;
		delete [] mang_fitness;
		delete [] mang_makespan;
		delete [] rankings;
		delete 	  dataT;
		delete 	  dataM;
		cout << "The population has been destroyed\n";
	}
	void GA_Evolution (int);
	void showNST() {
		for (int i = 300; i < 310; i++) {
			cout << "\nNhiem sac the thu " << i << "\n";
			for (int j = 0; j < numTask; j++) {
				cout << mang_NST[i*numTask + j] << "\t";
			}
		}
		cout<< "\n10 NST dau tien\n";
	}
	void printBest(char* ten_file) {
		// print to screen
/*
		cout<<"Simulation Results :\n";
		cout<<"TaskID\tHostID\tTime\tFinish\tCores\n";
		for (int i = 0; i < numTask; i++) {
			cout << i << "\t" << mang_NST[i] << "\t" << dataT->getStartTime(i) <<"\t" 
				 << dataT->getStartTime(i) + dataT->getDuration(i) << "\t" << dataT->getCores(i) << endl;
		}
*/
		cout << "Total energy consumption (Watt-hour): " << 1/mang_fitness[0] << endl;
		cout << "Fitness of ga solution: " << mang_fitness[0] << endl;

		// write to output file
		ofstream myfile(ten_file, ios::out);
		if (!myfile) {
			cout << "Error : Can't write output file best solution ";
			exit(EXIT_FAILURE);
		}

		myfile << "Fitness of ga solution: " << mang_fitness[0] << " - ";
		for( int i = 0; i < numTask; i++) {
			myfile << mang_NST[i] << " ";
		}
		myfile << endl;
		myfile.close();

		cout << "Writed best solution to file done" << endl;
	}

	void xuat_file(char* ten_file) {
		ofstream myfile(ten_file, ios::out);
		if (!myfile) {
			cout << "Error : Can't write file " << ten_file;
			exit(EXIT_FAILURE);
		}
		for (int i = 0; i < size_pop; i++) {
			myfile << mang_fitness[i] << " - ";
			for(int j = 0; j < numTask-1; j++) {
				myfile << mang_NST[i*numTask +j] << " ";
			}
			myfile << mang_NST[(i+1)*numTask -1] << "\n";
		}
		myfile.close();
	}
	void printNST(char*, int*, float);
	void bestFit();
	void readFile(char *);
};

/*
* random number in [s,d]
*/
#pragma offload_attribute(push, target(mic))
#pragma omp declare simd
int my_host_rand (int s, int d) {
	return rand() % (d-s+1) + s;
}
int my_parallel_rand (int s, int d) {
    return (lrand48()) % (d-s+1) + s;
}
#pragma offload_attribute(pop) 

void Population::host_init_NST () {

	this->mang_fitness 	= new float[sizePop];
	this->mang_makespan = new int[sizePop];
	this->rankings 		= new int[2*sizePop];
	this->mang_NST 		= new int[numTask*sizePop];
// Allocate the host input array 
	
	if (this->mang_NST == NULL )
    {
        fprintf(stderr, "Failed to allocate host vectors!\n");
        exit(EXIT_FAILURE);
    }

    cout << "Begin init random NSTs\n";
    int *mang_NST = this->mang_NST; 
    #pragma offload target(mic) out(mang_NST : length(numTask*sizePop))
    {   
        #pragma omp parallel for simd
        for (int i = 0; i < numTask*sizePop; i++) { 
            mang_NST[i] = my_parallel_rand(0, numMachine-1);
        }   
    }
    
    /* Test random
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < numTask; i++) {
            cout << mang_NST[i + j*numTask] << ' ';
        }
        cout << "\nXXX\n";
    }
    cout << endl;
    */
    cout << "Creating NSTs done\n";
    
    /*  Test result NST
    for (int i = 0; i < 10; i++ ) {
        cout << mang_NST[i] << endl ;
    }
    */                      
    
}

void Population::cross_over(int *NSToffspring, int *mang_NST) {
    #pragma omp parallel for
	for (int i = 0; i < sizePop; i++){
		int cross_point = my_parallel_rand(0, numTask - numTask/2);
		int neighbor 	= my_parallel_rand(0, sizePop - 1);
		if (neighbor == i) {neighbor = (i+1)%sizePop;}
        #pragma omp simd
		for (int p = 0; p < numTask/2 - 1; p++) {
			NSToffspring[i*numTask + cross_point + p] = mang_NST[neighbor*numTask + cross_point +p];
		}
	}
}

void Population::mutan (int* NSToffspring) {
    #pragma omp parallel for simd // speedup from report < 1
	for (int i = 0; i < numTask*sizePop; i++) {
		if (((float)lrand48() / RAND_MAX) < rateMutan) {
			NSToffspring[i] = my_parallel_rand(0, numMachine-1);
		}
	}
}

void Population::selection(int *NSToffspring, float *fitness_offspring, int *mang_NST, float *mang_fitness, int *rankings) {
// reset rankings
    #pragma omp parallel for simd
	for (int i = 0; i < 2*sizePop; i++) {rankings[i] = i;}
// copy key sort 
	float keySort[2*sizePop];
	#pragma omp parallel for simd
    for (int i = 0; i < sizePop; i++) {
		keySort[i] 			= mang_fitness[i]; 
		keySort[i+sizePop] 	= fitness_offspring[i]; 
	}
// sort rankings[] by key
	sort_by_key_greater(keySort, 2*sizePop, rankings);
// move top NST
	int temNST[sizePop*numTask];

    #pragma omp parallel for
	for (int i = 0; i < sizePop; i++) {
		// copy to temNST
		if (rankings[i] < sizePop) {
            #pragma omp simd
			for (int p = 0; p < numTask; p++) {
				temNST[i * numTask + p] = mang_NST[rankings[i]*numTask +p];
			}
		} else {
            #pragma omp simd
			for (int p = 0; p < numTask; p++) {
				temNST[i * numTask + p] = NSToffspring[(rankings[i] - sizePop) * numTask + p];
			}	
		}
	}
	// copy to parentsNST
	#pragma omp parallel for simd
    for (int i = 0; i < sizePop*numTask; i++) {
		mang_NST[i] 	= temNST[i]; 
		NSToffspring[i] = temNST[i];
	}
	// copy to parents fitness
	#pragma omp parallel for simd
	for (int i = 0; i < sizePop; i++) {
		mang_fitness[i] 	 = keySort[i];
		fitness_offspring[i] = keySort[i];
	}
}

void Population::GA_Evolution(int galoop) {
	if (galoop < 1 )
    {
        fprintf(stderr, "so the he khong hop le\n");
        exit(EXIT_FAILURE);
    }
// bring pointers to local to offload to MIC
    int   	*mang_NST 		= this->mang_NST;
    float 	*mang_fitness 	= this->mang_fitness;
    int 	*rankings 		= this->rankings;
    int 	*ncore 			= dataT->ncore;
    int 	*nmips 			= dataT->nmips;
    int 	*mipsPerCore 	= dataM->mipsPerCore;
    int 	*cores 			= dataM->cores;
    int 	*maxPower 		= dataM->maxPower;
    int 	*basePower 		= dataM->basePower;
// offspring
	int 	NSToffspring[numTask*sizePop];
	float 	fitness_offspring[sizePop];
	// create clone
	#pragma omp parallel for simd // run slower
	for (int i = 0; i < sizePop*numTask; i++) {NSToffspring[i] = mang_NST[i];}
// event list , create for eval
	cout << "Begin create list\n";
	int startList[numTask * 2], endList[numTask * 2];
    #pragma omp parallel for simd
	for(int i = 0; i < numTask; i++){
		startList[i] 		 = dataT->getStartTime(i);
		endList[i + numTask] = startList[i + numTask] = i;
		endList[i] 			 = startList[i] + dataT->getDuration(i);
	}
	sort_by_key_less(startList, numTask, startList + numTask);
	sort_by_key_less(endList, numTask, endList + numTask);
	cout << "Event list were created" << endl;

    // call eval funtion for parent
    #pragma offload target(mic) in(mang_NST : length(numTask*sizePop) ALLOC)\
    inout(mang_fitness : length(sizePop) ALLOC) \
    in(rankings : length(2*sizePop) ALLOC) \
    in(startList ) in(endList ) \
    in(ncore : length(numTask) ALLOC) \
    in(nmips : length(numTask) ALLOC) \
    in(mipsPerCore : length(numMachine) ALLOC) \
    in(cores : length(numMachine) ALLOC) \
    in(maxPower : length(numMachine) ALLOC) \
    in(basePower : length(numMachine) ALLOC)
    {
	    evalNST(mang_NST, mang_fitness, startList, endList, ncore, nmips, mipsPerCore, cores, maxPower, basePower);
    }
    // clone fitness
    #pragma omp parallel for simd // slower
	for (int i = 0; i< sizePop; i++) {fitness_offspring[i] = mang_fitness[i];}
// begin ga 
	clock_t t;
	//t = clock();
    #pragma offload target(mic) nocopy(mang_NST : length(numTask*sizePop) REUSE) \
    nocopy(mang_fitness : length(sizePop) REUSE)\
    nocopy(rankings : length(2*sizePop) REUSE) \
    in(startList) in(endList) in(NSToffspring) in(fitness_offspring)\
    nocopy(ncore : length(numTask) REUSE) \
    nocopy(nmips : length(numTask) REUSE) \
    nocopy(mipsPerCore : length(numMachine) REUSE) \
    nocopy(cores : length(numMachine) REUSE) \
    nocopy(maxPower : length(numMachine) REUSE) \
    nocopy(basePower : length(numMachine) REUSE)
    { 
        for (int g = 0; g < galoop; g++) {
            // call lai_ghep
            cross_over(NSToffspring, mang_NST);
        
            // call eval for offspring
            evalNST(NSToffspring, fitness_offspring, startList, endList, ncore, nmips, mipsPerCore, cores, maxPower, basePower);
            // call chon loc sort
            selection(NSToffspring, fitness_offspring, mang_NST, mang_fitness, rankings);
            // call dot bien
            mutan(NSToffspring);
            // call eval for offspring
            evalNST(NSToffspring, fitness_offspring, startList, endList, ncore, nmips, mipsPerCore, cores, maxPower, basePower);
            // call chon loc sort
            selection(NSToffspring, fitness_offspring, mang_NST, mang_fitness, rankings);
        //	t = clock() - t;
        }
    }
    #pragma offload target(mic) out(mang_NST : length(numTask*sizePop) FREE) \
    out(mang_fitness : length(sizePop) FREE)\
    nocopy(rankings : length(2*sizePop) FREE) \
    nocopy(ncore : length(numTask) FREE) \
    nocopy(nmips : length(numTask) FREE) \
    nocopy(mipsPerCore : length(numMachine) FREE) \
    nocopy(cores : length(numMachine) FREE) \
    nocopy(maxPower : length(numMachine) FREE) \
    nocopy(basePower : length(numMachine) FREE)
    {}
}

void Machines::host_init_Machine() {
	mipsPerCore = new int[nMachine];
	cores 		= new int[nMachine];
	maxPower 	= new int[nMachine];
	basePower 	= new int[nMachine];
	memory 		= new int[nMachine];
	storage 	= new int[nMachine];
	network 	= new int[nMachine];
    #pragma omp parallel for simd
	for(int i = 0; i < nMachine; i++){
		mipsPerCore[i] 	= 1;
		cores[i] 		= 1;
		maxPower[i] 	= 1;
		basePower[i] 	= 1;
		memory[i] 		= 1;
		storage[i] 		= 1;
		network[i] 		= 1;
	}
}

void Tasks::createTask(int id, int cores, int mips, int mem, int strg, int net, int startTime, int dur){
	cost[id] 		=	0;
	duration[id] 	= dur;
	start_time[id] 	= startTime;
	ncore[id] 		= cores;
	nmips[id] 		= mips;
	memory[id] 		= mem;		// my_host_rand(1, 4096); 	// Unit in MBytes
	storage[id] 	= strg; 	// my_host_rand(50, 5000);
	network[id] 	= net;		// my_host_rand(1, 10000);	//1 Mb/s - 10000 Mb/s
}

void Tasks::host_init_Task() {
	cost 		= new int[nTask];
	start_time 	= new int[nTask];
	ncore 		= new int[nTask];
	nmips 		= new int[nTask];
	duration 	= new int[nTask];
	memory 		= new int[nTask];	
	storage 	= new int[nTask];
	network 	= new int[nTask];
    #pragma omp parallel for simd
	for(int i = 0; i < nTask; i++){
		cost[i] 		= 0;
		start_time[i] 	= 0;
		ncore[i] 		= 0;
		nmips[i] 		= 0;
		duration[i] 	= 0;
		memory[i] 		= 0;
		storage[i] 		= 0;
		network[i] 		= 0;
	}
}


/*
Function name: bestFit

Function description: find a solution using best fit algorithm and write out into bestfit.txt file

Input: N/A
*/
void Population::bestFit(){
	//Code added: Create 2 list eventStartList and eventEndList
	int startList[numTask * 2], endList[numTask * 2];
	#pragma omp parallel for
    for (int i = 0; i < numTask; i++){
		startList[i] 		 = dataT->getStartTime(i);
		endList[i + numTask] = startList[i + numTask] = i;
		endList[i] 			 = startList[i] + dataT->getDuration(i);
	}

	sort_by_key_less(startList, numTask, startList + numTask);
	sort_by_key_less(endList, numTask, endList + numTask);

	int mips[numMachine], cores[numMachine], bestMachine;		//mpp = mips per power
	float bestmpp, mpp;

	//Initialize variables
	#pragma omp parallel for simd
	for(int i = 0; i < numMachine; i++){
		mips[i] = cores[i] = 0;
	}

	int startPtr = 0, endPtr = 0, task;

	while (endPtr != numTask){
		bestmpp 	= -1;
		bestMachine = -1;
		
		//case 1: started event begin sooner than ended event
		if(startPtr != numTask && startList[startPtr] <= endList[endPtr]){
			task = startPtr;
			startPtr++;
			for (int j = 0; j < numMachine; j++){
			//Checking for invalid assignment
				if(dataM->getMemory(j) < dataT->getMemory(task)
					|| dataM->getStorage(j) < dataT->getStorage(task)
					|| dataM->getNetwork(j) < dataT->getNetwork(task)
					|| dataM->getMipsPerCore(j) < dataT->getMips(task)
					|| dataM->getMipsPerCore(j) * dataM->getCores(j) < mips[j] + dataT->getMips(task) * dataT->getCores(task)
					|| dataM->getCores(j) < cores[j] + dataT->getCores(task)
					){
					continue;
				}
				//compare the machine with previous machine to choose out the best machine which has minimum MIPS/power
				mpp = (float) dataM->getMipsPerCore(j) * dataM->getCores(j) / dataM->getMaxPower(j);
				if (mpp > bestmpp){
					bestmpp 	= mpp;
					bestMachine = j;
				}
			}
			//assign the best machine to task
			mang_NST[task] = bestMachine;
			if (bestMachine < 0) break;
			mips[bestMachine] 	+= dataT->getMips(task) * dataT->getCores(task);
			cores[bestMachine] 	+= dataT->getCores(task);
		}
		//case 2: started event begin later than ended event
		else if (startPtr == numTask || endList[endPtr] < startList[startPtr]){
			task = endPtr;
			endPtr++;
			bestMachine = mang_NST[task];
			mips[bestMachine] 	-= dataT->getMips(task) * dataT->getCores(task);
			cores[bestMachine] 	-= dataT->getCores(task);
		}
	}

	bool skip = 0;

	for(int i = 0; i < numTask; i++){
		if(mang_NST[i] < 0){
			skip = 1;
			break;
		}
	}

	float fitness = 0;
	//Calculate fitness
	if (skip)
		fitness = -1;
	else{
		float power[numMachine], powerDatacenter = 0;
		#pragma omp parallel for simd
        for(int i = 0; i < numMachine; i++){
			power[i] = 0;
		}
        
		for(int i = 0; i < numTask; i++){
			power[mang_NST[i]] += (float)dataT->getMips(i) * dataT->getCores(i) 
								/ (dataM->getMipsPerCore(mang_NST[i]) * dataM->getCores(mang_NST[i])) 
								* (dataM->getMaxPower(mang_NST[i]) - dataM->getBasePower(mang_NST[i])) 
								* dataT->getDuration(i);
		}
		for(int j = 0; j < numMachine; j++){
			if (power[j] > 0)
				powerDatacenter += dataM->getBasePower(j);
			powerDatacenter += power[j];
		}
		fitness = (float) 1/powerDatacenter;
	}
	
	//Print out to file bestfit.txt
	ofstream myfile("bestfit.txt", ios::out);
		if (!myfile) {
			cout << "Error : Can't write file bestfit.txt";
			exit(EXIT_FAILURE);
		}
        cout << "Best fit solution: " << (float) fitness << endl;
		myfile << fitness << " - ";
		for(int j = 0; j < numTask; j++) {
			myfile << mang_NST[j] << " ";
		}
		myfile << endl;			
	myfile.close();
}

void Machines::createMachine(int id, int core, int mips, int mem, int strg, int net, int maxPwr, int basePwr){
	maxPower[id] 	= maxPwr;
	cores[id] 		= core;
	mipsPerCore[id] = mips;
	memory[id] 		= mem;
	storage[id] 	= strg;
	network[id] 	= net;
	basePower[id] 	= basePwr;
}

void Population::printNST(char *ten_file, int *nst, float fitness) {
		ofstream myfile(ten_file, ios::out);
		if (!myfile) {
			cout<< "Error : Can't write file nst.txt";
			exit(EXIT_FAILURE);
		}
		myfile << fitness <<" - ";
		for (int i = 0; i < numTask; i++) {
			myfile << nst[i] << " ";
		}
		myfile.close();
}
/*
void Population::evalNST(int *nst, float* fitness, int* eventStartList, int* eventEndList, 
                                    int *ncore, int *nmips, int * mipsPerCore, int * cores, int * maxPower, int * basePower){
*/
void Population::evalNST(int *nst, float* fitness, int* eventStartList, int* eventEndList, int* taskCore, int* taskMips, 
						 int* machineMipsPerCore, int* machineCore, int* machineMaxPower, int* machineBasePower){
	#pragma omp parallel for
    for (int i = 0; i < sizePop; i++ ) {
		float power[numMachine], powerDatacenter = 0;
		int mips[numMachine], cores[numMachine], times[numMachine];
        #pragma omp parallel for simd
		for(int j = 0; j < numMachine; j++){
			power[j] = 0;
			mips[j]  = 0;
			cores[j] = 0;
		}
		fitness[i] = 0;

		int startPtr = 0, endPtr = 0;		//startPtr and endPtr is 2 pointers for check the lists
		//endPtr is criteria because ended event always happen after started event
		while(endPtr != numTask){

			/*
			eventStartList[startPtr]: time start event
			eventStartList[startPtr+numTask]: task ID
			nst[i*numTask + eventStartList[startPtr+numTask]]: machine ID
			*/

			//case 1: started event begin sooner than ended event
			if(startPtr != numTask && eventStartList[startPtr] <= eventEndList[endPtr]){
				//Check condition, temp1 is mips after add new task, temp2 is total machine MIPS, temp3 is cores after add new task
				int temp1 = mips[nst[i*numTask + eventStartList[startPtr+numTask]]] 
						  	+ taskMips[eventStartList[startPtr+numTask]] * taskCore[eventStartList[startPtr+numTask]];
				int temp2 = machineMipsPerCore[nst[i*numTask+eventStartList[startPtr+numTask]]] 
						  	* machineCore[nst[i*numTask+eventStartList[startPtr+numTask]]];
				int temp3 = cores[nst[i*numTask + eventStartList[startPtr+numTask]]] 
							+ taskCore[eventStartList[startPtr+numTask]];
				if (temp1 > temp2 || temp3 > machineCore[nst[i*numTask+eventStartList[startPtr+numTask]]] 
					|| taskMips[eventStartList[startPtr+numTask]] > machineMipsPerCore[nst[i*numTask+eventStartList[startPtr+numTask]]]){
					if(fitness[i] > 0) fitness[i] = 0;
					fitness[i]--;
					startPtr++;
					continue;
				}
				
				power[nst[i*numTask+eventStartList[startPtr+numTask]]]  += (float)mips[nst[i*numTask + eventStartList[startPtr+numTask]]] / temp2 
																			* (machineMaxPower[nst[i*numTask+eventStartList[startPtr+numTask]]] - machineBasePower[nst[i*numTask+eventStartList[startPtr+numTask]]]) 
																			* (eventStartList[startPtr] - times[nst[i*numTask + eventStartList[startPtr+numTask]]]);
				times[nst[i*numTask+eventStartList[startPtr+numTask]]] 	 = eventStartList[startPtr];
				mips[nst[i*numTask + eventStartList[startPtr+numTask]]]  = temp1;
				cores[nst[i*numTask + eventStartList[startPtr+numTask]]] = temp3;
				startPtr++;
			}

			//case 2: started event begin later than ended event
			else if(startPtr == numTask || eventEndList[endPtr] < eventStartList[startPtr]){
				//Check condition, temp1 is mips after remove old task, temp2 is total machine MIPS, temp3 is cores after remove old task
				int temp1 = mips[nst[i*numTask + eventEndList[endPtr+numTask]]] 
							- taskMips[eventEndList[endPtr+numTask]] * taskCore[eventEndList[endPtr+numTask]];
				int temp2 = machineMipsPerCore[nst[i*numTask+eventEndList[endPtr+numTask]]] 
							* machineCore[nst[i*numTask+eventEndList[endPtr+numTask]]];
				int temp3 = cores[nst[i*numTask + eventEndList[endPtr+numTask]]] 
							- taskCore[eventEndList[endPtr+numTask]];

				power[nst[i*numTask+eventEndList[endPtr+numTask]]] 	+= (float)mips[nst[i*numTask + eventEndList[endPtr+numTask]]] / temp2 
																		* (machineMaxPower[nst[i*numTask+eventEndList[endPtr+numTask]]] - machineBasePower[nst[i*numTask+eventEndList[endPtr+numTask]]]) 
																		* (eventEndList[endPtr] - times[nst[i*numTask + eventEndList[endPtr+numTask]]]);
				times[nst[i*numTask+eventEndList[endPtr+numTask]]] 	 = eventEndList[endPtr];
				mips[nst[i*numTask + eventEndList[endPtr+numTask]]]  = temp1;
				cores[nst[i*numTask + eventEndList[endPtr+numTask]]] = temp3;
				endPtr++;
			}

		}

		if (fitness[i] >= 0){
			for (int j=0; j< numMachine; j++){
				if (power[j] > 0)
					powerDatacenter += machineBasePower[j];
				powerDatacenter += power[j];
			}
			fitness[i] = 1.0/powerDatacenter;
            //cout << "F " << fitness[i] << endl;
		}
	}
}


void Population::readFile (char *fname){
	cout << "Begin read file\n";
	string line;

	ifstream infile(fname, ifstream::in);
	if (!infile.is_open()){
		cout << "File open error\n";
		return;
	}
	cout << "File is opened\n";

	bool type = 0;
	getline(infile, line);
	
	int ntasks = 0, nmachines = 0;
	//Edit code from here 	xet tao task hay machine dua vao bien type
	while (!infile.eof()){
		if (line[0] != ';' && line[0] != NULL){
			if (line[0] == 'T'){
				type = 0;
			}
			else if (line[0] == 'M'){
				type = 1;
			}
			else if (type == 0 && ntasks < numTask){
				parseValue(line, type);
				ntasks++;
			}
			else if (type == 1 && nmachines < numMachine){
				parseValue(line, type);
				nmachines++;
			}
		}
		if (ntasks == numTask && nmachines == numMachine)
			break;
		getline(infile,line);
	}
	//End
	
	infile.close();
	//dataT->showTask();
	//dataM->showMachine();
}

void Population::parseValue(string line, bool type){
	//cout << "Begin parse value\n";
	float field[MAX_FIELD];

	//filter out space and extract field
	int index = 0;
	stringstream ss(line);
	while (ss.peek() != EOF){
		ss >> field[index];
		index++;
	}
	extractField(field, type);
}

void Population::extractField(float *fieldArray, bool type){
	//cout << "Begin extract field\n";
	//Create variable ready for task
	int id 			= fieldArray[ID];
	int cores 		= fieldArray[CORES];
	int mipsPerCore = fieldArray[MIPS];
	int memory 		= fieldArray[MEMORY];
	int storage 	= fieldArray[STORAGE];
	int network 	= fieldArray[NETWORK];

	if(type == 0){
		int startTime = fieldArray[TASK_START_TIME], duration = fieldArray[DURATION];
		dataT->createTask(id, cores, mipsPerCore, memory, storage, network, startTime, duration);
	}
	else{
		int maxPwr = fieldArray[MAX_POWER];
		int basePwr = fieldArray[BASE_POWER];
		dataM->createMachine(id, cores, mipsPerCore, memory, storage, network, maxPwr, basePwr);
	}
}

/*
int main() {
	srand(time(NULL));
	Population P;
    
	P.readFile("data.txt");
	P.bestFit();

	P.xuat_file("nst.txt");
	
	struct timeval timeRunOfGa;  
    gettimeofday(&timeRunOfGa, NULL);  
    double dTime1 = timeRunOfGa.tv_sec+(timeRunOfGa.tv_usec/1000000.0);  
     
	P.GA_Evolution(generation);
	
	gettimeofday(&timeRunOfGa, NULL);  
    double dTime2 = timeRunOfGa.tv_sec+(timeRunOfGa.tv_usec/1000000.0); 
	
	P.xuat_file("after_ga.txt");
	P.printBest("best_solution.txt");
	cout<<endl<<"total time run ga: "<< dTime2 - dTime1 << endl<<endl;
	cout<<"Done\n";
	
	return 0;
}
*/

int main() {
    srand(time(NULL));
    srand48(time(NULL));
#pragma omp parallel 
{   
    if (omp_get_thread_num() == 0)
        cout << "Num threads " << omp_get_num_threads() << endl;
}
    double avarage_time = 0;
    int i, temp = 0;
for (i = 0; i < 5; i++) {
    Population P;
    P.readFile("data.txt");

    struct timeval timeRunOfGa;  
    gettimeofday(&timeRunOfGa, NULL);  
    double dTime1 = timeRunOfGa.tv_sec+(timeRunOfGa.tv_usec/1000000.0);  
     
    P.GA_Evolution(generation);
    
    gettimeofday(&timeRunOfGa, NULL);  
    double dTime2 = timeRunOfGa.tv_sec+(timeRunOfGa.tv_usec/1000000.0); 
    
    P.xuat_file("after_ga.txt");
    P.printBest("best_solution.txt");
    //temp = P.mang_fitness[0];
    P.bestFit();
    P.xuat_file("nst.txt");
    cout<<endl<<"total time run ga: "<< dTime2 - dTime1 << endl<<endl;
    cout << "Test parameters: machines " << numMachine << " tasks " << numTask 
         << " population size " << sizePop << " generations " << generation  
         << " Mutation rate: " << rateMutan << endl;
    
    cout<<"Done\n";
    avarage_time += dTime2 - dTime1;
}
    cout << "Avarage time: " << avarage_time/(i) << endl;
    return 0;
}
