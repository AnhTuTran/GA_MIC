#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <string>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

using namespace std;

// Define
#ifndef numTask
#define numTask                 500         // number of Tasks, task 0, 1, ...
#endif
#ifndef numMachine
#define numMachine              500         // number of Machines, machine 0, 1, ...,       criteria: <=1000
#endif
#ifndef sizePop
#define sizePop                 1024        // size of population
#endif
#define edgeSizePop             64          // edge of population
#ifndef generation
#define generation              1000        // number of generations
#endif
#define timeStopGa              3600        // thoi gian toi da chay giai thuat ga, tinh bang giay
#define maxCost                 1000        // max cost excuting task, use in random create
#define minCost                 13          // min cost //
#ifndef rateMut
#define rateMut                 0.0155      // rate of mutation is 0.5%
#endif
#ifndef threadsPerBlock
#define threadsPerBlock         128
#endif
#define locateIndex             8           // index of locate, use to recombine
#define rateRe                  8           // so hang xom cua 1 thread

#define TaskCore                0
#define TaskMips                1
#define TaskMemory              2
#define TaskStorage             3
#define TaskNetwork             4
#define TaskStartTime           5
#define TaskDuration            6

#define MachineMipsPerCore      0
#define MachineCore             1
#define MachinePower            2
#define MachineBasePower        3
#define MachineMemory           4
#define MachineStorage          5
#define MachineNetwork          6

#define MAX_FIELD               8
#define ID                      0
#define CORES                   1
#define MIPS                    2
#define MEMORY                  3
#define STORAGE                 4
#define NETWORK                 5
#define TASK_START_TIME         6
#define DURATION                7
#define MAX_POWER               6
#define BASE_POWER              7

#define COMMENT                 ";"

#define outfile                 "nst.txt"

#define CUDA_CALL(x) do { if((x) != cudaSuccess) { \
    printf("Error at %s:%d\nError code %s\n",__FILE__,__LINE__, cudaGetErrorString(x)); \
    exit(EXIT_FAILURE);}} while(0)

class Tasks {
private :
    int * cost;         // time need to complete task
    int * start_time;   // 
    int * ncore;        // number of required cpu cores for this task
    int * nmips;        // number of required MIPS/core of this task
    int * duration;     // deadline of task 
    int   nTask;        // real number of tasks
    int * memory;
    int * storage;
    int * network;

    void host_init_Task();
public :
    Tasks ( ) {
        cost        = NULL;
        start_time  = NULL;
        ncore       = NULL;
        nmips       = NULL;
        duration    = NULL;
        memory      = NULL;
        storage     = NULL;
        network     = NULL;
        nTask       = numTask;
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
    int getCores    (int i) { return ncore[i];      }
    int getMips     (int i) { return nmips[i];      }
    int getMemory   (int i) { return memory[i];     }
    int getStorage  (int i) { return storage[i];    }
    int getNetwork  (int i) { return network[i];    }
    int getNumTask  (     ) { return nTask;         }
    int getStartTime(int i) { return start_time[i]; }
    int getDuration (int i) { return duration[i];   }
    
    void showTask() {
        cout<<" Cost    Start time  Cores   Mips    Dura    Mem Sto Net\n";
        for ( int i = 0; i<nTask; i++) {
            cout << i <<"  "<< cost[i]<< "   " << start_time[i]
                 << "     "<< ncore[i]<< "  " << nmips[i]<< "  " 
                 << duration[i] << "   " << memory[i]<< " " 
                 << storage[i] << "    "<< network[i] << "\n";
        }
    }

    int* importFromTask(char tenMang){
        int* copyVar = NULL;
        int* d_var   = NULL;
        switch (tenMang){
            case TaskCore:
                copyVar = ncore;
                break;
            case TaskMips:
                copyVar = nmips;
                break;
            case TaskMemory:
                copyVar = memory;
                break;
            case TaskStorage:
                copyVar = storage;
                break;
            case TaskNetwork:
                copyVar = network;
                break;
            case TaskStartTime:
                copyVar = start_time;
                break;
            case TaskDuration:
                copyVar = duration;
                break;
            default:
                return NULL;
        }
        CUDA_CALL( cudaMalloc(&d_var, sizeof(int) * nTask));
        CUDA_CALL( cudaMemcpy(d_var, copyVar, sizeof(int) * nTask, cudaMemcpyHostToDevice));
        return d_var;
    }
    void createTask(int, int, int, int, int, int, int, int);
};

class Machines {
private :
    int  * mipsPerCore;     // maximum of MIPS in each cpu core of this machine
    int  * cores;           // maximum of cpu cores in this machine
    int    nMachine;        // real number of machines
    int  * maxPower;        // maximum power for each machine
    int  * basePower;       // machine base power
    int  * memory;          // machine memory
    int  * storage;         // machine storage
    int  * network;         // machine network bandwitdth
    long * allocatedMips;   // number of allocated MIPS for assignment of each task   
    
    void host_init_Machine();
public :
    Machines () {
        mipsPerCore = NULL;
        cores       = NULL;
        maxPower    = NULL;
        basePower   = NULL;
        memory      = NULL;
        storage     = NULL;
        network     = NULL;
        nMachine    = numMachine;
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
        cout<<" Mip/C   Core    Max pwr Bs pwr  Mem Sto Net MPP\n";
        for (int i = 0; i< nMachine; i++) {
            cout << i << "  " << mipsPerCore[i] << "\t" << cores[i] << "\t" << maxPower[i] 
                 << "\t"<< basePower[i] << "\t" << memory[i] << "\t" << storage [i] << "\t" 
                 << network[i] << "\t" << (float)mipsPerCore[i]*cores[i]/maxPower[i] << "\n";
        }
    }
    int getNumMachine() { 
        return nMachine;
    }
    int getCores(int i){
        return cores[i];
    }
    int getMipsPerCore(int i){
        return mipsPerCore[i];
    }
    int getMemory(int i){
        return memory[i];
    }
    int getStorage(int i){
        return storage[i];
    }
    int getNetwork(int i){
        return network[i];
    }
    int getMaxPower(int i){
        return maxPower[i];
    }
    int getBasePower(int i){
        return basePower[i];
    }
    int* importFromMachine(char tenMang){
        int* copyVar = NULL;
        int* d_var   = NULL;
        switch (tenMang){
            case MachineMipsPerCore:
                copyVar = mipsPerCore;
                break;
            case MachineCore:
                copyVar = cores;
                break;
            case MachinePower:
                copyVar = maxPower;
                break;
            case MachineBasePower:
                copyVar = basePower;
                break;
            case MachineMemory:
                copyVar = memory;
                break;
            case MachineStorage:
                copyVar = storage;
                break;
            case MachineNetwork:
                copyVar = network;
                break;
            default:
                return NULL;
        }
        CUDA_CALL( cudaMalloc(&d_var, sizeof(int) * nMachine));
        CUDA_CALL( cudaMemcpy(d_var, copyVar, sizeof(int) * nMachine, cudaMemcpyHostToDevice));
        return d_var;
    }

};

class Population {
private :
    int         * mang_NST;         // set of NSTs
    float       * mang_fitness;     // do tuong thich
    int         * mang_makespan;    // thoi gian hoan thanh solution
    int           size_pop;         // so luong NST
    int         * rankings;         // rankings[0] = 5 meaning  NST 5 is the best solution.
    Tasks       * dataT;            // du lieu cac Tasks
    Machines    * dataM;            // du lieu cac Machines

    void  host_init_NST();
    void  vetcan(int[], int, float*);
    float evalNST(int*);
    void  parseValue(string, bool);
    void  extractField(float*, bool);
    void  sort_Pop( float*, float*);

public:
    Population() {
        mang_NST        = NULL;
        mang_fitness    = NULL;
        mang_makespan   = NULL;
        rankings        = NULL;
        dataT           = new Tasks();
        dataM           = new Machines();
        size_pop        = sizePop;
        host_init_NST();
        cout << "The population was created\n";
    }
    ~Population() {
        delete [] mang_NST;
        delete [] mang_fitness;
        delete [] mang_makespan;
        delete [] rankings;
        delete    dataT;
        delete    dataM;
        cout << "The population has been destroyed\n";
    }
    void GA_Evolution (int);
    void showNST() {
        for (int i = 300; i < 310; i++) {
            cout << "\nNhiem sac the thu "<< i << "\n";
            for (int j =0; j < numTask; j++ ) {
                cout<< mang_NST[i*numTask + j] << "  ";
            }
        }
        cout << "\n10 NST dau tien\n";
    }
    void printBest(char* ten_file) {
        // print to screen
        cout << "Simulation Results :\n";
        cout << "TaskID   HostID  Time    Finish  Cores\n";
        for (int i = 0; i< numTask ; i++ ) {
            cout << i << "  " << mang_NST[i] << "  " << dataT->getStartTime(i) << "   " 
                 << dataT->getStartTime(i) + dataT->getDuration(i) << "  " 
                 << dataT->getCores(i) << endl;
        }
        cout << "Total energy consumption (Watt-hour): "  << 1/mang_fitness[0] << endl;
        cout << "Fitness of ga solution: " << mang_fitness[0] << endl;

        // write to output file
        ofstream myfile(ten_file, ios::out);
        if (!myfile) {
            cout << "Error : Can't write output file best solution ";
            exit(EXIT_FAILURE);
        }

        myfile << "Fitness of ga solution: " << mang_fitness[0] << " - ";
        for(int i = 0; i < numTask; i++) {
            myfile << mang_NST[i] << " ";
        }
        myfile << endl;
        myfile.close();

        cout << "writed best solution to file done" << endl;
    }

    void xuat_file(char* ten_file) {
        ofstream myfile(ten_file, ios::out);
        if (!myfile ) {
            cout<< "Error : Can't write file "<< ten_file;
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
    void bruteForce();
    void readFile(char *);
};



// cuda kernel
//
//
//

void Population::host_init_NST() {
    int lengNST       = dataT->getNumTask();
    int soNST         = size_pop;
    int soMachine     = dataM->getNumMachine();

    mang_fitness      = new float[soNST];
    mang_makespan     = new int[soNST];
    rankings          = new int[2*sizePop];
    size_t sizeMemPop = sizeof(int)*lengNST*soNST;
    // Allocate the host input array h_A
    int *h_A = NULL;
    this->mang_NST = h_A = new int[lengNST*soNST];
    if (h_A == NULL) {
        fprintf(stderr, "Failed to allocate host vectors!\n");
        exit(EXIT_FAILURE);
    }
    // Allocate the device input array d_A
    int *d_A = NULL;
    CUDA_CALL( cudaMalloc((void **)&d_A, sizeMemPop) );

    // Call ramdom generator on GPU
    curandState *devStates;
    CUDA_CALL( cudaMalloc( (void **)&devStates, threadsPerBlock*sizeof(curandState) ) );

    // Launch kernel init
    int blocksPerGrid = (lengNST*soNST + threadsPerBlock - 1) / threadsPerBlock;
    printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
    
    device_init_NST<<<blocksPerGrid, threadsPerBlock>>>(d_A, soMachine , devStates);
    // CUDA_CALL( cudaGetLastError() );

    // Copy the device result vector in device memory to the host result vector in host memory.
    printf("Copy output data from the CUDA device to the host memory\n");
    CUDA_CALL( cudaMemcpy(h_A, d_A , sizeMemPop , cudaMemcpyDeviceToHost) );

    // Free device global memory
    CUDA_CALL( cudaFree(d_A) );
    CUDA_CALL( cudaFree(devStates) );
    // Reset the device and exit
    CUDA_CALL( cudaDeviceReset() );
}

void Population::GA_Evolution(int galoop) {
    if (galoop < 1 )
    {
        fprintf(stderr, "so the he khong hop le\n");
        exit(EXIT_FAILURE);
    }

    size_t A = sizeof(int)   * dataT->getNumTask() * size_pop;    // kich thuoc mang NST
    size_t B = sizeof(float) * size_pop;                    // kich thuoc mang fitness
    size_t R = sizeof(int)   * sizePop;

    int     *d_A, *d_temA, *d_temNST;       // mang NST
    float   *d_B, *d_temB, *d_temFitness; // mang fitness
    int     * d_charts;         // device array rankings

    CUDA_CALL(cudaMalloc((void**)&d_A,A));
    CUDA_CALL(cudaMalloc((void**)&d_temA,A));
    CUDA_CALL(cudaMalloc((void**)&d_temNST,A));
    CUDA_CALL(cudaMalloc((void**)&d_B,B));
    CUDA_CALL(cudaMalloc((void**)&d_temB,B));
    CUDA_CALL(cudaMalloc((void**)&d_temFitness,B));
    CUDA_CALL(cudaMalloc((void**)&d_charts,R)); 

    CUDA_CALL(cudaMemcpy(d_A, mang_NST, A, cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_temA, d_A, A, cudaMemcpyDeviceToDevice));
    

    // Phat viet, chuan bi du lieu cho ham tinh fitness
    int *tCore, *tMips, *tMemory, *tStorage, *tNetwork, *tStartTime, *tDuration;
    int *mCore, *mMips, *mPower, *mBasePower, *mMemory, *mStorage, *mNetwork;

    tCore = dataT->importFromTask(TaskCore);
    if (tCore == NULL)      {printf("Error at %s:%d\n",__FILE__,__LINE__); exit(EXIT_FAILURE);}
    
    tMips = dataT->importFromTask(TaskMips);
    if (tMips == NULL)      {printf("Error at %s:%d\n",__FILE__,__LINE__); exit(EXIT_FAILURE);}
    
    tMemory = dataT->importFromTask(TaskMemory);
    if (tMemory == NULL)    {printf("Error at %s:%d\n",__FILE__,__LINE__); exit(EXIT_FAILURE);}
    
    tStorage = dataT->importFromTask(TaskStorage);
    if (tStorage == NULL)   {printf("Error at %s:%d\n",__FILE__,__LINE__); exit(EXIT_FAILURE);}
    
    tNetwork = dataT->importFromTask(TaskNetwork);
    if (tNetwork == NULL)   {printf("Error at %s:%d\n",__FILE__,__LINE__); exit(EXIT_FAILURE);}
    
    tStartTime = dataT->importFromTask(TaskStartTime);
    if (tStartTime == NULL) {printf("Error at %s:%d\n",__FILE__,__LINE__); exit(EXIT_FAILURE);}
    
    tDuration = dataT->importFromTask(TaskDuration);
    if (tDuration == NULL)  {printf("Error at %s:%d\n",__FILE__,__LINE__); exit(EXIT_FAILURE);}

    mCore = dataM->importFromMachine(MachineCore);
    if (mCore == NULL)      {printf("Error at %s:%d\n",__FILE__,__LINE__); exit(EXIT_FAILURE);}
    
    mMips = dataM->importFromMachine(MachineMipsPerCore);
    if (mMips == NULL)      {printf("Error at %s:%d\n",__FILE__,__LINE__); exit(EXIT_FAILURE);}
    
    mPower = dataM->importFromMachine(MachinePower);
    if (mPower == NULL)     {printf("Error at %s:%d\n",__FILE__,__LINE__); exit(EXIT_FAILURE);}
    
    mBasePower = dataM->importFromMachine(MachineBasePower);
    if (mBasePower == NULL) {printf("Error at %s:%d\n",__FILE__,__LINE__); exit(EXIT_FAILURE);}
    
    mMemory = dataM->importFromMachine(MachineMemory);
    if (mMemory == NULL)    {printf("Error at %s:%d\n",__FILE__,__LINE__); exit(EXIT_FAILURE);}
    
    mStorage = dataM->importFromMachine(MachineStorage);
    if (mStorage == NULL)   {printf("Error at %s:%d\n",__FILE__,__LINE__); exit(EXIT_FAILURE);}
    
    mNetwork = dataM->importFromMachine(MachineNetwork);
    if (mNetwork == NULL)   {printf("Error at %s:%d\n",__FILE__,__LINE__); exit(EXIT_FAILURE);}
    // end

    int *eventStartList, *eventEndList;
    CUDA_CALL(cudaMalloc((void**)&eventStartList,sizeof(int)*2*numTask)); 
    CUDA_CALL(cudaMalloc((void**)&eventEndList,sizeof(int)*2*numTask)); 
    // Code added: Create 2 list eventStartList and eventEndList
    cout << "Begin create list\n";
    int startList[numTask * 2], endList[numTask * 2];
    for(int i = 0; i < numTask; i++){
        startList[i] = dataT->getStartTime(i);
        endList[i + numTask] = startList[i + numTask] = i;
        endList[i] = startList[i] + dataT->getDuration(i);
    }

    thrust::stable_sort_by_key(startList, startList + numTask, startList + numTask);
    thrust::stable_sort_by_key(endList, endList + numTask, endList + numTask);

    CUDA_CALL(cudaMemcpy(eventStartList, startList, sizeof(int) * 2 * numTask, cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(eventEndList, endList, sizeof(int) * 2 * numTask, cudaMemcpyHostToDevice));

    /*
    Test code for 2 lists
    cout << "Test 2 lists\n";
    for(int i = 0; i < numTask; i++){
        cout << startList[i] << " " << startList[i+numTask] << " " << endList[i] << " " << endList[i+numTask] << endl;
    }
    End code added
    */

    int blocksPerGrid =(size_pop + threadsPerBlock - 1) / threadsPerBlock;
    int blocksPerGridLarge = (size_pop*numTask + threadsPerBlock - 1) / threadsPerBlock;
    eval<<< blocksPerGrid, threadsPerBlock>>>(d_A, d_B, mCore, mPower, mMips, mBasePower, mMemory, mStorage, mNetwork, tCore, tMips, tMemory, tStorage, tNetwork, tStartTime, tDuration, eventStartList, eventEndList);         //khong tinh duoc
    // CUDA_CALL( cudaGetLastError() );

// Call ramdom generator on GPU
    curandState *devStates;
    CUDA_CALL( cudaMalloc( (void **)&devStates, threadsPerBlock*sizeof(curandState) ) );
    setup_curand<<<1,threadsPerBlock>>>(devStates);
    // CUDA_CALL( cudaGetLastError() );
 
    clock_t t;
    t = clock();
    for (int g = 0; g < galoop; g++) {
        lai_ghep<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_temA, devStates);
        // CUDA_CALL( cudaGetLastError() );
        eval<<< blocksPerGrid, threadsPerBlock>>>(d_temA, d_temB, mCore, mPower, mMips, mBasePower, mMemory, mStorage, mNetwork, tCore, tMips, tMemory, tStorage, tNetwork, tStartTime, tDuration, eventStartList, eventEndList);           //khong tinh duoc
        // CUDA_CALL( cudaGetLastError() );

        sort_Pop(d_B, d_temB);
        CUDA_CALL(cudaMemcpy(d_charts , rankings, R, cudaMemcpyHostToDevice));
        chon_loc_with_rankings<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_temA,d_temNST, d_B, d_temB, d_temFitness, d_charts);
        // CUDA_CALL( cudaGetLastError() );
        CUDA_CALL(cudaMemcpy(d_A, d_temNST, A, cudaMemcpyDeviceToDevice));
        CUDA_CALL(cudaMemcpy(d_temA, d_temNST, A, cudaMemcpyDeviceToDevice));
        CUDA_CALL(cudaMemcpy(d_B, d_temFitness, B, cudaMemcpyDeviceToDevice));
        CUDA_CALL(cudaMemcpy(d_temB, d_temFitness, B, cudaMemcpyDeviceToDevice));
        dot_bien<<<blocksPerGridLarge, threadsPerBlock>>>(d_temA, devStates);
        **CUDA_CALL( cudaGetLastError() );
        eval<<< blocksPerGrid, threadsPerBlock>>>(d_temA, d_temB, mCore, mPower, mMips, mBasePower, mMemory, mStorage, mNetwork, tCore, tMips, tMemory, tStorage, tNetwork, tStartTime, tDuration, eventStartList, eventEndList);           //khong tinh duoc
        // CUDA_CALL( cudaGetLastError() );
        sort_Pop(d_B, d_temB);
        CUDA_CALL(cudaMemcpy(d_charts , rankings, R, cudaMemcpyHostToDevice));
        chon_loc_with_rankings<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_temA,d_temNST, d_B, d_temB, d_temFitness, d_charts);
        // CUDA_CALL( cudaGetLastError() );
        CUDA_CALL(cudaMemcpy(d_A, d_temNST, A, cudaMemcpyDeviceToDevice));
        CUDA_CALL(cudaMemcpy(d_temA, d_temNST, A, cudaMemcpyDeviceToDevice));
        CUDA_CALL(cudaMemcpy(d_B, d_temFitness, B, cudaMemcpyDeviceToDevice));
        CUDA_CALL(cudaMemcpy(d_temB, d_temFitness, B, cudaMemcpyDeviceToDevice));
        t = clock() - t;
        if ( ((float)t)/CLOCKS_PER_SEC > timeStopGa ) {
            break;
        }
    }

// Copy the device result vector in device memory to the host result vector in host memory.
    printf("Copy output solution from the CUDA device to the host memory\n");
    CUDA_CALL( cudaMemcpy(mang_NST, d_A, A, cudaMemcpyDeviceToHost));
    CUDA_CALL( cudaMemcpy(mang_fitness, d_B, B,  cudaMemcpyDeviceToHost));
// Free device global memory
    CUDA_CALL(cudaFree(d_A));
    CUDA_CALL(cudaFree(d_B));
    CUDA_CALL(cudaFree(d_temA));
    CUDA_CALL(cudaFree(d_temB));
    CUDA_CALL( cudaFree(devStates) );
// Reset the device and exit
    CUDA_CALL( cudaDeviceReset() );
}

void Population::sort_Pop ( float* dev_fitnessOfParents, float* dev_fitnessOfOffspring) {
    // reset rankings 
    thrust::sequence(rankings, rankings + 2*sizePop);
    // copy key sort 
    float keySort[2*sizePop];
    // wrap device pointer
    thrust::device_ptr<float> dev_array1(dev_fitnessOfParents); 
    thrust::device_ptr<float> dev_array2(dev_fitnessOfOffspring);
    thrust::copy(  dev_array1,  dev_array1 + sizePop , keySort );
    thrust::copy(  dev_array2,  dev_array2 + sizePop, keySort + sizePop );
    // sort rankings[] by key
    thrust::stable_sort_by_key(keySort, keySort + 2*sizePop , rankings , thrust::greater<float>());
    // done
}

void Machines::host_init_Machine() {
    mipsPerCore = new int[nMachine];
    cores       = new int[nMachine];
    maxPower    = new int[nMachine];
    basePower   = new int[nMachine];
    memory      = new int[nMachine];
    storage     = new int[nMachine];
    network     = new int[nMachine];
    for(int i = 0; i < nMachine; i++){
        mipsPerCore[i] = 1;
        cores[i]       = 1;
        maxPower[i]    = 1;
        basePower[i]   = 1;
        memory[i]      = 1;
        storage[i]     = 1;
        network[i]     = 1;
    }
}

void Tasks::createTask(int id, int cores, int mips, int mem, int strg, int net, int startTime, int dur){
    cost[id]        = 0;
    duration[id]    = dur;
    start_time[id]  = startTime;
    ncore[id]       = cores;
    nmips[id]       = mips;
    memory[id]      = mem;  //my_host_rand(1, 4096); // Unit in MBytes
    storage[id]     = strg; //my_host_rand(50, 5000);
    network[id]     = net;  ///my_host_rand(1, 10000);        //1 Mb/s - 10000 Mb/s
}

void Tasks::host_init_Task() {
    cost        = new int[nTask];
    start_time  = new int[nTask];
    ncore       = new int[nTask];
    nmips       = new int[nTask];
    duration    = new int[nTask];
    memory      = new int[nTask];    
    storage     = new int[nTask];
    network     = new int[nTask];
    for(int i = 0; i < nTask; i++){
        cost[i]         = 0;
        start_time[i]   = 0;
        ncore[i]        = 0;
        nmips[i]        = 0;
        duration [i]    = 0;
        memory[i]       = 0;
        storage[i]      = 0;
        network[i]      = 0;
    }
}





/*
 * Function name: bestFit
 * Function description: find a solution using best fit algorithm and write out into bestfit.txt file
 * Input: N/A
*/

void Population::bestFit(){
    //Code added: Create 2 list eventStartList and eventEndList
    int startList[numTask * 2], endList[numTask * 2];
    for(int i = 0; i < numTask; i++) {
        startList[i] = dataT->getStartTime(i);
        endList[i + numTask] = startList[i + numTask] = i;
        endList[i] = startList[i] + dataT->getDuration(i);
    }

    thrust::stable_sort_by_key(startList, startList + numTask, startList + numTask);
    thrust::stable_sort_by_key(endList, endList + numTask, endList + numTask);

    int mips[numMachine], cores[numMachine], bestMachine;       //mpp = mips per power
    float bestmpp, mpp;
    int NST[numTask];

    // Initialize variables
    for(int i=0; i <numMachine; i++){
        mips[i] = cores[i] = 0;
    }
    for(int i=0; i <numTask; i++){
        NST[i] = -1;
    }

    int startPtr = 0, endPtr = 0, task;

    while(endPtr != numTask){
        bestmpp = -1;
        bestMachine = -1;
        
        // case 1: started event begin sooner than ended event
        if(startPtr != numTask && startList[startPtr] <= endList[endPtr]){
            task = startList[startPtr + numTask];
            startPtr++;
            for(int j=0; j<numMachine; j++){
            // Checking for invalid assignment
                if (dataM->getMipsPerCore(j) < dataT->getMips(task)
                    || dataM->getMipsPerCore(j) * dataM->getCores(j) < mips[j] + dataT->getMips(task) * dataT->getCores(task)
                    || dataM->getCores(j) < cores[j] + dataT->getCores(task)
                    ){
                        continue;
                }
                // compare the machine with previous machine to choose out the best machine which has minimum MIPS/power
                mpp = (float)dataM->getMipsPerCore(j) * dataM->getCores(j) / dataM->getMaxPower(j);
                if( mpp > bestmpp){
                    bestmpp = mpp;
                    bestMachine = j;
                }
            }
            // assign the best machine to task
            NST[task] = bestMachine;
            if(bestMachine < 0) break;
            mips[bestMachine] += dataT->getMips(task) * dataT->getCores(task);
            cores[bestMachine] += dataT->getCores(task);
        }
        // case 2: started event begin later than ended event
        else if(startPtr == numTask || endList[endPtr] < startList[startPtr]){
            task = endList[endPtr + numTask];
            endPtr++;
            bestMachine = NST[task];
            mips[bestMachine] -= dataT->getMips(task) * dataT->getCores(task);
            cores[bestMachine] -= dataT->getCores(task);
        }
    }

    bool skip = 0;

    for(int i = 0; i < numTask; i++){
        if(NST[i] < 0){
            skip = 1;
            break;
        }
    }

    float fitness = 0;
    // Calculate fitness
    if (skip)
        fitness = -1;
    else{
        float power[numMachine], powerDatacenter = 0;
        for(int i = 0; i< numMachine; i++){
            power[i] = 0;
        }

        for(int i = 0; i < numTask; i++){
            power[NST[i]] += (float)dataT->getMips(i) * dataT->getCores(i) 
                          / (dataM->getMipsPerCore(NST[i]) * dataM->getCores(NST[i])) * (dataM->getMaxPower(NST[i]) 
                          - dataM->getBasePower(NST[i])) * dataT->getDuration(i);
        }
        for(int j = 0; j < numMachine; j++){
            if(power[j] > 0)
                powerDatacenter += dataM->getBasePower(j);
            powerDatacenter += power[j];
        }
        fitness = (float) 1/powerDatacenter;
    }
    
    // Print out to file bestfit.txt
    ofstream myfile("bestfit.txt", ios::out);
        if (!myfile) {
            cout << "Error : Can't write file bestfit.txt";
            exit(EXIT_FAILURE);
        }
        cout << "Best fit solution - " << (float) fitness << endl;
        myfile << fitness << " - ";
        for( int j = 0; j < numTask; j++) {
            myfile << NST[j] << " ";
        }
        myfile << endl;         
    myfile.close();
}

void Machines::createMachine(int id, int core, int mips, int mem, int strg, int net, int maxPwr, int basePwr){
    maxPower[id]    = maxPwr;
    cores[id]       = core;
    mipsPerCore[id] = mips;
    memory[id]      = mem;
    storage[id]     = strg;
    network[id]     = net;
    basePower[id]   = basePwr;
}

void Population::printNST(char* ten_file, int *nst, float fitness) {
        ofstream myfile(ten_file, ios::out);
        if (!myfile) {
            cout<< "Error : Can't write file nst.txt";
            exit(EXIT_FAILURE);
        }
        myfile << fitness <<" - ";
        for(int i = 0; i < numTask; i++) {
            myfile << nst[i] << " ";
        }
        myfile.close();
}

void Population::readFile (char * fname){
    cout << "Begin read file\n";
    string line;

    ifstream infile(fname, ifstream::in);
    if(!infile.is_open()){
        cout << "File open error\n";
        return;
    }
    cout << "File is opened\n";

    bool type = 0;
    getline(infile, line);
    
    int ntasks = 0, nmachines = 0;
    // Edit code from here xet tao task hay machine dua vao bien type
    while(!infile.eof()){
        if(line[0] != ';' && line[0] != 0){
            if(line[0] == 'T'){
                type = 0;
            }
            else if(line[0] == 'M'){
                type = 1;
            }
            else if(type == 0 && ntasks < numTask){
                parseValue(line, type);
                ntasks++;
            }
            else if(type == 1 && nmachines < numMachine){
                parseValue(line, type);
                nmachines++;
            }
        }
        if(ntasks == numTask && nmachines == numMachine)
            break;
        getline(infile,line);
    }
    // End
    
    infile.close();
    dataT->showTask();
    dataM->showMachine();
}

void Population::parseValue(string line, bool type){
    // cout << "Begin parse value\n";
    float field[MAX_FIELD];

    // filter out space and extract field
    int index = 0;
    stringstream ss(line);
    while(ss.peek() != EOF){
        ss >> field[index];
        index++;
    }
    extractField(field, type);
}

void Population::extractField( float* fieldArray, bool type){
    // cout << "Begin extract field\n";
    // Create variable ready for task
    int id          = fieldArray[ID];
    int cores       = fieldArray[CORES];
    int mipsPerCore = fieldArray[MIPS];
    int memory      = fieldArray[MEMORY];
    int storage     = fieldArray[STORAGE];
    int network     = fieldArray[NETWORK];

    if(type == 0){
        int startTime = fieldArray[TASK_START_TIME], duration = fieldArray[DURATION];
        
    /*    
    curandGenerator_t generator;
    curandCreateGenerator (&generator, CURAND RNG PSEUDO MTGP32);
    curandSetPseudoRandomGeneratorSeed(generator, 1234ULL);
    unsigned int startTime, duration, stTemp[numTask], durTemp[numTask];
    curandGenerate (generator, stTemp, numTask);
    curandGenerate (generator, durTemp, numTask);
    startTime = stTemp[id]%5000;
    duration = durTemp[id]%5000 + 1;
    */

        dataT->createTask(id, cores, mipsPerCore, memory, storage, network, startTime, duration);
    }
    else{
        int maxPwr  = fieldArray[MAX_POWER];
        int basePwr = fieldArray[BASE_POWER];
        dataM->createMachine(id, cores, mipsPerCore, memory, storage, network, maxPwr, basePwr);
    }
}

int main() {
    float rateMutan = rateMut;
    rateMutan = rateMutan > 1 ? rateMutan/1000 : rateMutan;
    Population P;
    P.readFile("data.txt");
    // P.bruteForce();
    // P.xuat_file("nst.txt");
    
    struct timeval timeRunOfGa;  
    gettimeofday(&timeRunOfGa, NULL);  
    double dTime1 = timeRunOfGa.tv_sec+(timeRunOfGa.tv_usec/1000000.0);  
    P.GA_Evolution(generation);
    
    gettimeofday(&timeRunOfGa, NULL);  
    double dTime2 = timeRunOfGa.tv_sec+(timeRunOfGa.tv_usec/1000000.0); 
    
    // P.xuat_file("after_ga.txt");
    P.printBest("best_solution.txt");
    P.bestFit();
    cout << endl << "total time run ga: "<< dTime2 - dTime1 << endl<<endl;
    cout << "Test parameters: machines " << numMachine << " tasks " << numTask 
         << " population size " << sizePop << " generations " << generation << " threads " 
         << threadsPerBlock << " Mutation rate: " << rateMutan << endl;
    cout <<"Done\n";
    
    return 0;
}

