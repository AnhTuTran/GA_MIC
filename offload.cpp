#include <cstdio>
#include <stdlib.h>
#include <iostream>
#include <sys/time.h>

using namespace std;

#define numTask 10
#define sizePop 512

#pragma offload_attribute(push, target(mic)) 
int con = 111;
void a() {
    int num = 1, N =10;
    int *C = new int[N];
    #pragma offload_transfer target(mic: 0) in(C : length(N) free_if(0))
    printf("MIC %d %d\n", con, C[2]);
}
#pragma offload_attribute(pop)

int N = 10;
/*
int main() {
    printf("Hello World from host!\n"); 
    //int N = 10;
    int *A = (int*) malloc(sizeof(int)*N);
    int *B = (int*) malloc(sizeof(int)*N);
    
#pragma offload_transfer target(mic: 0) in(A : length(N) free_if(0))

#pragma offload target(mic: 0) in(A : length(N) alloc_if(0) free_if(0)) 
    {
        //printf("Hello World from coprocessor!\n"); fflush(stdout);
        for (int i = 0; i < N; i++)
            A[i] = 1;   
        a();
    }

    A[3] = 3;

#pragma offload target(mic: 0) nocopy(A : length(10))
    {
        printf("Hello World from coprocessor! %d\n", A[3]); fflush(stdout);
    }

    
    printf("CPU %d\n", con);
    printf("Hello World from cpu! %d\n", A[3]); fflush(stdout);

    printf("Bye\n"); 
}
*/
class Population {
public:
    int * mang_NST;      
    Population() {
        mang_NST = NULL;
    }
    ~Population() {  cout << "0 \n";
        delete [] mang_NST; cout << "1 \n";
    }
__attribute__ ((target (mic))) void print(int *mang_NST) {
        for (int i = 0; i < 10; i++)
            printf("%d\n", mang_NST[i]);
    }

    void func() {
        int *mang_NST = this->mang_NST;
    }
};

  
int* host_init_NST () {
     int *mang_NST = new int[numTask*sizePop];
    #pragma offload target(mic: 0) inout(mang_NST : length(numTask*sizePop))
    {
        for (int i = 0; i < numTask*sizePop; i++ ) {
            mang_NST[i] = i;
        }
    }

    cout<<"Creating NSTs done\n";
    return mang_NST;
}   
#pragma offload_attribute(push, target(mic))
int my_host_rand (int s, int d) {
    return rand() % (d-s+1) + s;
}
#pragma offload_attribute(pop) 

/*
int main() {
    int *A = new int[10];
    srand(time(NULL));
    #pragma offload target(mic: 0) inout(A : length(10)) 
    {
        for (int i = 0; i < 10; i++ ) { 
            A[i] =  my_host_rand(0, 100);
        }   
    }    

    for (int i = 0; i < 10; i++ ) { 
        cout << A[i] << endl ;
    }
  
    srand(time(NULL));
    Population p;     
    p.mang_NST = host_init_NST();
#pragma offload target(mic: 0) inout(p.mang_NST : length(numTask*sizePop))
    {    
        p.print(p.mang_NST);
    }
    p.func();
    return 0;
}
*/

int main() {
    int *A = new int[5];
    int B[5];
    B[:] = 3;
    for(int i = 0; i < 5; i++) A[i] = 5;
    #pragma offload target(mic) in(A : length(5) alloc_if(1) free_if(0)) in(B)
    {
    }

    #pragma offload target(mic) out(A : length(5) alloc_if(0) free_if(1))
    {
        for(int i = 0; i < 5; i++) { printf("MIC A %d\n", A[i]);  printf("MIC B %d\n", B[i]);}
        for(int i = 0; i < 5; i++) { A[i] = i; B[i] = 4 -i;}
    }
    for(int i = 0; i < 5; i++) { printf("Cpu A %d\n", A[i]);  printf("Cpu B %d\n", B[i]);}
    return 0;
}
