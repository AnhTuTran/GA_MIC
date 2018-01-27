#include <iostream>
#include <stdio.h>
#include "mpi.h"
#include <unistd.h>
#include <string>
#include <fstream>

using namespace std;

int main(int argc, char* argv[]) {
    int rank, size;
    int best_sol;
    int *solutions;  

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    //printf("I am %d of %d\n", rank, size);
    printf("Hello world! I have %ld logical processors, rank %d.\n",
            sysconf(_SC_NPROCESSORS_ONLN ), rank);

    string line, filename = "exp1.cpp";
    ifstream myfile ("/opt/share/tu/GA_MIC/GA_distributed/mpi_env_setup/exp1.cpp");
    if (myfile.is_open()) {
    
        while ( getline (myfile,line)) {
            cout << line << '\n';
        }
        myfile.close();
    } else cout << "Unable to open file\n";

    MPI_Finalize();
//    return 0;

    // execute GA
    // gather the best solutions from slaves
    best_sol = rank + 10;
    solutions = new int[size];
    MPI_Gather(&best_sol, 1, MPI_INTEGER, solutions, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

    // choose the best solutions among nodes
    if (rank == 0) {
    for (int i = 0; i < size; i++) {
        cout << i << " " << solutions[i] << endl;
    }
    }
}
