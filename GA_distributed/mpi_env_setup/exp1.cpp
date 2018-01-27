#include <iostream>
#include <stdio.h>
#include "mpi.h"
#include <unistd.h>
#include <string>
#include <fstream>

using namespace std;

int main(int argc, char* argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    //printf("I am %d of %d\n", rank, size);
    printf("Hello world! I have %ld logical processors.\n",
            sysconf(_SC_NPROCESSORS_ONLN ));

    string line, filename = "exp1.cpp";
    ifstream myfile ("/opt/share/tu/GA_MIC/GA_distributed/mpi_env_setup/exp1.cpp");
    if (myfile.is_open()) {
    
        while ( getline (myfile,line)) {
            cout << line << '\n';
        }
        myfile.close();
    } else cout << "Unable to open file\n";

    MPI_Finalize();
    return 0;
}
