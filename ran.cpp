#pragma offload_attribute(push, target(mic))
#include <stdlib.h>
#include <iostream>
#include <omp.h>
#pragma offload_attribute(pop)

using namespace std;

int main() {
    srand48(time(NULL));    

    double t1 = omp_get_wtime();
//#pragma offload target (mic: 0)
//{
#pragma omp parallel for
    for (int i = 0; i < 50000000; i++) {
        int a = lrand48();
        //cout << a << endl;
    }    
//}
    t1 = omp_get_wtime() - t1;
    cout << "T1 " << t1 << endl;
    return 0;
}
