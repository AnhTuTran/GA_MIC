#include <cstdio> 

int main() {
    printf("Hello World from host!\n"); 

#pragma offload target(mic: 1)
    {
        printf("Hello World from coprocessor!\n"); fflush(stdout);
    }

    printf("Bye\n"); 
}
