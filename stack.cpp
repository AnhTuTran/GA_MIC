#include <iostream>
#include <string>
#include <stdio.h>
#include <sys/resource.h>

using namespace std;

int main()
{
    
    const rlim_t kStackSize = 512L * 1024L * 1024L;   // min stack size = 64 Mb
    struct rlimit rl;
    int result;

    rl.rlim_cur = kStackSize;
    result = setrlimit(RLIMIT_STACK, &rl);    
    if (result != 0) {
        printf("setrlimit returned result = %d\n", result);
    }
    /*
    const rlim_t kStackSize = 128L * 1024L * 1024L;   // min stack size = 64 Mb
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {   
        cout << rl.rlim_cur << endl;
        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
        }
    }*/
    cout << rl.rlim_cur << endl;
    return 0;
}


