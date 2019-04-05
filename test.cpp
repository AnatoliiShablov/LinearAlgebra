#include "linearalgebra.h"

int main() {

    la::tensor a;

    std::cin >> a;

    std::printf("result tensor:\n");
    a.sym({true, false, true});
    std::cout << a;
}


