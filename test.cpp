#include "linearalgebra.h"

int main() {

    la::tensor a(2, 2, 2);

    std::vector<size_t> i(4);
    for (i[2] = 0; i[2] < 2; i[2]++) {
        for (i[3] = 0; i[3] < 2; i[3]++) {
            for (i[0] = 0; i[0] < 2; i[0]++) {
                for (i[1] = 0; i[1] < 2; i[1]++) {
                    std::scanf("%lf", &a(i));
                }
            }
        }
    }

    la::matrix T(2, 2);

    for (size_t j = 0; j < 2; j++) {
        for (size_t k = 0; k < 2; k++) {
            std::scanf("%lf", &T(j, k));
        }
    }

    std::printf("result tensor:\n");
    a.change_basis(T);
    for (i[2] = 0; i[2] < 2; i[2]++) {
        for (i[3] = 0; i[3] < 2; i[3]++) {
            for (i[0] = 0; i[0] < 2; i[0]++) {
                for (i[1] = 0; i[1] < 2; i[1]++) {
                    std::printf("%lf ", a(i));
                }
                std::printf("\n");
            }
            std::printf("----------\n");
        }
        std::printf("----------\n");
    }
}


