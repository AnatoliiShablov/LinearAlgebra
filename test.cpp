#include "linearalgebra.h"

int main() {

    la::tensor tensor(2, 2, 2);
    std::vector<size_t> i(4);
    std::printf("tensor:\n");
    for (i[2] = 0; i[2] < 2; i[2]++)
        for (i[3] = 0; i[3] < 2; i[3]++)
            for (i[0] = 0; i[0] < 2; i[0]++)
                for (i[1] = 0; i[1] < 2; i[1]++) {
                    std::scanf("%lf", &tensor(i));
                }
    la::matrix T(2, 2);
    std::printf("matrix:\n");
    for (size_t j = 0; j < 2; j++) {
        for (size_t k = 0; k < 2; k++) {
            std::scanf("%lf", &T(j, k));
        }
    }

    tensor = tensor.change_basis(T);
    std::printf("result tensor:\n");
    for (i[2] = 0; i[2] < 2; i[2]++) {
        for (i[3] = 0; i[3] < 2; i[3]++) {
            for (i[0] = 0; i[0] < 2; i[0]++) {
                for (i[1] = 0; i[1] < 2; i[1]++) {
                    std::printf("%lf ", tensor(i));
                }
                std::printf("\n");
            }
            std::printf("-----\n");
        }
        std::printf("-----\n");
    }
}