#include <iostream>
#include <vector>
#include <numeric>

void initLoop(std::vector<uint16_t>& nbhdRow, size_t i) {
    uint16_t k = 0;
    for (size_t j = 0; j < nbhdRow.size(); ++j, ++k) {
        nbhdRow[j] = (i == j) ? ++k : k;
    }
}

/**
 * Initializes elements of nbhdRow without loops.
 *
 * @param nbhdRow The container to initialize.
 * @param i       The reference index.
 */
void initIota(std::vector<uint16_t>& nbhdRow, size_t i) {
    iota(nbhdRow.begin(), nbhdRow.begin() + i, 0);
    iota(nbhdRow.begin() + i, nbhdRow.end(), i+1);
}

int main() {
    size_t n = 5; // The size of the row vector
    std::vector<uint16_t> nbhdRow(n);
    size_t i = 2;

    initLoop(nbhdRow, i);

    std::cout << "With Loop: " << std::endl;
    for (const auto& element : nbhdRow) {
        std::cout << element << " ";
    }

    initIota(nbhdRow, i);

    std::cout << "With Iota: " << std::endl;
    for (const auto& element : nbhdRow) {
        std::cout << element << " ";
    }

    return 0;
}
