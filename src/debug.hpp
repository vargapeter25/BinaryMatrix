#pragma once
#include "binary_matrix.hpp"
#include "utils.hpp"
#include <iostream>

template<uint32_t N>
void print(SmallBinaryMatrix<N> m){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++) {
            std::cout << m[i][j];
        }
        std::cout << '\n';
    }
}