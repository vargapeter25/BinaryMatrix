#pragma once
#include <stdint.h>

constexpr uint64_t factorial(uint64_t n) {
    return n <= 1 ? 1 : factorial(n - 1) * n;
}

constexpr uint64_t binomial(uint64_t n, uint64_t k){
    return n < k ? 0 : factorial(n) / factorial(k) / factorial(n - k);
}

template<uint32_t M, uint32_t K>
constexpr int sub_matrix_cnt() {
    if constexpr (K > 1){
        return sub_matrix_cnt<M, K - 1>() + binomial(M, K) * binomial(M, K);
    } else{
        return binomial(M, K) * binomial(M, K);
    }
}

template<uint32_t N, uint32_t K>
constexpr auto subset_mask = []() constexpr{
    constexpr size_t size = binomial(N, K);
    std::array<uint32_t, size> result{};
    int idx = 0;
    for(int i = 0; i < (1<<N); i++){
        if(__builtin_popcount(i) == K) {
            result[idx++] = i;
        }
    }
    return result;
}();