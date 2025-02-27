#pragma once
#include "binary_matrix.hpp"
#include "utils.hpp"
#include <bitset>
#include <algorithm>

template<uint32_t N>
SmallBinaryMatrix<N> get_normal_form(SmallBinaryMatrix<N> m){
    uint64_t best = m.get_data();
    uint32_t ord_col[N];
    for(int i = 0; i < N; i++) ord_col[i] = i;
    do {
        SmallBinaryMatrix<N> cur = m;
        cur.reorder_cols(ord_col);
        uint32_t ord_row[N];
        for(int i = 0; i < N; i++) ord_row[i] = i;
        std::sort(ord_row, ord_row + N, [&cur](uint32_t i, uint32_t j) { return (uint8_t)cur[i] < (uint8_t)cur[j]; });
        cur.reorder_rows(ord_row);
        best = std::max(best, cur.get_data()); 
    } while(std::next_permutation(ord_col, ord_col + N));
    return SmallBinaryMatrix<N>(best);
}

template<uint32_t N>
void sort_rows(SmallBinaryMatrix<N>& m){
    uint32_t ord[N];
    for(int i = 0; i < N; i++) ord[i] = i;
    std::sort(ord, ord + N, [&m](int i, int j) { return __builtin_popcount((uint8_t)m[i]) < __builtin_popcount((uint8_t)m[j]); });
    m.reorder_rows(ord);
}

template<uint32_t N>
void sort_cols(SmallBinaryMatrix<N>& m){
    m.transpose();
    sort_rows(m);
    m.transpose();
}

template<uint32_t N>
uint64_t get_data(SmallBinaryMatrix<N> m) {
    uint64_t res = 0;
    for(int i = 0; i < N; i++){
        res = res | ((m[i] & ((1ULL<<N)-1))<<(N*i));
    }
    return res;
}

template<uint32_t N>
SmallBinaryMatrix<N> mat_from_data(uint64_t d){
    SmallBinaryMatrix<N> res;
    for(int i = 0; i < N; i++){
        res[i] = (d>>(i*N)) & ((1<<N) - 1);
    }
    return res;
}

template<uint32_t N, uint32_t K>
std::bitset<sub_matrix_cnt<N, N>()> get_sub_matrix_mask(SmallBinaryMatrix<N> m){
    std::bitset<sub_matrix_cnt<N, N>()> res;
    uint32_t idx = 0;
    if constexpr (K > 1){
        res = get_sub_matrix_mask<N, K-1>(m);
        idx = sub_matrix_cnt<N, K - 1>();
    }
    // cout << K << " start: " << idx << '\n';
    for(uint32_t row_mask : subset_mask<N, K>) {
        for(uint32_t col_mask : subset_mask<N, K>) {
            auto sub_mat = m.template sub_matrix<K>(row_mask, col_mask);
            res[idx++] = !sub_mat.singular();
        }
    }
    return res;
}
