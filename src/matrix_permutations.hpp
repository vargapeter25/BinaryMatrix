#pragma once

#include <cstdint>
#include <vector>
#include <array>
#include <algorithm>
#include "binary_matrix.hpp"

template<uint32_t N>
struct matrix_permutation{
    matrix_permutation(const std::array<uint32_t, N>& row, const std::array<uint32_t, N>& col) {
        std::copy(row.begin(), row.end(), row_ord.begin());
        std::copy(col.begin(), col.end(), col_ord.begin());
    }
    matrix_permutation(uint32_t* row, uint32_t* col) {
        std::copy(row, row + N, row_ord.begin());
        std::copy(col, col + N, col_ord.begin());
    }
    std::array<uint32_t, N> row_ord;
    std::array<uint32_t, N> col_ord;

    bool operator<(const matrix_permutation<N>& mp) const {
        return row_ord != mp.row_ord ? row_ord < mp.row_ord : col_ord < mp.col_ord;
    }

    bool operator==(const matrix_permutation<N>& mp) const {
        return row_ord == mp.row_ord && col_ord == mp.col_ord;
    }
};  

template<uint32_t N>
struct permutation_data{
    std::vector<matrix_permutation<N>> permutations;

    void normalize() { sort(permutations.begin(), permutations.end()); }

    bool operator<=(const permutation_data& pd) const { // data must be normalized
        if(permutations.size() > pd.permutations.size()) return false;

        for(auto &mp : pd.permutations){
            std::array<uint32_t, N> row_trans = prod(inv_permutation(permutations[0].row_ord), mp.row_ord);
            std::array<uint32_t, N> col_trans = prod(inv_permutation(permutations[0].col_ord), mp.col_ord);
            std::vector<matrix_permutation<N>> tmp_permutations;

            assert(prod(permutations[0].row_ord, row_trans) == mp.row_ord);
            assert(prod(permutations[0].col_ord, col_trans) == mp.col_ord);

            for(auto &mp2 : permutations){
                tmp_permutations.push_back(matrix_permutation<N>(prod(mp2.row_ord, row_trans), prod(mp2.col_ord, col_trans)));
            }
            sort(tmp_permutations.begin(), tmp_permutations.end()); // normalize tmp data
            if(check_subset(tmp_permutations, pd.permutations)) return true;
        }   

        for(auto &mp : pd.permutations){ // flipped coordinates
            std::array<uint32_t, N> row_trans = prod(inv_permutation(permutations[0].row_ord), mp.col_ord);
            std::array<uint32_t, N> col_trans = prod(inv_permutation(permutations[0].col_ord), mp.row_ord);
            std::vector<matrix_permutation<N>> tmp_permutations;

            assert(prod(permutations[0].row_ord, row_trans) == mp.col_ord);
            assert(prod(permutations[0].col_ord, col_trans) == mp.row_ord);

            for(auto &mp2 : permutations){
                tmp_permutations.push_back(matrix_permutation<N>(prod(mp2.col_ord, col_trans), prod(mp2.row_ord, row_trans)));
            }
            sort(tmp_permutations.begin(), tmp_permutations.end()); // normalize tmp data
            if(check_subset(tmp_permutations, pd.permutations)) return true;
        }   

        return false;
    }

    bool contains(const matrix_permutation<N>& mp) const { // data must be normalized
        auto it = lower_bound(permutations.begin(), permutations.end(), mp);
        return it != permutations.end() && (*it) == mp;
    }

private:
    std::array<uint32_t, N> inv_permutation(const std::array<uint32_t, N>& p) const {
        std::array<uint32_t, N> res;
        for(int i = 0; i < N; i++) res[p[i]] = i;
        return res;
    }

    std::array<uint32_t, N> prod(const std::array<uint32_t, N>& p1, const std::array<uint32_t, N>& p2) const {
        std::array<uint32_t, N> res;
        for(int i = 0; i < N; i++) res[i] = p2[p1[i]];
        return res;
    }

    bool check_subset(const std::vector<matrix_permutation<N>>& p1, const std::vector<matrix_permutation<N>> &p2) const { // !!! data must be normalized
        // sort(p1.begin(), p1.end());
        // sort(p2.begin(), p2.end());
        int j = 0;
        for(auto &mp : p2) {
            if(j < (int)p1.size() && p1[j] == mp) j++;
        }
        return j == (int)p1.size();
    }
};

template<uint32_t N>
permutation_data<N> get_permutation_data(SmallBinaryMatrix<N> m) { // Returns normalized data
    uint32_t row_ord[N], col_ord[N];
    permutation_data<N> res;
    for(int i = 0; i < N; i++) row_ord[i] = i;
    do {
        for(int i = 0; i < N; i++) col_ord[i] = i;
        do {
            SmallBinaryMatrix<N> tmp(m);
            tmp.reorder_rows(row_ord);
            tmp.reorder_cols(col_ord);
            if(tmp.main_minors_non_singular() && tmp.anti_minors_non_singular()) res.permutations.push_back(matrix_permutation<N>(row_ord, col_ord));
        } while(std::next_permutation(col_ord, col_ord + N));
    } while(std::next_permutation(row_ord, row_ord + N));
    sort(res.permutations.begin(), res.permutations.end());
    return res;
}

template<uint32_t N, uint32_t K>
SmallBinaryMatrix<N> _reorder(SmallBinaryMatrix<N> m) {
    const int MASK_INV = (1<<N) - 1;
    const int mask = (1<<K) - 1;
    if constexpr (K >= N - 1) {
        return m;
    } else{
        for(int i = K; i < N; i++){
            for(int j = K; j < N; j++){
                int row_mask = mask | (1<<i);
                int col_mask = mask | (1<<j);
                auto sub_m = m.sub_matrix<K + 1>(row_mask, col_mask);
                auto sub_m2 = m.sub_matrix<N - K - 1>(row_mask^MASK_INV, col_mask^MASK_INV);
                if(!sub_m.singular() && !sub_m2.singular()) {
                    m.swap_rows(K, i);
                    m.swap_cols(K, j);
                    return reorder<K + 1>(m);
                }
            }
        }
    }
    return m;
}

template<uint32_t N>
SmallBinaryMatrix<N> reorder(SmallBinaryMatrix<N> m) {
    return _reorder<N, 0>(m);
}