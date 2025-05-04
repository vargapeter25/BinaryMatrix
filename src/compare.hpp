#pragma once

#include <cstdint>
#include <functional>
#include "binary_matrix.hpp"

template<uint32_t N, uint32_t K>
bool _check_prec_exact1(SmallBinaryMatrix<N> m1, SmallBinaryMatrix<N> m2){
    
    if(m1.get_data() == m2.get_data()) return false;

    if constexpr (K <= 1) {
        return (m2.get_data() & m1.get_data()) == m1.get_data();
    } else {
        if(!_check_prec_exact1<N, K-1>(m1, m2)) return false;
        for(auto rows : subset_mask<N, K>) {
            for(auto cols : subset_mask<N, K>) {
                auto m12 = m1.template sub_matrix<K>(rows, cols);
                auto m22 = m2.template sub_matrix<K>(rows, cols);
                if(!m12.singular() && m22.singular()) return false;
            }
        }
        return true;
    }
}

template<uint32_t N>
bool check_prec_exact1(SmallBinaryMatrix<N> m1, SmallBinaryMatrix<N> m2){
    return _check_prec_exact1<N, N-1>(m1, m2); 
}

template<uint32_t N, uint32_t K>
bool _check_prec_exact2(SmallBinaryMatrix<N> m1, SmallBinaryMatrix<N> m2){
    if constexpr (K > 1) {
        if(!_check_prec_exact2<N, K-1>(m1, m2)) return false;
    } 
    for(auto rows : subset_mask<N, K>) {
        for(auto cols : subset_mask<N, K>) {
            auto m12 = m1.template sub_matrix<K>(rows, cols);
            auto m22 = m2.template sub_matrix<K>(rows, cols);
     
            const auto mask_inv = (1<<N) - 1;
            auto m_comp = m1.template sub_matrix<N - K>(mask_inv ^ rows, mask_inv ^ cols);
            if(!m_comp.singular() && (!m12.singular() && m22.singular())) return false;
        }
    }
    return true;
}

template<uint32_t N>
bool check_prec_exact2(SmallBinaryMatrix<N> m1, SmallBinaryMatrix<N> m2){
    return _check_prec_exact2<N, N-1>(m1, m2);
}

template<uint32_t N>
bool check_prec(SmallBinaryMatrix<N> m1, SmallBinaryMatrix<N> m2, std::function<bool(SmallBinaryMatrix<N>, SmallBinaryMatrix<N>)> comp, SmallBinaryMatrix<N>* out_prec = nullptr) {
    uint32_t row_ord[N], col_ord[N];
    for(uint32_t i = 0; i < N; i++) row_ord[i] = i;
    do {
        for(uint32_t i = 0; i < N; i++) col_ord[i] = i;
        do {
            SmallBinaryMatrix<N> tmp1 = m1;
            SmallBinaryMatrix<N> tmp2 = m1;
            
            tmp1.reorder_rows(row_ord);
            tmp1.reorder_cols(col_ord);
            
            tmp2.transpose();
            tmp2.reorder_rows(row_ord);
            tmp2.reorder_cols(col_ord);
            
            if (comp(tmp1, m2)) {
                if(out_prec) *out_prec = tmp1;
                return true;
            }

            if (comp(tmp2, m2)) { // check transposed
                if(out_prec) *out_prec = tmp2;
                return true;
            }
        } while(std::next_permutation(col_ord, col_ord + N));
    } while(std::next_permutation(row_ord, row_ord + N));
    return false;
}

template<uint32_t N>
bool check_prec1(SmallBinaryMatrix<N> m1, SmallBinaryMatrix<N> m2, SmallBinaryMatrix<N>* out_prec = nullptr) {
    return check_prec<N>(m1, m2, check_prec_exact1<N>, out_prec);
}

template<uint32_t N>
bool check_prec2(SmallBinaryMatrix<N> m1, SmallBinaryMatrix<N> m2, SmallBinaryMatrix<N>* out_prec = nullptr) {
    return check_prec<N>(m1, m2, check_prec_exact2<N>, out_prec);
}