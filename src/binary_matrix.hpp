#pragma once
#include <stdint.h>
#include <array>
#include <bitset>
#include <cstring>
#include <bit>
#include <cassert>

#include <iostream>

template<uint32_t N>
requires requires() { N <= 8; }
class SmallBinaryMatrix{
public:
    struct CellWrapper{
        uint8_t *ptr;
        uint32_t i;
        operator bool() const { return ((*ptr)>>i)&1; }
        CellWrapper& operator=(const CellWrapper& cw) { *this = (bool)cw; return *this; }
        CellWrapper& operator=(bool val) {
            if(val != (((*ptr)>>i)&1)) (*ptr) ^= (1U<<i);
            return *this;
        }
    };
    struct RowWrapper{
        uint8_t *ptr;
        operator uint8_t() { return *ptr; }
        CellWrapper operator[](uint32_t i) { return CellWrapper{ptr, i}; }
        RowWrapper& operator+=(RowWrapper r) { (*ptr) ^= (*r.ptr); return *this; }
        RowWrapper& operator=(uint8_t d) { *ptr = d; return *this; }
        RowWrapper& operator=(const RowWrapper& rw) { return *this = (uint8_t)rw; }
    };

    SmallBinaryMatrix(uint64_t _data = 0) : data(_data & MASK) {}

    RowWrapper operator[](uint32_t i) { return get_row(i); }
    RowWrapper get_row(uint32_t i) { return RowWrapper{reinterpret_cast<uint8_t*>(&data) + i}; }

    void swap_rows(uint32_t i, uint32_t j) { if(i == j) return; auto r1 = get_row(i), r2 = get_row(j); r1 += r2; r2 += r1; r1 += r2; }
    void swap_cols(uint32_t i, uint32_t j) { transpose(); swap_rows(i, j); transpose(); }

    void reorder_rows(uint32_t *ord) {
        SmallBinaryMatrix new_mat;
        for(uint32_t i = 0; i < N; i++){
            new_mat[i] += get_row(ord[i]);
        }
        data = new_mat.data;
    }

    void reorder_cols(uint32_t *ord) {
        transpose();
        reorder_rows(ord);
        transpose();
    }

    void flip_horizontal() { data = __builtin_bswap64(data); data >>= (8 - N) * 8; }
    void flip_vertical() { transpose(); flip_horizontal(); transpose(); }

    void transpose() {
        uint64_t new_data = 0;
        for(uint32_t i = 0; i < N; i++){
            new_data |= row_transpose[*get_row(i).ptr]<<i;
        }
        data = new_data;
    }
    void anti_transpose() { flip_horizontal(); transpose(); flip_horizontal(); }

    inline uint32_t first_row_with_one(uint32_t r, uint32_t c) {
        return std::countr_zero((row_transpose[256 - (1<<r)]<<c) & data) >> 3;
    }

    uint64_t get_data() { return data; }

    bool singular() const {
        SmallBinaryMatrix tmp(data);
        for(uint32_t i = 0; i < N; i++){
            uint32_t r = tmp.first_row_with_one(i, i);
            if(r == 8) return true;
            if(r > i) {
                tmp.swap_rows(i, r);
            }
            auto cur_row = tmp[i];
            tmp.data ^= (((row_transpose[256 - (1<<(i+1))]<<i) & tmp.data) >> i) * (uint8_t)cur_row;
        }
        return false;
    }
    
    bool main_minors_non_singular() const {
        SmallBinaryMatrix tmp(data);
        for(uint32_t i = 0; i < N; i++){
            uint32_t r = tmp.first_row_with_one(i, i);
            if(r == 8) return false;
            if(r > i) return false;
            auto cur_row = tmp[i];
            tmp.data ^= (((row_transpose[256 - (1<<(i+1))]<<i) & tmp.data) >> i) * (uint8_t)cur_row;
        }
        return true;
    }

    bool anti_minors_non_singular() const {
        SmallBinaryMatrix tmp(data);
        tmp.flip_horizontal();
        tmp.flip_vertical();
        return tmp.main_minors_non_singular();
    }

    // requires requires() { M <= N; }
    template<uint32_t M>
    SmallBinaryMatrix<M> sub_matrix(uint32_t row_mask, uint32_t col_mask) {
        assert(__builtin_popcount(row_mask) == M); // must be square submatrix
        assert(__builtin_popcount(col_mask) == M); // must be square submatrix
        uint32_t ord_row[N], ord_col[N];
        int idx_row = 0, idx_col = 0;
        for(uint32_t i = 0; i < N; i++){
            if((row_mask>>i)&1) ord_row[idx_row++] = i;
            if((col_mask>>i)&1) ord_col[idx_col++] = i;
        }
        for(uint32_t i = 0; i < N; i++){
            if(!((row_mask>>i)&1)) ord_row[idx_row++] = i;
            if(!((col_mask>>i)&1)) ord_col[idx_col++] = i;
        }
        SmallBinaryMatrix tmp(data);
        tmp.reorder_rows(ord_row);
        tmp.reorder_cols(ord_col);
        return SmallBinaryMatrix<M>(tmp.data);
    }

    // template<uint32_t M>
    // requires requires() { M < N; }
    // SmallBinaryMatrix<M> sub_matrix(const std::array<uint32_t, M>& rows, const std::array<uint32_t, M>& cols) {
    //     uint32_t row_mask = 0, col_mask = 0;
    //     for(uint32_t r : rows) row_mask |= (1u<<r);
    //     for(uint32_t c : cols) col_mask |= (1u<<c);
    //     return sub_matrix<M>(row_mask, col_mask);
    // }

    bool operator==(SmallBinaryMatrix m) const { return data == m.data; }

    SmallBinaryMatrix operator*(SmallBinaryMatrix m) {
        SmallBinaryMatrix res;
        m.transpose();
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                res[i][j] = __builtin_popcount((uint8_t)get_row(i) & (uint8_t)m.get_row(j)) & 1;
            }
        }
        return res;
    }

    constexpr static const uint64_t row_transpose[256] = {0ULL, 1ULL, 256ULL, 257ULL, 65536ULL, 65537ULL, 65792ULL, 65793ULL, 16777216ULL, 16777217ULL, 16777472ULL, 16777473ULL, 16842752ULL, 16842753ULL, 16843008ULL, 16843009ULL, 4294967296ULL, 4294967297ULL, 4294967552ULL, 4294967553ULL, 4295032832ULL, 4295032833ULL, 4295033088ULL, 4295033089ULL, 4311744512ULL, 4311744513ULL, 4311744768ULL, 4311744769ULL, 4311810048ULL, 4311810049ULL, 4311810304ULL, 4311810305ULL, 1099511627776ULL, 1099511627777ULL, 1099511628032ULL, 1099511628033ULL, 1099511693312ULL, 1099511693313ULL, 1099511693568ULL, 1099511693569ULL, 1099528404992ULL, 1099528404993ULL, 1099528405248ULL, 1099528405249ULL, 1099528470528ULL, 1099528470529ULL, 1099528470784ULL, 1099528470785ULL, 1103806595072ULL, 1103806595073ULL, 1103806595328ULL, 1103806595329ULL, 1103806660608ULL, 1103806660609ULL, 1103806660864ULL, 1103806660865ULL, 1103823372288ULL, 1103823372289ULL, 1103823372544ULL, 1103823372545ULL, 1103823437824ULL, 1103823437825ULL, 1103823438080ULL, 1103823438081ULL, 281474976710656ULL, 281474976710657ULL, 281474976710912ULL, 281474976710913ULL, 281474976776192ULL, 281474976776193ULL, 281474976776448ULL, 281474976776449ULL, 281474993487872ULL, 281474993487873ULL, 281474993488128ULL, 281474993488129ULL, 281474993553408ULL, 281474993553409ULL, 281474993553664ULL, 281474993553665ULL, 281479271677952ULL, 281479271677953ULL, 281479271678208ULL, 281479271678209ULL, 281479271743488ULL, 281479271743489ULL, 281479271743744ULL, 281479271743745ULL, 281479288455168ULL, 281479288455169ULL, 281479288455424ULL, 281479288455425ULL, 281479288520704ULL, 281479288520705ULL, 281479288520960ULL, 281479288520961ULL, 282574488338432ULL, 282574488338433ULL, 282574488338688ULL, 282574488338689ULL, 282574488403968ULL, 282574488403969ULL, 282574488404224ULL, 282574488404225ULL, 282574505115648ULL, 282574505115649ULL, 282574505115904ULL, 282574505115905ULL, 282574505181184ULL, 282574505181185ULL, 282574505181440ULL, 282574505181441ULL, 282578783305728ULL, 282578783305729ULL, 282578783305984ULL, 282578783305985ULL, 282578783371264ULL, 282578783371265ULL, 282578783371520ULL, 282578783371521ULL, 282578800082944ULL, 282578800082945ULL, 282578800083200ULL, 282578800083201ULL, 282578800148480ULL, 282578800148481ULL, 282578800148736ULL, 282578800148737ULL, 72057594037927936ULL, 72057594037927937ULL, 72057594037928192ULL, 72057594037928193ULL, 72057594037993472ULL, 72057594037993473ULL, 72057594037993728ULL, 72057594037993729ULL, 72057594054705152ULL, 72057594054705153ULL, 72057594054705408ULL, 72057594054705409ULL, 72057594054770688ULL, 72057594054770689ULL, 72057594054770944ULL, 72057594054770945ULL, 72057598332895232ULL, 72057598332895233ULL, 72057598332895488ULL, 72057598332895489ULL, 72057598332960768ULL, 72057598332960769ULL, 72057598332961024ULL, 72057598332961025ULL, 72057598349672448ULL, 72057598349672449ULL, 72057598349672704ULL, 72057598349672705ULL, 72057598349737984ULL, 72057598349737985ULL, 72057598349738240ULL, 72057598349738241ULL, 72058693549555712ULL, 72058693549555713ULL, 72058693549555968ULL, 72058693549555969ULL, 72058693549621248ULL, 72058693549621249ULL, 72058693549621504ULL, 72058693549621505ULL, 72058693566332928ULL, 72058693566332929ULL, 72058693566333184ULL, 72058693566333185ULL, 72058693566398464ULL, 72058693566398465ULL, 72058693566398720ULL, 72058693566398721ULL, 72058697844523008ULL, 72058697844523009ULL, 72058697844523264ULL, 72058697844523265ULL, 72058697844588544ULL, 72058697844588545ULL, 72058697844588800ULL, 72058697844588801ULL, 72058697861300224ULL, 72058697861300225ULL, 72058697861300480ULL, 72058697861300481ULL, 72058697861365760ULL, 72058697861365761ULL, 72058697861366016ULL, 72058697861366017ULL, 72339069014638592ULL, 72339069014638593ULL, 72339069014638848ULL, 72339069014638849ULL, 72339069014704128ULL, 72339069014704129ULL, 72339069014704384ULL, 72339069014704385ULL, 72339069031415808ULL, 72339069031415809ULL, 72339069031416064ULL, 72339069031416065ULL, 72339069031481344ULL, 72339069031481345ULL, 72339069031481600ULL, 72339069031481601ULL, 72339073309605888ULL, 72339073309605889ULL, 72339073309606144ULL, 72339073309606145ULL, 72339073309671424ULL, 72339073309671425ULL, 72339073309671680ULL, 72339073309671681ULL, 72339073326383104ULL, 72339073326383105ULL, 72339073326383360ULL, 72339073326383361ULL, 72339073326448640ULL, 72339073326448641ULL, 72339073326448896ULL, 72339073326448897ULL, 72340168526266368ULL, 72340168526266369ULL, 72340168526266624ULL, 72340168526266625ULL, 72340168526331904ULL, 72340168526331905ULL, 72340168526332160ULL, 72340168526332161ULL, 72340168543043584ULL, 72340168543043585ULL, 72340168543043840ULL, 72340168543043841ULL, 72340168543109120ULL, 72340168543109121ULL, 72340168543109376ULL, 72340168543109377ULL, 72340172821233664ULL, 72340172821233665ULL, 72340172821233920ULL, 72340172821233921ULL, 72340172821299200ULL, 72340172821299201ULL, 72340172821299456ULL, 72340172821299457ULL, 72340172838010880ULL, 72340172838010881ULL, 72340172838011136ULL, 72340172838011137ULL, 72340172838076416ULL, 72340172838076417ULL, 72340172838076672ULL, 72340172838076673ULL};
    constexpr static const uint64_t MASK = ((1ULL<<(N*8)) - 1ULL) & (row_transpose[(1<<N)-1] * ((1ULL<<N) - 1));
    friend SmallBinaryMatrix reorder_principal_minors_non_singular<>(SmallBinaryMatrix);
    friend std::pair<SmallBinaryMatrix, SmallBinaryMatrix> LU_decomposition<>(SmallBinaryMatrix); 
    
private:
    uint64_t data;
};

template<uint32_t N>
SmallBinaryMatrix<N> reorder_principal_minors_non_singular(SmallBinaryMatrix<N> m) {
    assert(!m.singular());
    uint32_t row_ord[N];
    for(int i = 0; i < N; i++) row_ord[i] = i;
    // Gauss elimination
    SmallBinaryMatrix<N> tmp(m);
    for(uint32_t i = 0; i < N; i++){
        uint32_t r = tmp.first_row_with_one(i, i);
        if(r > i) {
            tmp.swap_rows(i, r);
            std::swap(row_ord[i], row_ord[r]);
        }
        auto cur_row = tmp[i];
        tmp.data ^= (((SmallBinaryMatrix<N>::row_transpose[256 - (1<<(i+1))]<<i) & tmp.data) >> i) * (uint8_t)cur_row;
    }
    SmallBinaryMatrix<N> res(m);
    res.reorder_rows(row_ord);
    return res;
}

template<uint32_t N>
std::pair<SmallBinaryMatrix<N>, SmallBinaryMatrix<N>> LU_decomposition(SmallBinaryMatrix<N> m) {
    assert(m.main_minors_non_singular());
    SmallBinaryMatrix<N> tmp(m), res;
    for(uint32_t i = 0; i < N; i++){
        uint32_t r = tmp.first_row_with_one(i, i);
        if(r > i) {
            tmp.swap_rows(i, r);
        }
        auto cur_row = tmp[i];
        res.data ^= ((SmallBinaryMatrix<N>::row_transpose[256 - (1<<(i+1))]<<i) & tmp.data);
        tmp.data ^= (((SmallBinaryMatrix<N>::row_transpose[256 - (1<<(i+1))]<<i) & tmp.data) >> i) * (uint8_t)cur_row;
        res[i][i] = 1;
    }
    return {res, tmp};
}