#include "binary_matrix.hpp"
#include "matrix_operations.hpp"
#include "utils.hpp"
#include "debug.hpp"
#include "compare.hpp"
#include "matrix_permutations.hpp"
#include <bits/stdc++.h>
#include <stdint.h>
#include <set>

template<uint32_t N>
struct ord_mat_generator{
    uint8_t rows[N];
    uint8_t rows_indep[N];
    ord_mat_generator() { for(int i = 0; i < N; i++) rows[i] = rows_indep[i] = 1<<i; }
    
    uint8_t check(uint8_t x, int n){
        int mask = x;
        for(int i = 0; i < n; i++){
            int j = __builtin_ctz(rows_indep[i]);
            if(x & (1<<j)) x ^= rows_indep[i];
            mask |= rows[i];
        }

        return __builtin_popcount(mask + 1) == 1 ? x : 0;
    }

    bool nxt(int i = N - 1){
        rows[i]++;
        while (rows[i] != (1<<N) && !(rows_indep[i] = check(rows[i], i))) ++rows[i];
        if(rows[i] == (1<<N)) rows[i] = 0;
        if(rows[i] == 0 && i == 0) return false;
        if(rows[i] == 0){
            if(!nxt(i-1)) return false;
            rows[i] = rows[i-1];
            return nxt(i);
        }
        return true;
    }

    SmallBinaryMatrix<N> get_mat(){
        SmallBinaryMatrix<N> res;
        for(int i = 0; i < N; i++) res[i] = rows[i];
        return res;
    }
};

const int N = 4;
using mat = SmallBinaryMatrix<N>;
using mat_perm = permutation_data<N>;

int main(){

    ord_mat_generator<N> gen;
    std::set<uint64_t> checked;
    std::vector<std::pair<mat_perm, mat>> minimals;
    do {
        mat m = gen.get_mat();
        uint64_t key = get_normal_form(m).get_data();
        
        if(checked.count(key)) continue;
        checked.insert(key);

        bool is_min = true;
        mat_perm cur_mp = get_permutation_data<N>(m);

        for(auto [mp, _] : minimals) {
            if(mp <= cur_mp) {
                is_min = false;
                break;
            }
        }
        if(is_min){
            for(int i = (int)minimals.size() - 1; i >= 0; i--){ // deleting smaller elements
                if(cur_mp <= minimals[i].first) {
                    minimals.erase(minimals.begin() + i);
                }
            }
            minimals.emplace_back(cur_mp, m);
        }

    } while(gen.nxt());

    std::cout << "1. ordering count: " << minimals.size() << std::endl;

    for(int i = (int)minimals.size() - 1; i >= 0; i--){
        int j = 0;
        while(j < (int)minimals.size() && (j == i || !check_prec2(minimals[j].second, minimals[i].second))) ++j;
        if(j < (int)minimals.size()) {
            minimals.erase(minimals.begin() + i);
        }
    }
    
    std::cout << "2. ordering count: " << minimals.size() << std::endl;
    
    for(int i = (int)minimals.size() - 1; i >= 0; i--){
        int j = 0;
        while(j < (int)minimals.size() && (j == i || !check_prec1(minimals[j].second, minimals[i].second))) ++j;
        if(j < (int)minimals.size()) {
            minimals.erase(minimals.begin() + i);
        }
    }

    std::cout << "3. ordering count: " << minimals.size() << std::endl;

    for(int i = 0; i < (int)minimals.size(); i++){
        std::cout << "id: " << i << std::endl;
        print(minimals[i].second);
        std::cout << '\n';
    }

    return 0;
}