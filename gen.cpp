#include "binary_matrix.hpp"
#include "matrix_operations.hpp"
#include "utils.hpp"
#include <iostream>
#include <bitset>
#include <random>
using namespace std;

constexpr const int N = 6;
using mat = SmallBinaryMatrix<N>;

template<uint32_t N>
void print(SmallBinaryMatrix<N> m){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++) {
            cout << m[i][j];
        }
        cout << '\n';
    }
}

using sub_matrix_mask = bitset<sub_matrix_cnt<N, N>()>;

template<uint32_t N>
bool check_matrix(SmallBinaryMatrix<N> m){
    for(uint32_t i = 0; i < N; i++){
        for(uint32_t j = i + 1; j < N; j++){
            uint8_t ri = (uint8_t)m[i];
            uint8_t rj = (uint8_t)m[j];
            bool subset = (ri&rj) == ri || (ri&rj) == rj;
            if(!subset && __builtin_popcount(ri^rj) == 2) return true;
            else if(subset && __builtin_popcount(ri^rj) == 1) return true;
        }
    }
    m.transpose();
    for(uint32_t i = 0; i < N; i++){
        for(uint32_t j = i + 1; j < N; j++){
            uint8_t ri = (uint8_t)m[i];
            uint8_t rj = (uint8_t)m[j];
            bool subset = (ri&rj) == ri || (ri&rj) == rj;
            if(!subset && __builtin_popcount(ri^rj) == 2) return true;
            else if(subset && __builtin_popcount(ri^rj) == 1) return true;
        }
    }
    return false;
}

template<uint32_t N>
bool check_matrix_row_col(SmallBinaryMatrix<N> m){
    for(uint32_t i = 0; i < N; i++){
        uint8_t r = (uint8_t)m[i];
        if(__builtin_popcount(r) == 1) return true; 
    }
    m.transpose();
    for(uint32_t i = 0; i < N; i++){
        uint8_t r = (uint8_t)m[i];
        if(__builtin_popcount(r) == 1) return true; 
    }
    return false;
}

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

int main(){
    int counter = 0;
    vector<pair<mat, sub_matrix_mask>> normal_forms;
    
    ord_mat_generator<N> gen;
    vector<uint64_t> norm_form_data;
    do{
        if(counter%1000 == 0) cerr << "iteration: " << (counter / 1000) << endl;
        
        counter++;
        mat m = gen.get_mat();
        if(m.singular()) {
            cerr << "wrong!!!!!!!!!!!!" << endl;
            print(m);
            cerr << "-----------------" << endl;
        }

        norm_form_data.push_back(get_normal_form(m).get_data());
    } while(gen.nxt());

    sort(norm_form_data.begin(), norm_form_data.end());
    norm_form_data.erase(unique(norm_form_data.begin(), norm_form_data.end()), norm_form_data.end());
    // cout << norm_form_data.size() << '\n';
    for(uint64_t d : norm_form_data) normal_forms.emplace_back(mat(d), get_sub_matrix_mask<N, N>(mat(d)));

    
    counter = 0;
    for(int i = 0; i < (int)normal_forms.size(); i++){
        bool fst = true;
        bool ok = true;
        for(int j = 0; j < (int)normal_forms.size(); j++){
            ok = ok && (i == j || ((normal_forms[j].second & normal_forms[i].second) != normal_forms[j].second));
        }
        if(ok){
            
            cout << "idx: " << counter << '\n';
            print(normal_forms[i].first);
            ++counter;
        }
    }
    cout << "Number of normal forms: " << norm_form_data.size() << endl;
    cout << "Number of minimal normal forms: " << counter << '\n';

    return 0;
}