#include "binary_matrix.hpp"
#include "matrix_operations.hpp"
#include "utils.hpp"
#include <iostream>
#include <bitset>
#include <random>
#include <set>
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

template<uint32_t N, uint32_t K>
bool comp(SmallBinaryMatrix<N> m1, SmallBinaryMatrix<N> m2){
    if constexpr (K <= 1) {
        return (m2.get_data() & m1.get_data()) == m1.get_data();
    } else {
        if(!comp<N, K-1>(m1, m2)) return false;
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

bool is_minimal(mat m1) {
    // cout << bitset<64>(d) << '\n';
    uint64_t d = m1.get_data();
    for(uint64_t s = (d - 1)&d; s > 0; s = ((s - 1) & d)) {
        mat m2(s);
        if(!m2.singular() && comp<N, N>(m2, m1)) return false;
    }
    return true;
}
    

// bool is_minimal(uint64_t d) {
//     // cout << bitset<64>(d) << '\n';
//     auto mask1 = get_sub_matrix_mask<N, N>(mat(d));
//     for(uint64_t s = (d - 1)&d; s > 0; s = ((s - 1) & d)) {
//         mat m2(s);
//         if(!m2.singular()) {
//             auto mask2 = get_sub_matrix_mask<N, N>(m2);
//             if((mask1 & mask2) == mask2){
//                 // cout << "counter" << endl;
//                 // cout << mask1 << endl;
//                 // cout << mask2 << endl;
//                 // print(m2);
//                 return false;
//             }
//         }
//     }
//     return true;
// }

template<uint32_t N, uint32_t K>
bool check_condition(SmallBinaryMatrix<N> m) {
    auto m1 = m.template sub_matrix<K>((1<<K) - 1, (1<<K) - 1);
    auto m2 = m.template sub_matrix<K>(((1<<K) - 1)<<(N-K), ((1<<K) - 1)<<(N-K));
    if constexpr (K > 1) {
        return check_condition<N, K-1>(m) && !m1.singular() && !m2.singular(); 
    } else {
        return !m1.singular() && !m2.singular();
    }
}

mat get_good_form(mat m){
    uint32_t ord_col[N], ord_row[N];
    for(int i = 0; i < N; i++) ord_col[i] = ord_row[i] = i;
    do {
        do {
            mat tmp = m;
            tmp.reorder_cols(ord_col);
            tmp.reorder_rows(ord_row);
            if(check_condition<N, N>(tmp)) return tmp;
        } while(next_permutation(ord_row, ord_row + N));
    } while(next_permutation(ord_col, ord_col + N));
    cerr << "Cannot find good form of the matrix!" << endl;
    return mat{};
}

bool check_mat_for_rows(mat m){
    for(int i = 0; i < N; i++){
        for(int j = i + 1; j < N; j++) {
            uint8_t r1 = m[i];
            uint8_t r2 = m[j];
            if(__builtin_popcount(r1^r2) == 1 || (__builtin_popcount(r1^r2) == 2 && (r1&r2) != min(r1, r2))) {
                return true;
            }
        }
    }
    return false;
}

int main(){

    // mat tmp;

    // for(int i = 0; i < N; i++){
    //     for(int j = 0; j < N; j++){
    //         if(j != N - i - 2) {
    //             tmp[i][j] = 1;
    //         }
    //     }
    // }


    // tmp = get_normal_form(tmp);

    // print(tmp);


    // cout << (is_minimal(tmp.get_data()) && !tmp.singular() ? "GOOD" : " BAD") << endl;

    // return 0;

    int counter1 = 0, counter2 = 0;
    int counter = 0;
    vector<uint64_t> normal_forms;

    ord_mat_generator<N> gen;

    std::set<uint64_t> used;

    do{
        mat m = gen.get_mat();
        // uint64_t d = m.get_data();

        // cout << "print" << endl;
        // print(m);

        uint64_t normal_form_data = get_normal_form(m).get_data();
        if(used.count(normal_form_data)) {
            continue;
        }

        used.insert(normal_form_data);

        counter1++;
        // cout << "GOOD" << endl;

        // auto mask1 = get_sub_matrix_mask<N, N>(m);

        // bool ok = true;
        // for(int s = (d - 1)&d; s > 0 && ok; s = ((s - 1) & d)) {
        //     // mat m2(s);
        //     // if(!m2.singular()) {
        //     //     auto mask2 = get_sub_matrix_mask<N, N>(m2);
        //     //     ok = ok && (mask1 & mask2) != mask2;
        //     // }
        // }

        if(is_minimal(m)){
            if(!check_mat_for_rows(m)) {
                counter2++;
                cout << "special: " << counter2 << '\n';
                print(m);
            }
            normal_forms.push_back(normal_form_data);
            ++counter;
            cerr << "found: " << counter << " | " << counter2 << endl;
        }
    } while(gen.nxt());

    sort(normal_forms.begin(), normal_forms.end());
    normal_forms.erase(unique(normal_forms.begin(), normal_forms.end()), normal_forms.end());

    for(int i = 0; i < (int)normal_forms.size(); i++){
        if(check_mat_for_rows(mat(normal_forms[i]))) continue;
        cout << "id: " << i << endl;
        mat m = get_good_form(mat(normal_forms[i]));
        print(m);
    }

    cout << "minimal: " << normal_forms.size() << endl;
    cout << "all: " << counter1 << '\n';

    return 0;
}