#pragma once
#include <vector>

struct IndicesResult {
    std::vector<int> cols;
    std::vector<int> rows;
    int nz_count;
};
IndicesResult computeIndices(const int N, const std::vector<int>& patt) {
    std::vector<int> cols;
    std::vector<int> rows;
    int nz_cnt = 0;

    // Estimate the maximum size for cols and reserve space
    cols.reserve(N * patt.size());
    rows.reserve(N * patt.size());  // Assuming each pattern entry could theoretically be valid

    for (int k = 1; k <= N; k++) {
        for (size_t i = 0; i < patt.size(); i++) {
            int val = k + patt[i];
            if (val > 0 && val <= N) {
                nz_cnt++;
                cols.push_back(k - 1);
                rows.push_back(val - 1);
            }
        }
    }

    return { cols, rows, nz_cnt };
}