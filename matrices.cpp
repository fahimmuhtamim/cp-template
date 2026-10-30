#include <bits/stdc++.h>
using namespace std;

namespace MatrixAlgo {
using Matrix = vector<vector<long double>>;
constexpr long double EPS = 1e-12L;

inline size_t rows(const Matrix& a) {
    return a.size();
}

inline size_t cols(const Matrix& a) {
    return a.empty() ? 0 : a[0].size();
}

inline void validate(const Matrix& a) {
    if (a.empty()) {
        return;
    }
    const size_t columnCount = a[0].size();
    for (const auto& row : a) {
        if (row.size() != columnCount) {
            throw invalid_argument("Matrix rows must have equal lengths");
        }
    }
}

inline void validateSquare(const Matrix& a) {
    validate(a);
    if (rows(a) != cols(a)) {
        throw invalid_argument("Matrix must be square");
    }
}

inline Matrix identity(size_t n) {
    Matrix result(n, vector<long double>(n, 0));
    for (size_t i = 0; i < n; ++i) {
        result[i][i] = 1;
    }
    return result;
}

inline Matrix transpose(const Matrix& a) {
    validate(a);
    Matrix result(cols(a), vector<long double>(rows(a)));
    for (size_t i = 0; i < rows(a); ++i) {
        for (size_t j = 0; j < cols(a); ++j) {
            result[j][i] = a[i][j];
        }
    }
    return result;
}

inline Matrix add(const Matrix& a, const Matrix& b) {
    validate(a);
    validate(b);
    if (rows(a) != rows(b) || cols(a) != cols(b)) {
        throw invalid_argument("Matrix dimensions must match for addition");
    }

    Matrix result = a;
    for (size_t i = 0; i < rows(a); ++i) {
        for (size_t j = 0; j < cols(a); ++j) {
            result[i][j] += b[i][j];
        }
    }
    return result;
}

inline Matrix subtract(const Matrix& a, const Matrix& b) {
    validate(a);
    validate(b);
    if (rows(a) != rows(b) || cols(a) != cols(b)) {
        throw invalid_argument("Matrix dimensions must match for subtraction");
    }

    Matrix result = a;
    for (size_t i = 0; i < rows(a); ++i) {
        for (size_t j = 0; j < cols(a); ++j) {
            result[i][j] -= b[i][j];
        }
    }
    return result;
}

inline Matrix multiply(const Matrix& a, const Matrix& b) {
    validate(a);
    validate(b);
    if (cols(a) != rows(b)) {
        throw invalid_argument("Incompatible matrix dimensions for multiplication");
    }

    Matrix result(rows(a), vector<long double>(cols(b), 0));
    for (size_t i = 0; i < rows(a); ++i) {
        for (size_t k = 0; k < cols(a); ++k) {
            for (size_t j = 0; j < cols(b); ++j) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return result;
}

inline Matrix scale(const Matrix& a, long double factor) {
    validate(a);
    Matrix result = a;
    for (auto& row : result) {
        for (long double& value : row) {
            value *= factor;
        }
    }
    return result;
}

inline long double trace(const Matrix& a) {
    validateSquare(a);
    long double result = 0;
    for (size_t i = 0; i < rows(a); ++i) {
        result += a[i][i];
    }
    return result;
}

inline size_t matrixRank(Matrix a, long double eps = EPS) {
    validate(a);
    const size_t m = rows(a);
    const size_t n = cols(a);
    size_t pivotRow = 0;

    for (size_t col = 0; col < n && pivotRow < m; ++col) {
        size_t pivot = pivotRow;
        for (size_t row = pivotRow + 1; row < m; ++row) {
            if (fabsl(a[row][col]) > fabsl(a[pivot][col])) {
                pivot = row;
            }
        }
        if (fabsl(a[pivot][col]) <= eps) {
            continue;
        }

        swap(a[pivot], a[pivotRow]);
        for (size_t row = pivotRow + 1; row < m; ++row) {
            const long double factor = a[row][col] / a[pivotRow][col];
            for (size_t j = col; j < n; ++j) {
                a[row][j] -= factor * a[pivotRow][j];
            }
        }
        ++pivotRow;
    }
    return pivotRow;
}

inline Matrix power(Matrix base, unsigned long long exponent) {
    validateSquare(base);
    Matrix result = identity(rows(base));
    while (exponent > 0) {
        if (exponent & 1ULL) {
            result = multiply(result, base);
        }
        exponent >>= 1;
        if (exponent > 0) {
            base = multiply(base, base);
        }
    }
    return result;
}

inline long double determinant(Matrix a, long double eps = EPS) {
    validateSquare(a);
    const size_t n = rows(a);
    long double result = 1;

    for (size_t col = 0; col < n; ++col) {
        size_t pivot = col;
        for (size_t row = col + 1; row < n; ++row) {
            if (fabsl(a[row][col]) > fabsl(a[pivot][col])) {
                pivot = row;
            }
        }
        if (fabsl(a[pivot][col]) <= eps) {
            return 0;
        }
        if (pivot != col) {
            swap(a[pivot], a[col]);
            result = -result;
        }

        const long double pivotValue = a[col][col];
        result *= pivotValue;
        for (size_t row = col + 1; row < n; ++row) {
            const long double factor = a[row][col] / pivotValue;
            for (size_t j = col + 1; j < n; ++j) {
                a[row][j] -= factor * a[col][j];
            }
        }
    }
    return result;
}

inline optional<Matrix> inverse(Matrix a, long double eps = EPS) {
    validateSquare(a);
    const size_t n = rows(a);
    Matrix result = identity(n);

    for (size_t col = 0; col < n; ++col) {
        size_t pivot = col;
        for (size_t row = col + 1; row < n; ++row) {
            if (fabsl(a[row][col]) > fabsl(a[pivot][col])) {
                pivot = row;
            }
        }
        if (fabsl(a[pivot][col]) <= eps) {
            return nullopt;
        }
        swap(a[pivot], a[col]);
        swap(result[pivot], result[col]);

        const long double pivotValue = a[col][col];
        for (size_t j = 0; j < n; ++j) {
            a[col][j] /= pivotValue;
            result[col][j] /= pivotValue;
        }
        for (size_t row = 0; row < n; ++row) {
            if (row == col) {
                continue;
            }
            const long double factor = a[row][col];
            for (size_t j = 0; j < n; ++j) {
                a[row][j] -= factor * a[col][j];
                result[row][j] -= factor * result[col][j];
            }
        }
    }
    return result;
}

enum class SolveStatus {
    UNIQUE,
    INFINITE,
    NONE
};

struct SolveResult {
    SolveStatus status;
    size_t rank;
    vector<long double> particular;
    vector<vector<long double>> nullspace;
};

// Solves A*x=b. For infinitely many solutions, nullspace contains a basis
// for all homogeneous solutions; every solution is particular + a combination
// of those basis vectors.
inline SolveResult solve(const Matrix& a, const vector<long double>& b,
                         long double eps = EPS) {
    validate(a);
    if (rows(a) != b.size()) {
        throw invalid_argument("A and b dimensions do not match");
    }

    const size_t m = rows(a);
    const size_t n = cols(a);
    Matrix augmented(m, vector<long double>(n + 1));
    for (size_t i = 0; i < m; ++i) {
        copy(a[i].begin(), a[i].end(), augmented[i].begin());
        augmented[i][n] = b[i];
    }

    vector<size_t> pivotColumns;
    size_t pivotRow = 0;
    for (size_t col = 0; col < n && pivotRow < m; ++col) {
        size_t pivot = pivotRow;
        for (size_t row = pivotRow + 1; row < m; ++row) {
            if (fabsl(augmented[row][col]) > fabsl(augmented[pivot][col])) {
                pivot = row;
            }
        }
        if (fabsl(augmented[pivot][col]) <= eps) {
            continue;
        }

        swap(augmented[pivot], augmented[pivotRow]);
        const long double pivotValue = augmented[pivotRow][col];
        for (size_t j = col; j <= n; ++j) {
            augmented[pivotRow][j] /= pivotValue;
        }
        for (size_t row = 0; row < m; ++row) {
            if (row == pivotRow) {
                continue;
            }
            const long double factor = augmented[row][col];
            for (size_t j = col; j <= n; ++j) {
                augmented[row][j] -= factor * augmented[pivotRow][j];
            }
        }

        pivotColumns.push_back(col);
        ++pivotRow;
    }

    for (size_t row = 0; row < m; ++row) {
        bool allZero = true;
        for (size_t col = 0; col < n; ++col) {
            if (fabsl(augmented[row][col]) > eps) {
                allZero = false;
                break;
            }
        }
        if (allZero && fabsl(augmented[row][n]) > eps) {
            return {SolveStatus::NONE, pivotColumns.size(), {}, {}};
        }
    }

    vector<long double> particular(n, 0);
    for (size_t row = 0; row < pivotColumns.size(); ++row) {
        particular[pivotColumns[row]] = augmented[row][n];
    }

    vector<bool> isPivot(n, false);
    for (size_t col : pivotColumns) {
        isPivot[col] = true;
    }
    vector<vector<long double>> nullspace;
    for (size_t freeCol = 0; freeCol < n; ++freeCol) {
        if (isPivot[freeCol]) {
            continue;
        }
        vector<long double> basis(n, 0);
        basis[freeCol] = 1;
        for (size_t row = 0; row < pivotColumns.size(); ++row) {
            basis[pivotColumns[row]] = -augmented[row][freeCol];
        }
        nullspace.push_back(move(basis));
    }

    const SolveStatus status = nullspace.empty()
        ? SolveStatus::UNIQUE
        : SolveStatus::INFINITE;
    return {status, pivotColumns.size(), move(particular), move(nullspace)};
}
}  // namespace MatrixAlgo
