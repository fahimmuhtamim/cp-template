#include <bits/stdc++.h>
using namespace std;

namespace FFT {
using Complex = complex<double>;
const double PI = acos(-1.0);

// In-place FFT. The input size must be a power of two.
void transform(vector<Complex>& values, bool inverse) {
    const size_t n = values.size();

    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            swap(values[i], values[j]);
        }
    }

    for (size_t length = 2; length <= n;) {
        const double angle = 2.0 * PI / length * (inverse ? -1.0 : 1.0);
        const Complex root(cos(angle), sin(angle));

        for (size_t start = 0; start < n; start += length) {
            Complex factor(1.0, 0.0);
            for (size_t offset = 0; offset < length / 2; ++offset) {
                const Complex even = values[start + offset];
                const Complex odd = values[start + offset + length / 2] * factor;
                values[start + offset] = even + odd;
                values[start + offset + length / 2] = even - odd;
                factor *= root;
            }
        }

        if (length == n) {
            break;
        }
        length <<= 1;
    }

    if (inverse) {
        for (Complex& value : values) {
            value /= static_cast<double>(n);
        }
    }
}

// Returns the coefficients of the product of two integer polynomials.
// Floating-point FFT rounding is suitable when coefficient magnitudes are moderate.
vector<long long> multiply(const vector<long long>& a, const vector<long long>& b) {
    if (a.empty() || b.empty()) {
        return {};
    }

    const size_t resultSize = a.size() + b.size() - 1;
    size_t n = 1;
    while (n < resultSize) {
        n <<= 1;
    }

    vector<Complex> fa(n), fb(n);
    for (size_t i = 0; i < a.size(); ++i) {
        fa[i] = static_cast<double>(a[i]);
    }
    for (size_t i = 0; i < b.size(); ++i) {
        fb[i] = static_cast<double>(b[i]);
    }

    transform(fa, false);
    transform(fb, false);
    for (size_t i = 0; i < n; ++i) {
        fa[i] *= fb[i];
    }
    transform(fa, true);

    vector<long long> result(resultSize);
    for (size_t i = 0; i < resultSize; ++i) {
        result[i] = llround(fa[i].real());
    }
    return result;
}
}  // namespace FFT
