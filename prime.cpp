#include <bits/stdc++.h>
using namespace std;

namespace PrimeTools {
using ll = long long;
using i128 = __int128_t;

inline ll multiplyMod(ll a, ll b, ll mod) {
    return static_cast<ll>(static_cast<i128>(a) * b % mod);
}

inline ll powerMod(ll base, ll exponent, ll mod) {
    ll result = 1;
    base %= mod;
    while (exponent > 0) {
        if (exponent & 1) result = multiplyMod(result, base, mod);
        base = multiplyMod(base, base, mod);
        exponent >>= 1;
    }
    return result;
}

// Miller-Rabin deterministic for 64-bit integers.
// Returns true if n is prime, false otherwise.
inline bool isPrime(ll n) {
    if (n < 2) return false;
    if (n % 2 == 0) return n == 2;
    if (n % 3 == 0) return n == 3;

    static const vector<ll> small = {5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (ll p : small) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }

    ll d = n - 1;
    int s = 0;
    while ((d & 1) == 0) {
        ++s;
        d >>= 1;
    }

    auto witness = [&](ll a) {
        ll x = powerMod(a, d, n);
        if (x == 1 || x == n - 1) return false;

        for (int r = 1; r < s; ++r) {
            x = multiplyMod(x, x, n);
            if (x == n - 1) return false;
        }
        return true;
    };

    static const vector<ll> bases = {
        2, 325, 9375, 28178, 450775, 9780504, 179526140
    };

    for (ll a : bases) {
        if (a % n == 0) continue;
        if (witness(a)) return false;
    }
    return true;
}

// Returns all prime factors with multiplicity, e.g. 12 -> {2, 2, 3}.
inline vector<ll> primeFactors(ll n) {
    vector<ll> result;
    if (n < 2) return result;

    for (ll p = 2; p <= n / p; ++p) {
        while (n % p == 0) {
            result.push_back(p);
            n /= p;
        }
    }
    if (n > 1) result.push_back(n);
    return result;
}

// Returns prime factorization with exponents, e.g. 12 -> {{2,2},{3,1}}.
inline vector<pair<ll, int>> primePowerFactorization(ll n) {
    vector<pair<ll, int>> result;
    if (n < 2) return result;

    for (ll p = 2; p <= n / p; ++p) {
        if (n % p == 0) {
            int cnt = 0;
            while (n % p == 0) {
                ++cnt;
                n /= p;
            }
            result.push_back({p, cnt});
        }
    }
    if (n > 1) result.push_back({n, 1});
    return result;
}

// Returns the number of distinct prime factors and the total number of prime factors
// with multiplicity if needed. Example: 12 -> 2 distinct primes, 3 total factors.
inline pair<int, int> factorStats(ll n) {
    auto fac = primePowerFactorization(n);
    int distinct = fac.size();
    int totalMult = 0;
    for (auto [p, e] : fac) {
        totalMult += e;
    }
    return {distinct, totalMult};
}

// Fast integer sqrt for 64-bit values.
inline ll sqrt64(ll n) {
    if (n < 0) {
        throw invalid_argument("Integer square root is undefined for negative values");
    }
    ll x = static_cast<ll>(sqrtl((long double)n));
    while (x < n / (x + 1)) ++x;
    while (x > 0 && x > n / x) --x;
    return x;
}

// Extended Euclidean algorithm: returns gcd(a,b) and finds x,y with ax+by=gcd.
inline ll euclid(ll a, ll b, ll& x, ll& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    ll nextX, nextY;
    const ll g = euclid(b, a % b, nextX, nextY);
    x = nextY;
    y = nextX - (a / b) * nextY;
    return g;
}

// Combines x = a (mod m) and x = b (mod n).
// Supports non-coprime moduli; throws if no solution or the result modulus overflows ll.
inline ll crt(ll a, ll m, ll b, ll n) {
    if (m <= 0 || n <= 0) {
        throw invalid_argument("CRT moduli must be positive");
    }

    ll x, y;
    const ll g = euclid(m, n, x, y);
    const i128 difference = static_cast<i128>(b) - a;
    if (difference % g != 0) {
        throw domain_error("Inconsistent CRT congruences");
    }

    const i128 reducedModulus = n / g;
    const i128 modulus = static_cast<i128>(m / g) * n;
    if (modulus > numeric_limits<ll>::max()) {
        throw overflow_error("CRT result modulus does not fit in long long");
    }

    i128 multiplier = (difference / g * x) % reducedModulus;
    if (multiplier < 0) multiplier += reducedModulus;
    i128 result = (static_cast<i128>(a) + static_cast<i128>(m) * multiplier) % modulus;
    if (result < 0) result += modulus;
    return static_cast<ll>(result);
}

constexpr int LIM = 5000000;
inline vector<int> phi(LIM);

inline void calculatePhi() {
    for (int i = 0; i < LIM; ++i) {
        phi[i] = i;
    }
    for (int p = 2; p < LIM; ++p) {
        if (phi[p] == p) {
            for (int multiple = p; multiple < LIM; multiple += p) {
                phi[multiple] -= phi[multiple] / p;
            }
        }
    }
}
}

// Example usage:
// int main() {
//     using namespace PrimeTools;
//     cout << isPrime(2147483647) << "\n"; // 1
//     auto fac = primeFactors(360);
//     for (auto x : fac) cout << x << ' '; // 2 2 2 3 3 5
//     cout << "\n";
//     auto pf = primePowerFactorization(360);
//     for (auto [p, e] : pf) cout << p << '^' << e << ' '; // 2^3 3^2 5^1
//     cout << "\n";
//     return 0;
// }
