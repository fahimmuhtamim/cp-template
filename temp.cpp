#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <map>
#include <set>
#include <queue>
#include <stack>
#include <deque>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <numeric>
using namespace std;

#define int long long

using ll = long long;
using pii = pair<int, int>;
using vi = vector<int>;
using vll = vector<ll>;
const int INF = 1e9;
const ll LINF = 1e18;
const int MOD = 1e9 + 7;

typedef vector<bool> vb;
typedef vector<vi> vvi;
#define FOR(i,a,b) for(int i=a;i<b;i++)
#define all(x) (x).begin(), (x).end()
#define rall(x) (x).rbegin(), (x).rend()
#define pb push_back
#define sz(x) (int)(x).size()

void fast_io() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
}

template<typename T>
istream& operator>>(istream& is, vector<T>& v) {
    for (T& x : v) is >> x;
    return is;
}

template<typename T>
ostream& operator<<(ostream& os, const vector<T>& v) {
    for (int i = 0; i < sz(v); ++i) os << v[i] << (i == sz(v) - 1 ? "" : " ");
    return os;
}

ll binpow(ll a, ll b, ll mod = MOD) {
    ll res = 1;
    a %= mod;
    while (b) {
        if (b & 1) res = res * a % mod;
        a = a * a % mod;
        b >>= 1;
    }
    return res;
}

ll inv(ll a) {
    return binpow(a, MOD - 2);
}

const int N = 2e5;
ll fact[N], invfact[N];

void init_comb() {
    fact[0] = 1;
    for (int i = 1; i < N; i++) fact[i] = fact[i-1] * i % MOD;
    invfact[N-1] = inv(fact[N-1]);
    for (int i = N-2; i >= 0; i--)
        invfact[i] = invfact[i+1] * (i+1) % MOD;
}

ll nCr(int n, int r) {
    if (r < 0 || r > n) return 0;
    return fact[n] * invfact[r] % MOD * invfact[n-r] % MOD;
}

struct StringHash {
    static const ll mod1 = 1000000007;
    static const ll mod2 = 1000000009;
    static const ll base = 911382323; // random large base

    vector<ll> h1, h2, p1, p2;

    StringHash(const string& s) {
        int n = sz(s);
        h1.resize(n + 1, 0);
        h2.resize(n + 1, 0);
        p1.resize(n + 1, 1);
        p2.resize(n + 1, 1);

        for (int i = 0; i < n; i++) {
            h1[i + 1] = (h1[i] * base + s[i]) % mod1;
            h2[i + 1] = (h2[i] * base + s[i]) % mod2;
            p1[i + 1] = (p1[i] * base) % mod1;
            p2[i + 1] = (p2[i] * base) % mod2;
        }
    }

    // 0-based, inclusive
    pair<ll,ll> get_hash(int l, int r) {
        ll x1 = (h1[r + 1] - h1[l] * p1[r - l + 1]) % mod1;
        ll x2 = (h2[r + 1] - h2[l] * p2[r - l + 1]) % mod2;
        if (x1 < 0) x1 += mod1;
        if (x2 < 0) x2 += mod2;
        return {x1, x2};
    }
};

void solve() {
    int n, m, k;
    cin >> n >> m >> k;
    vll a(n), b(m);
    cin >> a >> b;
    string s; cin >> s;
    sort(all(b));
}

int main() {
    fast_io();

    int t;
    cin >> t;
    while (t--) {
        solve();
    }
    
    return 0;
}
