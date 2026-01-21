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

using ll = long long;
using pii = pair<int, int>;
using vi = vector<int>;
using vll = vector<ll>;
const int INF = 1e9;
const ll LINF = 1e18;
const int MOD = 1e9 + 7;

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

template<typename T>
class SegmentTree {
    int n;
    vector<T> tree;
    T neutral;
    
    // Merge function (Sum/Min/Max/XOR)
    T merge(T a, T b) {
        return a + b; 
    }

    void build(const vector<T>& data, int node, int start, int end) {
        if (start == end) {
            tree[node] = data[start];
        } else {
            int mid = (start + end) / 2;
            build(data, 2 * node, start, mid);
            build(data, 2 * node + 1, mid + 1, end);
            tree[node] = merge(tree[2 * node], tree[2 * node + 1]);
        }
    }

    void update(int node, int start, int end, int idx, T val) {
        if (start == end) {
            tree[node] = val;  // For add: tree[node] += val;
        } else {
            int mid = (start + end) / 2;
            if (start <= idx && idx <= mid)
                update(2 * node, start, mid, idx, val);
            else
                update(2 * node + 1, mid + 1, end, idx, val);
            tree[node] = merge(tree[2 * node], tree[2 * node + 1]);
        }
    }

    T query(int node, int start, int end, int l, int r) { // NOTE: l, r are 1-indexed
        if (r < start || end < l) return neutral;
        if (l <= start && end <= r) return tree[node];
        int mid = (start + end) / 2;
        T p1 = query(2 * node, start, mid, l, r);
        T p2 = query(2 * node + 1, mid + 1, end, l, r);
        return merge(p1, p2);
    }

public:
    SegmentTree(const vector<T>& data, T neutral_elem = 0) {
        n = sz(data);
        neutral = neutral_elem;
        tree.resize(4 * n, neutral);
        build(data, 1, 0, n - 1);
    }

    // 0-based index update
    void update(int idx, T val) {
        update(1, 0, n - 1, idx, val);
    }

    // 0-based range [l, r] query
    T query(int l, int r) {
        return query(1, 0, n - 1, l, r);
    }
};

const int MAX_VAL = 1e5 + 5;

struct WaveletTree {
    int lo, hi;
    WaveletTree *l = nullptr, *r = nullptr;
    vector<int> b;

    WaveletTree(vector<int>::iterator from, vector<int>::iterator to, int x, int y) {
        lo = x, hi = y;
        if (lo == hi || from >= to) return;
        int mid = lo + (hi - lo) / 2;
        auto f = [mid](int x) { return x <= mid; };
        b.reserve(to - from + 1);
        b.pb(0);
        for (auto it = from; it != to; it++)
            b.pb(b.back() + f(*it));
        auto pivot = stable_partition(from, to, f);
        l = new WaveletTree(from, pivot, lo, mid);
        r = new WaveletTree(pivot, to, mid + 1, hi);
    }

    int kth(int l, int r, int k) {
        if (l > r) return 0;
        if (lo == hi) return lo;
        int inLeft = b[r] - b[l - 1];
        int lb = b[l - 1]; 
        int rb = b[r];
        if (k <= inLeft) return this->l->kth(lb + 1, rb, k);
        return this->r->kth(l - lb, r - rb, k - inLeft);
    }

    int lte(int l, int r, int k) {
        if (l > r || k < lo) return 0;
        if (hi <= k) return r - l + 1;
        int lb = b[l - 1], rb = b[r];
        return this->l->lte(lb + 1, rb, k) + this->r->lte(l - lb, r - rb, k);
    }
    
    ~WaveletTree() {
        delete l; delete r;
    }
};

template <typename T>
vector<T> sos_dp(vector<T> A, int num_bits) {
    int size = 1 << num_bits;
    for (int i = 0; i < num_bits; ++i) {
        for (int mask = 0; mask < size; ++mask) {
            if (mask & (1 << i)) {
                A[mask] += A[mask ^ (1 << i)];
            }
        }
    }
    return A;
}

struct DSU {
    vector<int> p, r;
    DSU(int n = 0) { init(n); }

    void init(int n) {
        p.resize(n);
        r.assign(n, 0);
        iota(all(p), 0);
    }

    int find(int x) {
        if (p[x] == x) return x;
        return p[x] = find(p[x]);
    }

    bool unite(int a, int b) {
        a = find(a); b = find(b);
        if (a == b) return false;
        if (r[a] < r[b]) swap(a, b);
        p[b] = a;
        if (r[a] == r[b]) r[a]++;
        return true;
    }
};

void bfs(int src, vector<vector<int>>& g, vector<int>& dist) {
    queue<int> q;
    q.push(src);
    dist[src] = 0;
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : g[u]) {
            if (dist[v] == -1) {
                dist[v] = dist[u] + 1;
                q.push(v);
            }
        }
    }
}

void dfs(int u, const vector<vector<int>>& g, vector<int>& vis) {
    vis[u] = 1;
    for (int v : g[u]) {
        if (!vis[v]) {
            dfs(v, g, vis);
        }
    }
}

vector<ll> dijkstra(int s, vector<vector<pair<int,ll>>>& g) {
    vector<ll> dist(sz(g), LINF);
    priority_queue<pair<ll,int>, vector<pair<ll,int>>, greater<>> pq;
    dist[s] = 0;
    pq.push({0, s});
    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (d > dist[u]) continue;
        for (auto [v, w] : g[u]) {
            if (dist[v] > d + w) {
                dist[v] = d + w;
                pq.push({dist[v], v});
            }
        }
    }
    return dist;
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

int main() {
    fast_io();

    int t;
    cin >> t;
    while (t--) {
        int n;
        // cin >> n;
        // vector<long long> arr(n);
        // cin >> arr;

    }
    
    return 0;

}
