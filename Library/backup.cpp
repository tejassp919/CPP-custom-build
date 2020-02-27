
int counter, droot, rootchild;
vi dnum, dlow, dpar;

void dfs(int v, vvi &adj, vector<bool> &artvert)
{
    dlow[v] = dnum[v] = counter++;
    forr(x,adj[v])
    {
        if (dnum[x] == -1)
        {
            dpar[x] = v;
            if (v == droot) rootchild++;
            dfs(x,adj,artvert);
            if (dlow[x] >= dnum[v]) artvert[v] = 1;
//if(dlow[x]>dnum[u]) THIS EDGE IS BRIDGE;
            dlow[v] = min(dlow[v], dlow[x]);
        }
        else if (x != dpar[v])
            dlow[v] = min(dlow[v], dnum[x]);
    }
}

vector<bool> artpoint(int num, vvi &adj)
{
    int n=adj.size();
    counter=0;
    dnum=vi(n,-1); dlow=vi(n,-1); dpar=vi(n,0);
    vector<bool> artvert(n,0);
    fol(i,1,num+1)
    {
        if(dnum[i]!=-1) continue;
        droot=i; rootchild=0;
        dfs(i,adj,artvert);
        artvert[i]= (rootchild>1);
    }
    return artvert;
}


//Use addmax and querymax for max type implementation
//Use addmin and querymin for min type implementation

struct Line {
	lli m, c, p;
	bool operator<(const Line& o) const { return m < o.m; }
	bool operator<(lli x) const { return p < x; }
};

struct LineContainer : multiset<Line, less<>> {//less
	// (for doubles, use linf = 1/.0, div(a,b) = a/b)
	lli linf = INF;
	lli div(lli a, lli b) { // floored division
		return a / b - ((a ^ b) < 0 && a % b); }
	bool isect(iterator x, iterator y) {
		if (y == end()) { x->p = linf; return false; }
		if (x->m == y->m) x->p = x->c > y->c ? linf : -linf;
		else x->p = div(y->c - x->c, x->m - y->m);
		return x->p >= y->p;
	}
	void addmax(lli m, lli c) {
		auto z = insert({m, c, 0}), y = z++, x = y;
		while (isect(y, z)) z = erase(z);
		if (x != begin() && isect(--x, y)) isect(x, y = erase(y));
		while ((y = x) != begin() && (--x)->p >= y->p)
			isect(x, erase(y));
	}
	lli querymax(lli x) {
		if(empty()) return linf;
		auto l = *lower_bound(x);
		return l.m * x + l.c;
	}
	void addmin(lli m, lli c) { addmax(-m,-c); }
	lli querymin(lli x) { return -querymax(x); }
};


void depape(int s, vi &d, vi &p)
{
    vi m(sz(d),2);
    d.assign(sz(d),mod);
    p.assign(sz(p),-1);
    d[s]=0;
    deque<int> qu;
    qu.pb(s);
    while(!qu.empty())
    {
        int x = qu.front();
        qu.pop_front();
        m[x] = 0;
        forr(y,adj[x])
        {
            if(cap[x][y]>0&&(d[y]>(d[x]+cost[x][y])))
			{
				d[y]=d[x]+cost[x][y];
				p[y]=x;
				if(m[y]==2) qu.pb(y);
				else if(m[y]==0) qu.pf(y);
				m[y]=1;
			}
        }
    }
}

void dijkstra(int s, vi &d, vi &p)
{
    d.assign(sz(d),mod);
    p.assign(sz(p),-1);
    d[s] = 0;
    using pii = pair<int, int>;
    priority_queue<pii, vector<pii>, greater<pii>> q;
    q.push({0, s});
    while (!q.empty()) {
        int v = q.top().second;
        int d_v = q.top().first;
        q.pop();
        if (d_v != d[v])    continue;

        for (auto to : adj[v]) {
            if (cap[v][to]&&d[v] + cost[v][to] < d[to]) {
                d[to] = d[v] + cost[v][to];
                p[to] = v;
                q.push({d[to], to});
            }
        }
    }
}

const int MOD = 998244353, ROOT = 3;
vector<int> rev, roots{0, 1};
void dft(vector<int> &a) {
    int n = a.size();
    if (int(rev.size()) != n) {
        int k = __builtin_ctz(n) - 1;
        rev.resize(n);
        fol(i,0,n) rev[i] = rev[i>>1]>>1 | (i&1)<<k;
    }
    fol(i,0,n)  if (rev[i] < i) swap(a[i], a[rev[i]]);
    if (int(roots.size()) < n) {
        int k = __builtin_ctz(roots.size());
        roots.resize(n);
        while ((1 << k) < n) {
            int e = mpow(ROOT, (MOD-1)>>(k+1), MOD);
            for (int i = 1<<(k-1); i < (1<<k); ++i) {
                roots[2 * i] = roots[i];
                roots[2 * i + 1] = 1LL * roots[i] * e % MOD;
            }
            ++k;
        }
    }
    for (int k = 1; k < n; k *= 2)
        for (int i = 0; i < n; i += 2 * k)
            for (int j = 0; j < k; ++j) {
                int u = a[i + j];
                int v = 1LL * a[i + j + k] * roots[k + j] % MOD;
                a[i + j] = (u + v)%MOD;
                a[i + j + k] = (u - v + MOD)%MOD;
            }
}

void idft(vector<int> &a) {
    int n = a.size();
    reverse(a.begin() + 1, a.end());
    dft(a);
    int inv = mpow(n, MOD - 2, MOD);
    fol(i,0,n) a[i] = 1LL * a[i] * inv % MOD;
}

vector<int> operator*(vector<int> a, vector<int> b) {
    int sz = 1, tot = a.size() + b.size() - 1;
    while (sz < tot) sz <<= 1;
    a.resize(sz);   b.resize(sz);
    dft(a); dft(b);
    fol(i,0,sz) a[i] = 1LL * a[i] * b[i] % MOD;
    idft(a);
    a.resize(tot);
    return a;
}

//For max type operation change all min to max and '<' to '>'
//and initialize line vector with {0,-INF}

const int XMAX = 2e5, XMIN = -2e5; //Domain of X
typedef complex<lli> point;
vector<point> line(4 * (XMAX-XMIN),{0,INF});

lli dot(point a, point b)   {   return (conj(a) * b).real(); }
lli f(point a,  int x)      {   return dot(a, {x, 1});       }

void add_line(point nw, int v = 1, int l = XMIN, int r = XMAX) {
    int m = (l + r) / 2;
    bool lef = f(nw, l) < f(line[v], l);
    bool mid = f(nw, m) < f(line[v], m);
    if(mid) swap(line[v], nw);
    if(r - l == 1)          return;
    else if(lef != mid)     add_line(nw, 2 * v, l, m);
    else                    add_line(nw, 2 * v + 1, m, r);
}

lli get(int x, int v = 1, int l = XMIN, int r = XMAX) {
    int m = (l + r) / 2;
    if(r - l == 1)  return f(line[v], x);
    else if(x < m)  return min(f(line[v], x), get(x, 2 * v, l, m));
    else            return min(f(line[v], x), get(x, 2 * v + 1, m, r));
}

int counter;
vi dnum, dlow;
vvi SCC;
deque<int> sta;

void tarjanSCC(int v, vvi &adj, vector<bool> &vis)
{
    dlow[v] = dnum[v] = counter++;
    vis[v]=1;
    sta.pb(v);
    forr(x,adj[v])
    {
        if (dnum[x] == -1) tarjanSCC(x,adj,vis);
        if (vis[x]) dlow[v] = min(dlow[v], dlow[x]);
    }
    if(dlow[v]==dnum[v])
    {
        vi scc;
        while(sta.back()!=v)
        {
            scc.pb(sta.back());
            vis[sta.back()]=0;
            sta.ppb;
        }
        scc.pb(sta.back()); sta.ppb;
        vis[v]=0;
        SCC.pb(scc);
    }
}

void storeSCC(int num, vvi &adj)
{
    int n=adj.size();
    counter=0;
    dnum=vi(n,-1); dlow=vi(n,-1);
    vector<bool> vis(n,0);
    fol(i,1,num+1) if(dnum[i]==-1) tarjanSCC(i,adj,vis);
}

#define MAF 150010
lli tre[MAF];
lli get(int x) {
    lli res = 0;
    for (int i=x+2; i > 0; i -= (i&(-i))) res += tre[i];
    return res;
}
void add(int l, int r, lli val) {
    for (int i=l+2; i<MAF; i+=(i&(-i)))
		tre[i]+=val;
	for (int i=r+3; i<MAF; i+=(i&(-i)))
		tre[i]-=val;
}

void fft(vector<cd> &a) {
    int n = a.size();
    if (n == 1) return;
    vector<cd> a0(n / 2), a1(n / 2);
    for (int i = 0; i < n/2; i++)
    {
        a0[i] = a[2*i];
        a1[i] = a[2*i+1];
    }
    fft(a0);    fft(a1);
    long double ang = 2 * PI / n;
    cd w(1), wn(cos(ang), sin(ang));
    for (int i = 0; i < n/2 ; i++)
    {
        a[i] = a0[i] + w * a1[i];
        a[i + n/2] = a0[i] - w * a1[i];
        w *= wn;
    }
}

vector<lli> multiply(vector<lli> &a, vector<lli> &b)
{
    vector<cd> fa(all(a)), fb(all(b));
    vector<lli> result;
    int n = 1, m = sz(a)+sz(b)-1;;
    while (n < (sz(a) + sz(b))) n <<= 1;
    fa.resize(n);   fb.resize(n);   result.resize(n);
    fft(fa);    fft(fb);
    fol(i,0,n) fa[i]=conj(fa[i]*fb[i]);
    fft(fa);
    fol(i,0,n) result[i] = round((fa[i].real()/n));
    result.resize(m);
    return result;
}

int generator (int p) {
    vector<int> fact;
    int phi = p-1,  n = phi;
    for (int i=2; i*i<=n; ++i)
        if (n % i == 0) {
            fact.pb(i);
            while (n % i == 0)  n /= i;
        }
    if (n > 1)  fact.pb(n);
    for (int res=2; res<=p; ++res) {
        bool ok = true;
        for (size_t i=0; i<fact.size() && ok; ++i)
            ok &= mpow (res, phi / fact[i], p) != 1;
        if (ok)  return res;
    }
    return -1;
}

bool is_prime_once(ulli n, ulli val)
{
    ulli m=n-1; m/=2;
    while(m%2==0)
    {
        if(modpow(val,m,n)==(n-1)) return true;
        m/=2;
    }
    m=modpow(val,m,n);
    if(m==(n-1)||m==1) return true;
    else return false;
}

bool is_prime(ulli n, ulli k)
{
    if(n<=1) return false;
    else if(n<4) return true;
    bool ans=true;
    fol(i,0,min(k,n-2))
    {
        ulli val = uid<ulli>(2,n-1)(rng);
        ans=ans&is_prime_once(n,val);
    }
    return ans;
}
#pragma GCC optimize("Ofast")
#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace std;
using namespace __gnu_pbds;

#define uid uniform_int_distribution
mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
#define bpop(x) __builtin_popcountll((ulli)x)
#define blead(x) __builtin_clzll((ulli)x)
#define btrail(x) __builtin_ctzll((ulli)x)
const double PI = acos(-1);
template<class my>
using ind_set = tree<my,null_type,less<my>,rb_tree_tag,tree_order_statistics_node_update>;

#define ff first
#define ss second
#define pb push_back
#define pf push_front
#define ppb pop_back()
#define ppf pop_front()
#define all(vec) vec.begin(), vec.end()
#define fol(i,a,b) for(int i=a;i<b;i++)
#define loop(i,a,b) for(int i=a;i>=b;i--)
#define forr(x,arr) for(auto& x:arr)
#define mod 1000000007
#define INF 0x3f3f3f3f3f3f3f3f
#define EPS 1e-7
#define sz(x) (int)(x).size()

using lli   = long long;
using ulli  = unsigned long long int;
using pll   = pair<lli, lli>;
using ttt   = pair<lli, pll>;
using vttt  = vector<ttt>;
using vll   = vector<pll>;
using vl    = vector<lli>;
using vi    = vector<int>;
using vvi   = vector<vector<int>>;
using cd    = complex<long double>;

#ifdef tejasp
template<typename T>
void __p(T a) { cout << a << " "; }
template<typename T, typename F>
void __p(pair<T, F> a) { cout << "{ "; __p(a.ff); __p(a.ss); cout << "} "; }
template<typename Arg1>
void __f(const char *name, Arg1 &&arg1) {
	cout<<name<<" : ";__p(arg1); cout<<endl;
}
template<typename Arg1>
void __t(const char *name, Arg1 &&arg1)
{ cout<<name<<" : { "; for (auto p : arg1) __p(p); cout<<"}"<<endl; }
template<typename Arg1, typename ... Args>
void __f(const char *names, Arg1 &&arg1, Args &&... args) {
	int bracket=0,i=0;
	for(; ;i++)
		if(names[i]==','&&bracket==0)
			break;
		else if(names[i]=='(')
			bracket++;
		else if(names[i]==')')
			bracket--;
	cout.write(names,i)<<" : ";
	__p(arg1);  cout<<"| ";
	__f(names+i+1,args...);
}
template<typename Arg1, typename Arg2>
void __f(const char *names, Arg1 arg1[], Arg2 &&arg2){
    int i=0;
	for(; ;i++) if(names[i]==',') break;
	cout.write(names,i)<<" : { ";
	fol(i,0,arg2) __p(arg1[i]);
	cout << "} "<<endl;
}
#define trace(...) { cout<<"Line:"<<__LINE__<<" | "; __f(#__VA_ARGS__, __VA_ARGS__); }
#define cotra(...) { cout<<"Line:"<<__LINE__<<" | "; __t(#__VA_ARGS__, __VA_ARGS__); }
#else
#define endl '\n'
#define trace(...)
#define cotra(...)
#endif

void tejas_919(int kkkk)
{
    lli n, m, k, q, u, v, temp=0;

}

int main()
{
    #ifndef tejasp
        ios_base::sync_with_stdio(false);cin.tie(NULL);cout.tie(NULL);
    #endif
    cout << fixed << setprecision(10);
    int t=1;
    //cin>>t;
    fol(i,0,t) { tejas_919(i+1); }
}

int max_flow(int s, int t, int n, int req_flow)
{
    int flow=0, co=0;
    vi dist(n+2), par(n+2);
    while(flow<req_flow)
    {
        //Note: if negative edges use D'Esopo algo
        SSP_ALGO(s,dist,par);
        if(dist[t]==mod) break;
		int f=req_flow-flow,cur=t;
		while(cur!=s)
		{
			f=min(f,cap[par[cur]][cur]);
			cur=par[cur];
		}
		cur=t;
		flow+=f;
		co+=(dist[t]*f);
		while(cur!=s)
		{
			cap[par[cur]][cur]-=f;
			cap[cur][par[cur]]+=f;
			cur=par[cur];
		}
    }
    if(flow<req_flow) return -1;
    return co;
}

inline lli mpow(int a, int b, int m=mod) {
    int ans = 1;
    while (b) {
        if (b&1) ans=(1LL*ans*a)%m;
        a=(1LL*a*a)%m; b/=2;
    }
    return ans;
}

#define MAXN 200010
lli fact[MAXN], invf[MAXN];

void __prec()
{
    fact[0]=1;
    fol(i,1,MAXN) fact[i]=(fact[i-1]*i)%mod;
    fol(i,0,MAXN) invf[i]=mpow(fact[i],mod-2);
}

lli nck (lli n, lli k)
{
    if(k<0||n<k) return 0;
    return (((fact[n]*invf[n-k])%mod)*invf[k])%mod;
}

struct suffix
{
    int index;
    int rank[2];
};

int cmp(suffix &a, suffix &b)
{
    return (a.rank[0] == b.rank[0])? (a.rank[1] < b.rank[1] ?1: 0):
               (a.rank[0] < b.rank[0] ?1: 0);
}

void buildsuff(string &txt, int n, int suffarray[])
{
    suffix suff[n];
    int ind[n];

    fol(i,0,n)
    {
        suff[i].index = i;
        suff[i].rank[0] = txt[i] - 'a';
        suff[i].rank[1] = ((i+1) < n)? (txt[i + 1] - 'a'): -1;
    }
    sort(suff, suff+n, cmp);

    for (int k = 4; k < 2*n; k = k*2)
    {
        int rank = 0;
        int prev_rank = suff[0].rank[0];
        suff[0].rank[0] = rank;
        ind[suff[0].index] = 0;

        fol(i,1,n)
        {
            if (suff[i].rank[0] == prev_rank &&
                    suff[i].rank[1] == suff[i-1].rank[1])
            {
                prev_rank = suff[i].rank[0];
                suff[i].rank[0] = rank;
            }
            else
            {
                prev_rank = suff[i].rank[0];
                suff[i].rank[0] = ++rank;
            }
            ind[suff[i].index] = i;
        }

        fol(i,0,n)
        {
            int nind = suff[i].index + k/2;
            suff[i].rank[1] = (nind < n) ? suff[ind[nind]].rank[0] : -1;
        }
        sort(suff, suff+n, cmp);
    }

    fol(i,0,n)  suffixArr[i] = suff[i].index;
}

ulli multi(ulli a, ulli b, ulli modi)
{
    a%=modi;
    ulli res=0;
    while(b)
    {
        if(b&1){ res+=a; res%=modi;}
        a*=2; b/=2;
        a%=modi;
    }
    return res;
}

ulli modpow(ulli a, ulli b, ulli modi)
{
    a=a%modi;
    ulli res=1;
    while(b)
    {
        if(b&1) res=multi(res,a,modi);
        b/=2;
        a=multi(a,a,modi);
    }
    return res;
}

struct DSU
{
    vi p, siz;
    DSU(int N) {
        N+=5; siz.assign(N, 1); p.assign(N, 0);
        fol(i,0,N) p[i] = i;
    }
    int findset(int i) { return (p[i]==i) ? i : (p[i]=findset(p[i])); }
    void unite(int x, int y) {
        x=findset(x); y=findset(y);
        if (x!=y) {
            if(siz[x]<siz[y]) swap(x,y);
            p[y]=x; siz[x]+=siz[y];
        }
    }
};
