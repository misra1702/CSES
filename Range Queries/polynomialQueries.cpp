#include <bits/stdc++.h>
using namespace std;

using ll = long long;
using ld = long double;
using pll = pair<ll,ll>;
using pld = pair<ld,ld>;
using vbb = vector<bool>;
using vll = vector<ll>;
using vld = vector<ld>; 
using vss = vector<string>;
using vpll = vector<pll>; 
using vpld = vector<pld>;
 
#define tcT template<class T
#define tcTU tcT, class U
// ^ lol this makes everything look weird but I'll try it
tcT> using V = vector<T>; 
tcT, size_t SZ> using AR = array<T,SZ>; 
tcT> using PR = pair<T,T>;
 
// pairs
#define f first
#define s second
// vectors
#define all(x) begin(x), end(x)
// loops
#define FOR(i,a,b) for (ll i = (a); i < (b); ++i)
#define ROF(i,a,b) for (ll i = (b)-1; i >= (a); --i)
#define trav(a,x) for (auto& a: x)
template<class T> using pqg = priority_queue<T,vector<T>,greater<T>>;			//min heap

 
const ll MOD = 1e9+7; // 998244353;
const ll MX = 2e5+5;
const ll INF = 1e18; // not too close to LLONG_MAX
const ld PI = acos((ld)-1);
const ll dx[4] = {1,0,-1,0}, dy[4] = {0,1,0,-1}; // for every grid problem!!
map<pll,char> dirFromDiff = {{{1,0},'R'},{{0,1},'U'},{{-1,0},'L'} , {{0,-1},'D'}};
map<char,pll> diffFromDir = {{'R',{1,0}} , {'U',{0,1}}, {'L',{-1,0}},{'D',{0,-1}}};
mt19937 rng((uint32_t)chrono::steady_clock::now().time_since_epoch().count()); 
vll dijkstra(ll x, vector<vpll>& v){	ll n = v.size();	vll dis(n,INF),vis(n);	dis[x] = 0;	pqg<pll> q;	q.push({0,x});
	while(!q.empty()){	pll a = q.top();	q.pop();	ll ind = a.s;	if(vis[ind] == 1)	continue;	vis[ind] = 1;
	for(auto i:v[ind]){	if(vis[i.s]==1)	continue;	dis[i.s] = min(dis[i.s],dis[ind]+i.f);	q.push({dis[i.s],i.s});}} return dis;
}	//for graph
 
// helper funcs
constexpr int pct(int x) { return __builtin_popcount(x); } // # of bits set
constexpr int bits(int x) { return 31-__builtin_clz(x); } // floor(log2(x)) 
ll cdiv(ll a, ll b) { return a/b+((a^b)>0&&a%b); } // divide a by b rounded up
ll fdiv(ll a, ll b) { return a/b-((a^b)<0&&a%b); } // divide a by b rounded down
 
 
tcTU> T fstTrue(T lo, T hi, U f) {hi ++; assert(lo <= hi); while (lo < hi) {T mid = lo+(hi-lo)/2;
		f(mid) ? hi = mid : lo = mid+1; } return lo;}
tcTU> T lstTrue(T lo, T hi, U f) {lo --; assert(lo <= hi);while (lo < hi) {T mid = lo+(hi-lo+1)/2;	
		f(mid) ? lo = mid : hi = mid-1;} return lo;}

#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
template<class T> 
using oset = tree<T, null_type, less<T>, rb_tree_tag,tree_order_statistics_node_update>;
struct Lnode{
	ll sum,l,r,lazy,k;
	Lnode(): sum(0),l(0),r(0),lazy(0),k(0){;}
};
ll cl , cr;
tcT> class SEGL{
public:
	vector<Lnode> t;	ll n;
	ll defaultReturn = 0;
	
	void proc(ll ind, ll lz, ll x){
		ll cnt = (t[ind].r-t[ind].l+1);

		ll tsumlz = lz*(cnt*(cnt+1))/2;
		ll tsumk = x*cnt;

		t[ind].sum += tsumk+tsumlz;
		t[ind].lazy += lz;
		t[ind].k += x;
	}
	void push(ll ind){
		if(t[ind].lazy==0)	return;
		proc(2*ind,t[ind].lazy,t[ind].k);
		ll tm = (t[ind].l+t[ind].r)/2;
		ll fcnt = (tm-t[ind].l+1);
		proc(2*ind+1,t[ind].lazy,t[ind].k+fcnt*t[ind].lazy);
		t[ind].lazy = 0;	t[ind].k = 0;	
	}
	SEGL(vector<T>& v){	n = v.size();	t.resize(4*n);	build(1,0,n-1,v);}
	void build(ll ind , ll tl ,ll tr ,vector<T>& v){	
		if(tl==tr){	
			t[ind].sum = v[tl];	
			t[ind].l = tl;
			t[ind].r = tr;
			return;
		}
		ll tm = (tl+tr)/2;	build(2*ind,tl,tm,v);	build(2*ind+1,tm+1,tr,v);	
		t[ind].sum = t[2*ind].sum + t[2*ind+1].sum;
		t[ind].l = t[2*ind].l;	t[ind].r = t[2*ind+1].r;
	}
	void update(ll ind ,ll tl , ll tr ,ll l, ll r){	
		if(l>r){return;}
		if(tl==l && tr==r){	
			proc(ind,1,tl-cl);
			return;
		}	
		push(ind);	
		ll tm = (tl+tr)/2;
		update(2*ind,tl,tm,l,min(tm,r));	
		update(2*ind+1,tm+1,tr,max(l,tm+1),r);
		t[ind].sum = t[2*ind].sum+t[2*ind+1].sum;
	}
	T query(ll ind ,ll tl ,ll tr , ll l , ll r){
		if(l>r){return defaultReturn;}
		if(tl==l && tr==r){	
			return t[ind].sum;
		}
		push(ind);	
		ll tm = (tl+tr)/2;	
		T a = query(2*ind,tl,tm,l,min(tm,r));	
		T b = query(2*ind+1,tm+1,tr,max(l,tm+1),r);
		return a+b;
	}
};

int main(){
	cin.sync_with_stdio(0);     cin.tie(0);     cout.tie(0);
	ll n,q;	cin>>n>>q;
	vll v(n);
	FOR(i,0,n)
		cin>>v[i];
	SEGL<ll> tr(v);
	FOR(i,0,q){
		ll x;
		cin>>x>>cl>>cr;
		--cl;	--cr;
		if(x==1){
			tr.update(1,0,n-1,cl,cr);
		}
		else{
			cout<<tr.query(1,0,n-1,cl,cr)<<endl;
		}
	}
}


