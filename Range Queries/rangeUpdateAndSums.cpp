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
template<class T> using oset = tree<T, null_type, less<T>, rb_tree_tag,tree_order_statistics_node_update>;

struct Lnode{
	ll sum,set,setAll,inc,cnt;
	Lnode(){
		sum =0;
		set = 0;
		setAll = 0;
		inc = 0;
		cnt = 0;
	}
};

tcT> class SEGL{
public:
	vector<Lnode> t;		ll n;
	ll defaultReturn = 0;
	void proc(ll ind , Lnode b){
		if(b.setAll ==1){
			t[ind].set = b.set;
			t[ind].inc = b.inc;
			t[ind].setAll = b.setAll;
		}
		else{
			t[ind].inc += b.inc;
		}
		if(t[ind].setAll == 1)
			t[ind].sum = (t[ind].set+t[ind].inc)*t[ind].cnt;
		else{
			t[ind].sum += b.inc*t[ind].cnt;
		}
	}

	void push(ll ind){
		proc(2*ind,t[ind]);
		proc(2*ind+1,t[ind]);
		t[ind].set = 0;
		t[ind].setAll = 0;
		t[ind].inc = 0;
		return;
	}

	SEGL(vector<T>& v){
		n = v.size();	
		t.resize(4*n);
		build(1,0,n-1,v);
	}

	void build(ll ind , ll tl ,ll tr ,vector<T>& v){	
		if(tl==tr){	
			t[ind].sum = v[tl];	
			t[ind].cnt = 1;
			return;
		}
		ll tm = (tl+tr)/2;	build(2*ind,tl,tm,v);	build(2*ind+1,tm+1,tr,v);	
		t[ind].sum = t[2*ind].sum + t[2*ind+1].sum;
		t[ind].cnt = t[2*ind].cnt + t[2*ind+1].cnt;
	}
	
	void update(ll ind ,ll tl , ll tr ,ll l, ll r, Lnode val){	
		if(l>r){
			return;
		}
		if(tl==l && tr==r){	
			proc(ind,val);
			return;
		}
		push(ind);
		ll tm = (tl+tr)/2;	
		update(2*ind,tl,tm,l,min(tm,r),val);	
		update(2*ind+1,tm+1,tr,max(l,tm+1),r,val);
		t[ind].sum =t[2*ind].sum + t[2*ind+1].sum;
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
		ll x,l,r;	cin>>x>>l>>r;
		--l;	--r;
		if(x==3){
			cout<<tr.query(1,0,n-1,l,r)<<"\n";
			continue;
		}
		ll val;	cin>>val;
		Lnode a;
		if(x==1){
			a.inc = val;
			tr.update(1,0,n-1,l,r,a);
		}
		if(x==2){
			a.set = val;	a.setAll  = 1;
			tr.update(1,0,n-1,l,r,a);
		}
	}
}











