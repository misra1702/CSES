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
class SEGL{
public:
	ll n,defaultReturn = -INF;	vll t , mark , lazy;
	ll combine(const ll& a,const ll& b){	return max(a,b);}
	SEGL(const vll& v){	n = v.size();	t.resize(4*n);	mark.resize(4*n);	lazy.resize(4*n);	build(1,0,n-1, v);}
	void build(ll ind,ll tl ,ll tr , const vll& v){	if(tl==tr){	t[ind] = v[tl];	return;	}	ll tm = (tl+tr)/2;	
		build(2*ind,tl,tm,v);	build(2*ind+1,tm+1,tr,v);	t[ind] = combine(t[2*ind],t[2*ind+1]);}
	//push defined here for addition
	void push(ll x){
		if(!mark[x])	return;	
		t[2*x] += lazy[x];	lazy[2*x] += lazy[x];	mark[2*x] =1;
		t[2*x+1] += lazy[x];	lazy[2*x+1] += lazy[x];	mark[2*x+1] = 1;	lazy[x] = 0;	mark[x] = 0;}

	void updateRange(ll ind,ll tl ,ll tr, ll l ,ll r ,ll val){if(l>r)	return;
		if(tl==l && tr==r){	t[ind] += val;	lazy[ind] += val;	mark[ind] = 1;	return ;}push(ind);	ll tm = (tl+tr)/2;	
		updateRange(2*ind,tl,tm,l,min(r,tm),val);	updateRange(2*ind+1,tm+1,tr,max(tm+1,l),r,val);
		t[ind] = combine(t[2*ind],t[2*ind+1]);}
	ll query(ll ind ,ll tl ,ll tr, ll l ,ll r){	if(l>r){	return defaultReturn;	}	if(tl==l && tr==r){	return t[ind];}
		push(ind);	ll tm = (tl+tr)/2;	return combine(query(2*ind,tl,tm,l,min(tm,r)), query(2*ind+1,tm+1,tr,max(l,tm+1),r));}
};


struct Node{
	ll sum =0 ,pref =0 ,suf =0 ,ans= 0;
	Node(){;}
	Node(ll a,ll b,ll c,ll d):sum(a),pref(b),suf(c),ans(d){
		;
	}
	Node(ll x){
		sum = x;
		// pref  = max(0ll,x);
		// suf = max(0ll,x);
		// ans = max(0ll,x);
		pref = suf = ans = x;
	}
};

class SEG{
public:
	vector<Node>t;ll n;	Node defaultReturn = Node(0);
	Node combine(const Node & a , const Node & b){ 
		Node x(0);  
		x.sum = a.sum+b.sum;
		x.pref = max(a.pref,a.sum+b.pref);
		x.suf = max(b.suf,a.suf+b.sum);
		x.ans = max(a.ans,max(b.ans,a.suf+b.pref));
		return x;
	}
	SEG(const vll& v){n = v.size();	t.resize(v.size()*4);	build(1,0,n-1,v);}
	void build(ll ind , ll tl ,ll tr ,const vll& v){if(tl==tr)	{	t[ind] = Node(v[tl]);	return;}
		ll tm = (tl+tr)/2;	build(2*ind,tl,tm,v);	build(2*ind+1,tm+1,tr,v);	t[ind] = combine(t[2*ind],t[2*ind+1]);}
	void update(ll ind ,ll tl , ll tr ,ll pos, ll val){	if(tl==tr){	t[ind] = Node(val);return;}
		ll tm = (tl+tr)/2;	if(pos<=tm)		update(2*ind,tl,tm,pos,val);	else update(2*ind+1,tm+1,tr,pos,val);
		t[ind] = combine(t[2*ind],t[2*ind+1]);}
	Node query(ll ind ,ll tl ,ll tr , ll l , ll r){if(l>r){return defaultReturn;}
		if(tl==l && tr==r){	return t[ind];}
		ll tm = (tl+tr)/2;	Node a = query(2*ind,tl,tm,l,min(tm,r));	Node b = query(2*ind+1,tm+1,tr,max(l,tm+1),r);return combine(a,b);}
};

int main(){
	ll n,q;	cin>>n>>q;
	vll v(n);
	FOR(i,0,n)	cin>>v[i];
	SEG tr(v);
	FOR(i,0,q){
		ll pos,val;	cin>>pos>>val;
		--pos;
		tr.update(1,0,n-1,pos,val);
		cout<<max(0ll,tr.query(1,0,n-1,0,n-1).ans)<<endl;
	}
}











