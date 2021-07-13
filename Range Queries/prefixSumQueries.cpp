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
class SEG{
public:
	vll t;ll n,defaultReturn = -INF;
	ll combine(const ll & a , const ll & b){	return max(a,b);	}
	SEG(const vll& v){n = v.size();	t.resize(v.size()*4);	build(1,0,n-1,v);}
	void build(ll ind , ll tl ,ll tr ,const vll& v){if(tl==tr)	{	t[ind] = v[tl];	return;}
		ll tm = (tl+tr)/2;	build(2*ind,tl,tm,v);	build(2*ind+1,tm+1,tr,v);	t[ind] = combine(t[2*ind],t[2*ind+1]);}
	void update(ll ind ,ll tl , ll tr ,ll pos, ll val){	if(tl==tr){	t[ind] = val;return;}
		ll tm = (tl+tr)/2;	if(pos<=tm)		update(2*ind,tl,tm,pos,val);	else update(2*ind+1,tm+1,tr,pos,val);
		t[ind] = combine(t[2*ind],t[2*ind+1]);}
	ll query(ll ind ,ll tl ,ll tr , ll l , ll r){if(l>r){return defaultReturn;}
		if(tl==l && tr==r){	return t[ind];}
		ll tm = (tl+tr)/2;	ll a = query(2*ind,tl,tm,l,min(tm,r));	ll b = query(2*ind+1,tm+1,tr,max(l,tm+1),r);return combine(a,b);}
};
class SEGL{
public:
	ll n,defaultReturn = -INF;	vll t , mark , lazy;
	ll combine(const ll& a,const ll& b){	return max(a,b);}
	SEGL(const vll& v){	n = v.size();	t.resize(4*n);	mark.resize(4*n);	lazy.resize(4*n);	build(1,0,n-1, v);}
	void build(ll ind,ll tl ,ll tr , const vll& v){	if(tl==tr){	t[ind] = v[tl];	return;	}	ll tm = (tl+tr)/2;	
		build(2*ind,tl,tm,v);	build(2*ind+1,tm+1,tr,v);	t[ind] = combine(t[2*ind],t[2*ind+1]);}
	//push defined here for addition
	void push(ll x){
		if(mark[x]==0)	return;	
		t[2*x] += lazy[x];	lazy[2*x] += lazy[x];	mark[2*x] =1;
		t[2*x+1] += lazy[x];	lazy[2*x+1] += lazy[x];	mark[2*x+1] = 1;	lazy[x] = 0;	mark[x] = 0;}

	void updateRange(ll ind,ll tl ,ll tr, ll l ,ll r ,ll val){if(l>r)	return;
		if(tl==l && tr==r){	t[ind] += val;	lazy[ind] += val;	mark[ind] = 1;	return ;}push(ind);	ll tm = (tl+tr)/2;	
		updateRange(2*ind,tl,tm,l,min(r,tm),val);	updateRange(2*ind+1,tm+1,tr,max(tm+1,l),r,val);
		t[ind] = combine(t[2*ind],t[2*ind+1]);}
	ll query(ll ind ,ll tl ,ll tr, ll l ,ll r){	if(l>r){	return defaultReturn;	}	if(tl==l && tr==r){	return t[ind];}
		push(ind);	ll tm = (tl+tr)/2;	return combine(query(2*ind,tl,tm,l,min(tm,r)), query(2*ind+1,tm+1,tr,max(l,tm+1),r));}
};
class FEN{
public: vll v;	ll n;
		FEN(const vll & t): n(t.size()){	v.resize(n);	FOR(i,0,n){	add(i,t[i]);}}
		ll sumq(ll l , ll r){	return sum(r) - sum(l-1);}
		void add(ll pos, ll val){	while(pos<n){	v[pos]+=val;	pos = pos|(pos+1);}}
private:	ll sum(ll r){ll x = 0;	while(r>=0){	x += v[r];	r = (r&(r+1))-1;}	return x;}
};

ll binpow(ll x,ll n){if(n==0)	return 1;ll temp = binpow(x,n/2);temp  =(temp*temp)%MOD;
			if(n&1) {temp = (temp*x)%MOD;}return temp;}


int main(){
	cin.sync_with_stdio(0);     cin.tie(0);     cout.tie(0);
	ll n,q;	cin>>n>>q;
	vll temp(n);
	FOR(i,0,n)	cin>>temp[i];
	vll v(n);
	FOR(i,0,n){
		v[i] = temp[i];
		v[i] += (i==0)?0:v[i-1];
	}
	SEGL tr(v);

	FOR(i,0,q){
		ll x;	cin>>x;
		if(x==2){
			ll l,r;	cin>>l>>r;	--l;	--r;
			// cout<<l<<" "<<r<<endl;
			ll ans = tr.query(1,0,n-1,l,r);
			if(l!=0)	ans -= tr.query(1,0,n-1,l-1,l-1);
			cout<<max(ans,0ll)<<endl;
		}
		else{
			ll pos,nval;	cin>>pos>>nval;	--pos;
			ll cval = temp[pos];
			ll toBeadded = nval-cval;
			// cout<<pos<<" "<<toBeadded<<endl;
			tr.updateRange(1,0,n-1,pos,n-1,toBeadded);
			temp[pos] = nval;
			// FOR(j,0,n){
			// 	cout<<tr.query(1,0,n-1,j,j)<<" ";
			// }
			// cout<<endl;
		}
	}

}











