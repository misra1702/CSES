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

struct Node{
	ll sum;
	Node *l , *r;
	Node(){
		sum = 0;
		l = nullptr;
		r = nullptr;
	}
};

class SEG{
public:
	ll n,defaultReturn = 0;
	Node * root ;
	ll combine(const ll & a , const ll & b){	return a+b;	}
	SEG(const vll& v){
		n = v.size();
		root = build(0,n-1,v);
	}
	Node * build(ll tl ,ll tr ,const vll& v){
		Node * a  = new Node;
		if(tl==tr)	{	
			a->sum = v[tl];
			return a;
		}
		ll tm = (tl+tr)/2;	
		a->l = build(tl,tm,v);	
		a->r = build(tm+1,tr,v);	
		a->sum = a->l->sum + a->r->sum;
		return a;
	}
	Node * update(Node * x ,ll tl , ll tr ,ll pos, ll val){	
		Node * temp = new Node;
		if(tl==tr){
			temp->sum = val;
			return temp;
		}
		ll tm = (tl+tr)/2;
		if(pos<=tm){
			temp->l = update(x->l,tl,tm,pos,val);	
			temp->r = x->r;
		}	
		else{
			temp->l = x->l;
			temp->r = update(x->r,tm+1,tr,pos,val);
		}
		temp->sum = temp->l->sum + temp->r->sum;
		return temp;
	}
	ll query(Node * x,ll tl ,ll tr , ll l , ll r){
		if(l>r){
			return defaultReturn;
		}
		if(tl==l && tr==r){	
			return x->sum;
		}
		ll tm = (tl+tr)/2;	
		ll a = query(x->l,tl,tm,l,min(tm,r));	
		ll b = query(x->r,tm+1,tr,max(l,tm+1),r);
		return a+b;
	}
};

int main(){
	cin.sync_with_stdio(0);     cin.tie(0);     cout.tie(0);
	ll n,q;	cin>>n>>q;
	vll v(n);
	FOR(i,0,n){
		cin>>v[i];
	}
	SEG tr(v);
	vector<Node*> arr;
	arr.push_back(tr.root);
	FOR(i,0,q){
		ll x;	cin>>x;
		if(x==1){
			ll k,pos,val;	cin>>k>>pos>>val;
			k--;	pos--;
			arr[k] = tr.update(arr[k],0,n-1,pos,val);
		}
		else if(x==2){
			ll k,l,r;	cin>>k>>l>>r;
			--k;	--l;	--r;
			cout<<tr.query(arr[k],0,n-1,l,r)<<endl;
		}
		else{
			ll k;	cin>>k;
			k--;
			arr.push_back(arr[k]);
		}
	}
}











