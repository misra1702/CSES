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

class SEGY{
public:
	ll n,defaultReturn = 0;
	vll t;

	SEGY(const vll& v){
		n = v.size();	t.resize(4*n);	build(1,0,n-1,v);
	}

	void build(ll ind ,ll tl, ll tr, const vll& v){
		if(tl==tr){
			t[ind] = v[tl];
			return;
		}
		ll tm = (tl+tr)/2;
		build(2*ind,tl,tm,v);
		build(2*ind+1,tm+1,tr,v);
		t[ind] = t[2*ind]+t[2*ind+1];
	}

	void update(ll ind,ll tl ,ll tr, ll pos){
		if(tl==tr){
			t[ind] ^= 1;
			return;
		}
		ll tm = (tl+tr)/2;
		if(pos<=tm)		
			update(2*ind,tl,tm,pos);
		else
			update(2*ind+1,tm+1,tr,pos);
		t[ind] = t[2*ind]+t[2*ind+1];
	}

	ll query(ll ind,ll tl,ll tr,ll l,ll r){
		// cout<<"inside"<<" "<<ind<<" "<<tl<<" "<<tr<<" "<<l<<" "<<r<<endl;
		if(l>r){
			return defaultReturn;
		}
		if(tl==l && tr==r){
			return t[ind];
		}
		ll tm = (tl+tr)/2;
		ll a= query(2*ind,tl,tm,l,min(tm,r));
		ll b= query(2*ind+1,tm+1,tr,max(tm+1,l),r);
		return a+b;
	}
};


class SEGX{
public:
	ll n,defualtReturn = 0 ;
	vector<SEGY> t;
	SEGX(const vector<vll> & v){
		n = v.size();	
		vll temp(n);
		t.resize(4*n,SEGY(temp));
		build(1,0,n-1,v);
	}

	void combine(ll ind){
		ll segsz = t[ind].t.size();
		FOR(i,0,segsz){
			t[ind].t[i] = t[2*ind].t[i] + t[2*ind+1].t[i];
		}
	}
	void combine(ll xind,ll ind ,ll tl , ll tr,ll pos){
		if(tl==tr){
			t[xind].t[ind] = t[2*xind].t[ind] + t[2*xind+1].t[ind];
			return;
		}
		ll tm = (tl+tr)/2;
		if(pos<=tm)
			combine(xind,2*ind,tl,tm,pos);
		else
			combine(xind,2*ind+1,tm+1,tr,pos);
		t[xind].t[ind] = t[2*xind].t[ind] + t[2*xind+1].t[ind];
	}

	void build(ll ind,ll tl,ll tr,const vector<vll>&v){
		if(tl==tr){
			t[ind].build(1,0,n-1,v[tl]);
			return;
		}
		ll tm = (tl+tr)/2;
		build(2*ind,tl,tm,v);
		build(2*ind+1,tm+1,tr,v);
		combine(ind);
	}

	void update(ll ind,ll tl, ll tr, ll x , ll y ){
		if(tl==tr){
			t[ind].update(1,0,n-1,y);
			return;
		}
		ll tm = (tl+tr)/2;
		if(x<=tm)
			update(2*ind,tl,tm,x,y);
		else
			update(2*ind+1,tm+1,tr,x,y);
		combine(ind,1,0,n-1,y);

	}
	ll query(ll ind , ll tl ,ll tr,ll lx ,ll rx ,ll ly,ll ry){
		if(lx>rx){
			return defualtReturn;
		}
		if(tl==lx && tr == rx){
			return t[ind].query(1,0,n-1,ly,ry);
		}
		ll tm = (tl+tr)/2;
		ll a = query(2*ind,tl,tm,lx,min(rx,tm),ly,ry);
		ll b= query(2*ind+1,tm+1,tr,max(lx,tm+1),rx,ly,ry);
		return a+b;
	}
};


int main(){
	cin.sync_with_stdio(0);     cin.tie(0);     cout.tie(0);
	ll n,q;	cin>>n>>q;
	vector<vll> v(n,vll(n));
	FOR(i,0,n){
		FOR(j,0,n){
			char c;	cin>>c;
			if(c=='*')	v[i][j] = 1;
		}
	}
	SEGX tr(v);

	FOR(i,0,q){
		ll temp;	cin>>temp;
		if(temp==1){
			ll	x,y;	cin>>x>>y;
			--x;	--y;
			tr.update(1,0,n-1,x,y);
		}
		if(temp==2){
			ll x1,y1,x2,y2;
			cin>>x1>>y1>>x2>>y2;
			--x1;	--x2;	--y1;	--y2;
			cout<<tr.query(1,0,n-1,x1,x2,y1,y2)<<"\n";
		}
		
	}

}











