// Whole-parameter-cell certificates for marked binary zippers.
// Compile with g++ -O2 -std=c++17 -fno-fast-math -ffp-contract=off.
// Every elementary interval operation is widened by nextafter. Inputs are
// exact dyadic grid endpoints. A successful cell proves all its parameters.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

static double dn(double x){return std::nextafter(x,-INFINITY);}
static double up(double x){return std::nextafter(x,INFINITY);}
struct I{double l,h; I(double x=0):l(x),h(x){} I(double a,double b):l(a),h(b){}};
static I operator+(I a,I b){return {dn(a.l+b.l),up(a.h+b.h)};}
static I operator-(I a,I b){return {dn(a.l-b.h),up(a.h-b.l)};}
static I operator-(I a){return {-a.h,-a.l};}
static I operator*(I a,I b){
  double x[]={a.l*b.l,a.l*b.h,a.h*b.l,a.h*b.h};
  return {dn(*std::min_element(x,x+4)),up(*std::max_element(x,x+4))};
}
static I sq(I a){
  double h=std::max(std::abs(a.l),std::abs(a.h));
  double l=a.l<=0&&a.h>=0?0:std::min(std::abs(a.l),std::abs(a.h));
  return {l==0?0:dn(l*l),up(h*h)};
}
struct C{I x,y; C(I a=0,I b=0):x(a),y(b){}};
static C operator+(C a,C b){return {a.x+b.x,a.y+b.y};}
static C operator-(C a,C b){return {a.x-b.x,a.y-b.y};}
static C operator-(C a){return {-a.x,-a.y};}
static C operator*(C a,C b){return {a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x};}
static C conj(C a){return {a.x,-a.y};}
static double normhi(C a){return up(std::sqrt(up(sq(a.x).h+sq(a.y).h)));}
static double norm2lo(C a){return std::max(0.,dn(sq(a.x).l+sq(a.y).l));}
struct Map{C t,a; bool rev=false; double scale=1; unsigned depth=0;};
static C at(const Map& f,C z){return f.t+f.a*(f.rev?conj(z):z);}
static Map compose(const Map& f,const Map& g){
  return {at(f,g.t),f.a*(f.rev?conj(g.a):g.a),bool(f.rev^g.rev),up(f.scale*g.scale),f.depth+g.depth};
}
struct Pair{Map u,v;};
struct Family{std::string name; int k1,k2,e1,e2; std::string U,V,R1,R2;};
static std::vector<Family> families={
  {"DD01",0,0,0,1,"01","11","0","0"},
  {"DD10",0,0,1,0,"00","10","1","1"},
  {"DD11",0,0,1,1,"0","1","01","10"},
  {"DO10",0,1,1,0,"00","10","11","11"},
  {"DO11",0,1,1,1,"0","1","0101","1010"},
  {"OO01",1,1,0,1,"01","11","00","00"},
  {"OO10",1,1,1,0,"00","10","11","11"}
};
static Map word(const std::array<Map,2>& f,const std::string& w){
  Map out{C(0,0),C(1,0),false,1,0};
  for(char c:w)out=compose(out,f[c-'0']);
  return out;
}
static std::string bits(unsigned a,unsigned n){std::string s(n,'0');for(unsigned i=0;i<n;++i)s[n-1-i]='0'+((a>>i)&1);return s;}
struct Result{bool yes=false; unsigned visits=0, maxdepth=0; double margin=INFINITY;};
static Result check(const Family& spec,double x0,double x1,double y0,double y1,unsigned budget,unsigned maxdepth){
  C q(I(x0,x1),I(y0,y1)), b=C(1,0)-q;
  double ra=normhi(q), rb=normhi(b);
  if(!(ra<1&&rb<1))return {};
  double R=std::max(up(rb/dn(2*dn(1-ra))),up(ra/dn(2*dn(1-rb))));
  if(!std::isfinite(R))return {};
  std::array<Map,2> f={Map{spec.e1?q:C(0,0),spec.e1?-q:q,bool(spec.k1),ra,1},
    Map{spec.e2?C(1,0):q,spec.e2?-b:b,bool(spec.k2),rb,1}};
  unsigned d=spec.U.size(),m=spec.R1.size();
  std::vector<Pair> stack;stack.reserve(budget+512);
  // Initial first-level intersection, apart from its eventual endpoint seed.
  for(unsigned i=0;i<(1u<<(d-1));++i)for(unsigned j=0;j<(1u<<(d-1));++j){
    std::string u="0"+bits(i,d-1),v="1"+bits(j,d-1);
    if(u==spec.U&&v==spec.V)continue;
    stack.push_back({word(f,u),word(f,v)});
  }
  // The endpoint seed intersection equals its contracting return image if
  // every pair apart from the return pair is disjoint.
  for(unsigned i=0;i<(1u<<m);++i)for(unsigned j=0;j<(1u<<m);++j){
    std::string a=bits(i,m),b=bits(j,m);
    if(a==spec.R1&&b==spec.R2)continue;
    stack.push_back({word(f,spec.U+a),word(f,spec.V+b)});
  }
  Result out;
  while(!stack.empty()){
    Pair p=stack.back();stack.pop_back();
    if(++out.visits>budget)return out;
    out.maxdepth=std::max(out.maxdepth,std::max(p.u.depth,p.v.depth));
    C delta=at(p.u,C(.5,0))-at(p.v,C(.5,0));
    double radius=up(up(p.u.scale*R)+up(p.v.scale*R));
    double gap=dn(norm2lo(delta)-up(radius*radius));
    if(gap>0){out.margin=std::min(out.margin,gap);continue;}
    // Split only the larger enclosing ball; this is still an exact cover.
    if(p.u.scale>=p.v.scale){
      if(p.u.depth>=maxdepth)return out;
      stack.push_back({compose(p.u,f[0]),p.v});
      stack.push_back({compose(p.u,f[1]),p.v});
    }else{
      if(p.v.depth>=maxdepth)return out;
      stack.push_back({p.u,compose(p.v,f[0])});
      stack.push_back({p.u,compose(p.v,f[1])});
    }
  }
  out.yes=true;return out;
}
int main(int argc,char**argv){
  if(argc<3){std::cerr<<"Usage: certify_zipper_cells FAMILY OUT.csv [N=256] [budget=5000] [depth=30] [xmin xmax ymin ymax]\n";return 2;}
  std::string name=argv[1];auto it=std::find_if(families.begin(),families.end(),[&](auto&f){return f.name==name;});
  if(it==families.end()){std::cerr<<"Unknown family\n";return 2;}
  unsigned N=argc>3?std::stoul(argv[3]):256,budget=argc>4?std::stoul(argv[4]):5000,depth=argc>5?std::stoul(argv[5]):30;
  double xmin=0,xmax=1,ymin=0,ymax=.5;
  if(argc>9){xmin=std::stod(argv[6]);xmax=std::stod(argv[7]);ymin=std::stod(argv[8]);ymax=std::stod(argv[9]);}
  unsigned NY=argc>9?N:N/2;
  std::ofstream out(argv[2]);out<<"ix,iy,x0,x1,y0,y1,status,visits,maxdepth,minimum_squared_gap\n"<<std::setprecision(17);
  unsigned yes=0,unk=0,outside=0;uint64_t totalvisits=0;
  for(unsigned j=0;j<NY;++j){
    for(unsigned i=0;i<N;++i){
      double x0=xmin+(xmax-xmin)*i/N,x1=xmin+(xmax-xmin)*(i+1)/N;
      double y0=ymin+(ymax-ymin)*j/NY,y1=ymin+(ymax-ymin)*(j+1)/NY;
      double dx=x0<=.5&&x1>=.5?0:std::min(std::abs(x0-.5),std::abs(x1-.5));
      if(dx*dx+y0*y0>=.25){++outside;continue;}
      auto r=check(*it,x0,x1,y0,y1,budget,depth);totalvisits+=r.visits;
      if(r.yes)++yes;else++unk;
      out<<i<<','<<j<<','<<x0<<','<<x1<<','<<y0<<','<<y1<<','<<(r.yes?1:0)<<','<<r.visits<<','<<r.maxdepth<<','<<(r.yes?r.margin:0)<<'\n';
    }
    if(j%16==15)std::cerr<<name<<" row "<<j+1<<"/"<<NY<<" yes="<<yes<<" unresolved="<<unk<<" visits="<<totalvisits<<'\n';
  }
  std::cerr<<name<<" complete yes="<<yes<<" unresolved="<<unk<<" outside="<<outside<<" visits="<<totalvisits<<'\n';
}
