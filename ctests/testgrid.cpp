#include "../src/lib/grid.h"
#include <iostream>
#include <iomanip>

double EPS=1.0e-10;

int testUniform(){
  grid::Uniform unif;
  unif.build({0.0,1.0},8);

  std::cout<<"start testing uniform grid"<<std::endl;
  std::cout<<"grid:"<<std::endl;
  for(int i=0;i<unif.size;i++){
    std::cout<<std::setw(9)<<i<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<unif.size;i++){
    std::cout<<std::setw(9)<<unif.grid[i]<<"\t";
  }
  std::cout<<std::endl;

  std::cout<<0.0<<"\t"<<unif.floor(0.0)<<std::endl;
  std::cout<<unif.grid[0]<<"\t"<<unif.floor(unif.grid[0])<<std::endl;
  std::cout<<1.0<<"\t"<<unif.floor(1.0)<<std::endl;
  std::cout<<unif.grid[unif.size-1]<<"\t"<<unif.floor(unif.grid[unif.size-1])<<std::endl;

  for(int i=1;i<unif.size-1;i++){
    std::cout<<"around "<<unif.grid[i]<<std::endl;
    std::cout<<unif.grid[i]-EPS<<"\t"<<unif.floor(unif.grid[i]-EPS)<<"\t"<<(i-1)<<std::endl;
    std::cout<<unif.grid[i]+EPS<<"\t"<<unif.floor(unif.grid[i]+EPS)<<"\t"<<(i)<<std::endl;
  }
  std::cout<<"end testing uniform"<<std::endl;

}

int testTau(){
  double beta=40;
  grid::Tau tau;
  tau.build(beta,64,6.0/4.0);

  std::cout<<"start testing tau grid"<<std::endl;
  std::cout<<"grid:"<<std::endl;
  for(int i=0;i<tau.size;i++){
    std::cout<<std::setw(9)<<i<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<tau.size;i++){
    std::cout<<std::setw(9)<<tau.grid[i]<<"\t";
  }
  std::cout<<std::endl;

  std::cout<<0.0<<"\t"<<tau.floor(0.0)<<std::endl;
  std::cout<<tau.grid[0]<<"\t"<<tau.floor(tau.grid[0])<<std::endl;
  std::cout<<beta<<"\t"<<tau.floor(beta)<<std::endl;
  std::cout<<tau.grid[tau.size-1]<<"\t"<<tau.floor(tau.grid[tau.size-1])<<std::endl;

  for(int i=1;i<tau.size-1;i++){
    std::cout<<"around "<<tau.grid[i]<<std::endl;
    std::cout<<tau.grid[i]-EPS<<"\t"<<tau.floor(tau.grid[i]-EPS)<<"\t"<<(i-1)<<std::endl;
    std::cout<<tau.grid[i]+EPS<<"\t"<<tau.floor(tau.grid[i]+EPS)<<"\t"<<(i)<<std::endl;
  }
  std::cout<<"end testing tau"<<std::endl;
}

int testTauUL(){
  double beta=40;
  grid::TauUL tau;
  tau.build(10,4,4,0.01);

  std::cout<<"start testing tau grid"<<std::endl;
  std::cout<<"grid:"<<std::endl;
  for(int i=0;i<tau.size;i++){
    std::cout<<std::setw(9)<<i<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<tau.size;i++){
    std::cout<<std::setw(9)<<tau.grid[i]<<"\t";
  }
  std::cout<<std::endl;

  std::cout<<0.0<<"\t"<<tau.floor(0.0)<<std::endl;
  std::cout<<tau.grid[0]<<"\t"<<tau.floor(tau.grid[0])<<std::endl;
  std::cout<<beta<<"\t"<<tau.floor(beta)<<std::endl;
  std::cout<<tau.grid[tau.size-1]<<"\t"<<tau.floor(tau.grid[tau.size-1])<<std::endl;

  for(int i=1;i<tau.size-1;i++){
    std::cout<<"around "<<tau.grid[i]<<std::endl;
    std::cout<<tau.grid[i]-EPS<<"\t"<<tau.floor(tau.grid[i]-EPS)<<"\t"<<(i-1)<<std::endl;
    std::cout<<tau.grid[i]+EPS<<"\t"<<tau.floor(tau.grid[i]+EPS)<<"\t"<<(i)<<std::endl;
  }
  std::cout<<"end testing tau"<<std::endl;
}


int testFermiK(){
  double kf=1;
  double maxk=3;
  double hl=2;
  int N=16;

  grid::FermiK K;
  K.build(kf,maxk,N,hl);

  std::cout<<"start testing fermiK grid"<<std::endl;
  std::cout<<"grid:"<<std::endl;
  for(int i=0;i<K.size;i++){
    std::cout<<std::setw(9)<<i<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<K.size;i++){
    std::cout<<std::setw(9)<<K.grid[i]<<"\t";
  }
  std::cout<<std::endl;

  std::cout<<0.0<<"\t"<<K.floor(0.0)<<std::endl;
  std::cout<<K.grid[0]<<"\t"<<K.floor(K.grid[0])<<std::endl;
  std::cout<<maxk<<"\t"<<K.floor(maxk)<<std::endl;
  std::cout<<K.grid[K.size-1]<<"\t"<<K.floor(K.grid[K.size-1])<<std::endl;
  std::cout<<K.grid[K.size-1]-EPS<<"\t"<<K.floor(K.grid[K.size-1]-EPS)<<std::endl;

  for(int i=1;i<K.size-1;i++){
    std::cout<<"around "<<K.grid[i]<<std::endl;
    std::cout<<K.grid[i]-EPS<<"\t"<<K.floor(K.grid[i]-EPS)<<"\t"<<(i-1)<<std::endl;
    std::cout<<K.grid[i]+EPS<<"\t"<<K.floor(K.grid[i]+EPS)<<"\t"<<(i)<<std::endl;
  }
  std::cout<<"end testing fermi K"<<std::endl;
}

int testFermiKUL(){
  double kf=1;
  double maxk=3;
  int m=4;
  int n=4;
  double scale=0.001;

  grid::FermiKUL K;
  K.build(kf,maxk,m,n,scale);

  std::cout<<"start testing fermiK grid"<<std::endl;
  std::cout<<"grid:"<<std::endl;
  for(int i=0;i<K.size;i++){
    std::cout<<std::setw(9)<<i<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<K.size;i++){
    std::cout<<std::setw(9)<<K.grid[i]<<"\t";
  }
  std::cout<<std::endl;

  std::cout<<0.0<<"\t"<<K.floor(0.0)<<std::endl;
  std::cout<<K.grid[0]<<"\t"<<K.floor(K.grid[0])<<std::endl;
  std::cout<<maxk<<"\t"<<K.floor(maxk)<<std::endl;
  std::cout<<K.grid[K.size-1]<<"\t"<<K.floor(K.grid[K.size-1])<<std::endl;
  std::cout<<K.grid[K.size-1]-EPS<<"\t"<<K.floor(K.grid[K.size-1]-EPS)<<std::endl;

  for(int i=1;i<K.size-1;i++){
    std::cout<<"around "<<K.grid[i]<<std::endl;
    std::cout<<K.grid[i]-EPS<<"\t"<<K.floor(K.grid[i]-EPS)<<"\t"<<(i-1)<<std::endl;
    std::cout<<K.grid[i]+EPS<<"\t"<<K.floor(K.grid[i]+EPS)<<"\t"<<(i)<<std::endl;
  }
  std::cout<<"end testing fermi K"<<std::endl;
}

int testBoseK(){
  double kf=1;
  double maxk=3;
  double hl=2;
  int N=16;

  grid::BoseK K;
  K.build(kf,maxk,N,hl);

  std::cout<<"start testing boseK grid"<<std::endl;
  std::cout<<"grid:"<<std::endl;
  for(int i=0;i<K.size;i++){
    std::cout<<std::setw(9)<<i<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<K.size;i++){
    std::cout<<std::setw(9)<<K.grid[i]<<"\t";
  }
  std::cout<<std::endl;

  std::cout<<0.0<<"\t"<<K.floor(0.0)<<std::endl;
  std::cout<<K.grid[0]<<"\t"<<K.floor(K.grid[0])<<std::endl;
  std::cout<<maxk<<"\t"<<K.floor(maxk)<<std::endl;
  std::cout<<K.grid[K.size-1]<<"\t"<<K.floor(K.grid[K.size-1])<<std::endl;
  std::cout<<K.grid[K.size-1]-EPS<<"\t"<<K.floor(K.grid[K.size-1]-EPS)<<std::endl;

  for(int i=1;i<K.size-1;i++){
    std::cout<<"around "<<K.grid[i]<<std::endl;
    std::cout<<K.grid[i]-EPS<<"\t"<<K.floor(K.grid[i]-EPS)<<"\t"<<(i-1)<<std::endl;
    std::cout<<K.grid[i]+EPS<<"\t"<<K.floor(K.grid[i]+EPS)<<"\t"<<(i)<<std::endl;
  }
  std::cout<<"end testing bose K"<<std::endl;
}

int testBoseKUL(){
  double kf=1;
  double maxk=3;
  double scale=0.001;
  int m=4;
  int n=4;

  grid::BoseKUL K;
  K.build(kf,maxk,m,n,scale);

  std::cout<<"start testing boseK grid"<<std::endl;
  std::cout<<"grid:"<<std::endl;
  for(int i=0;i<K.size;i++){
    std::cout<<std::setw(9)<<i<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<K.size;i++){
    std::cout<<std::setw(9)<<K.grid[i]<<"\t";
  }
  std::cout<<std::endl;

  std::cout<<0.0<<"\t"<<K.floor(0.0)<<std::endl;
  std::cout<<K.grid[0]<<"\t"<<K.floor(K.grid[0])<<std::endl;
  std::cout<<maxk<<"\t"<<K.floor(maxk)<<std::endl;
  std::cout<<K.grid[K.size-1]<<"\t"<<K.floor(K.grid[K.size-1])<<std::endl;
  std::cout<<K.grid[K.size-1]-EPS<<"\t"<<K.floor(K.grid[K.size-1]-EPS)<<std::endl;

  for(int i=1;i<K.size-1;i++){
    std::cout<<"around "<<K.grid[i]<<std::endl;
    std::cout<<K.grid[i]-EPS<<"\t"<<K.floor(K.grid[i]-EPS)<<"\t"<<(i-1)<<std::endl;
    std::cout<<K.grid[i]+EPS<<"\t"<<K.floor(K.grid[i]+EPS)<<"\t"<<(i)<<std::endl;
  }
  std::cout<<"end testing bose K"<<std::endl;
}


int testUniLog(){
  double end=1;
  double init=0;
  int m=3;
  int n=3;
  double alpha=0.5;

  grid::UniLog ul;
  ul.init({init,end}, 0, m, n, -alpha);
  std::cout<<"start testing unilog grid"<<std::endl;
  std::cout<<"length"<<ul.length<<std::endl;
  std::cout<<"m"<<ul.m<<std::endl;
  std::cout<<"n"<<ul.n<<std::endl;
  std::cout<<"grid:"<<std::endl;
  for(int i=-1;i<(m+1)*n+1;i++){
    std::cout<<std::setw(9)<<i<<"\t";
  }
  std::cout<<std::endl;
  for(int i=-1;i<(m+1)*n+1;i++){
    std::cout<<std::setw(9)<<ul.grid(i)<<"\t";
  }
  std::cout<<std::endl;

}

bool testgrid(){
  // test log grid for tau
  //testTau();

  // testUniLog();
  // testTauUL();

  testFermiKUL();
  
  // test log grid for fermiK
  //testFermiK();

  // testBoseK();
  testBoseKUL();
  // testUniform();
  return 1;
}

int main(){
  bool result=testgrid();
  return !result;
}
