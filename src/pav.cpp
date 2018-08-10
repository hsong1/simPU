#include "RcppEigen.h"
#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
using namespace Rcpp;

class blocks{
public:
  int nblocks;
  vector<int> sidx;
  vector<int> size;
  
  blocks();
  blocks(int nblocks_, vector<int> sidx_, vector<int>size_)
    :nblocks(nblocks_),sidx(sidx_),size(size_){};
  void init_blocks(int nblocks_);
  void pool_blocks(int k);
  void split_blocks(int k, int l);
  void print();
};

blocks::blocks(){
  nblocks=0; sidx.resize(0);size.resize(0);
};

//initialize blocks such that each block is of size 1
void blocks::init_blocks(int nblocks_)
{
  nblocks=nblocks_;
  sidx.resize(nblocks); size.resize(nblocks);
  for(int i=0;i<nblocks;i++){
    sidx[i]=i; size[i]=1;
  }
}

//pool block k and k+1
void blocks::pool_blocks(int k){
  size[k] += size[k+1];
  sidx.erase(sidx.begin()+(k+1));
  size.erase(size.begin()+(k+1));
  nblocks-=1;
}

// split block k at lth position
// for example if l=1, we split blocks into (k),(k+1,k+2,...,)
// for example if l=2, we split blocks into (k,k+1),(k+2,...,)
// for example if l=3, we split blocks into (k,k+1,k+2), (k+3,k+4,..)
// l should be less than the block size
void blocks::split_blocks(int k, int l){
  
  sidx.insert(sidx.begin()+(k+1), sidx[k]+l); // sidx of new block = sidx[k]+l
  size.insert(size.begin()+(k+1), size[k]-l); // size of new block = size[k]-l
  size[k] = l; //size of block k =l
  nblocks+=1;
  
}

void blocks::print(){
  for (int i=0;i<nblocks;i++){
    cout<<"block "<<i<<":\n sidx:"<<sidx[i]<<", size"<<size[i]<<"\n";
  }
}

class pav {
private:
  VectorXd y;
  VectorXd ycum;
  VectorXd yfit;
  blocks b;
public:
  pav(VectorXd & y_);
  pav(VectorXd & y_, blocks b_);
  double block_mean(int k);
  bool adjacent_violators(int k);
  void pava_pre();
  void pava();
  VectorXd getFitted(){return yfit;}
  blocks getblocks(){return b;}
};

// pav constructor
pav::pav(VectorXd & y_)
  :y(y_),ycum(y_),yfit(y_){
  // Each observation == a block
  b.init_blocks(static_cast<int>(y.size()));
  // Save cumulative sum of y needed for future calculation
  for(int i=1;i<y.size();i++){
    ycum[i] = ycum[i-1]+y[i];
  }
}

// pav constructor--when block is provided
pav::pav(VectorXd & y_, blocks b_)
  :y(y_),ycum(y_),yfit(y_),b(b_){
  // Save cumulative sum of y needed for future calculation
  for(int i=1;i<y.size();i++){
    ycum[i] = ycum[i-1]+y[i];
  }
}

// return mean of y in block k
double pav::block_mean(int k){
  double meank(0);
  if(k>0){
    meank = (ycum(b.sidx[k]+b.size[k]-1)-ycum(b.sidx[k]-1))/b.size[k];
    
  }else{
    meank = ycum(b.sidx[k]+b.size[k]-1)/b.size[k];
  }
  //    cout<<"meank"<<meank<<"ref"<<y.segment(b.sidx[k], b.size[k]).array().mean()<<endl;
  return meank;
}

// return whether block k, k+1 are adjacent violators
bool pav::adjacent_violators(int k){
  double meank(0), meank1(0);
  meank = block_mean(k);
  meank1 = block_mean(k+1);
  return meank>meank1;
}

void pav::pava_pre(){
  int k(0), i(0); double meank; bool violation(false);
  VectorXd ycum_segmean;
  //For each block,
  while(k<b.nblocks){
    meank = block_mean(k);
    if(k>0){
      ycum_segmean = ycum.segment(b.sidx[k],b.size[k]).array() - ycum(b.sidx[k]-1);
    }else{
      ycum_segmean = ycum.segment(b.sidx[k],b.size[k]);
    }
    
    i=0; violation=false;
    while(!violation && i<ycum_segmean.size()){
      violation = (ycum_segmean(i)/(i+1)< meank);
      i+=1;
    }
    if(violation){
      b.split_blocks(k,i);
    }else{
      
      yfit.segment(b.sidx[k], b.size[k]).setConstant(block_mean(k));
      k+=1;
    }
  }
}
void pav::pava(){
  // cout<<"Execute PAVA"<<endl;
  int k(0);
  //block = (0,1,..,nblocks-1), k iterates up to nblocks-2
  //If k,k+1 adjacent_violators, pool two blocks and restart the proceduere
  //Otherwise move on to the next block
  bool flag(true);
  while(flag){
    
    flag = false; k=0;
    while(k<(b.nblocks-1)){
      // cout<<"adjacent_violators at "<<k<<":= "<<adjacent_violators(k)<<endl;
      if(adjacent_violators(k)){
        
        b.pool_blocks(k);
        yfit.segment(b.sidx[k], b.size[k]).setConstant(block_mean(k));
        flag = true;
      }else{
        k+=1;
      }
    }
  }
  
}
// #'@export
// [[Rcpp::export]]
List pava_main(Eigen::VectorXd & y_, int nblocks, std::vector<int> sidx, 
               std::vector<int> size, bool warmstart){
  blocks blk(nblocks,sidx,size);
  if(warmstart){
    
    pav pav1(y_,blk);
    pav1.pava_pre();
    pav1.pava();
    
    int nblocks = pav1.getblocks().nblocks;
    std::vector<int> sidx = pav1.getblocks().sidx;
    std::vector<int> size = pav1.getblocks().size;
    
    return Rcpp::List::create(Rcpp::Named("fit")=pav1.getFitted(),
                              Rcpp::Named("nblocks") = nblocks,
                              Rcpp::Named("sidx")= sidx,
                              Rcpp::Named("size")=size);
  }else{
    pav pav1(y_);
    pav1.pava();

    int nblocks = pav1.getblocks().nblocks;
    std::vector<int> sidx = pav1.getblocks().sidx;
    std::vector<int> size = pav1.getblocks().size;

    return Rcpp::List::create(Rcpp::Named("fit")=pav1.getFitted(),
                              Rcpp::Named("nblocks") = nblocks,
                              Rcpp::Named("sidx")= sidx,
                              Rcpp::Named("size")=size);

  }
  
}



