#include <RcppArmadillo.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;

// arma::uword samp_temp(const arma::vec &pvec) {
//   
//   // generate a vector of indices, so that 0 represents the largest 
//   
//   // element, 1 the second largest, and so on
//   
//   arma::uvec indx = arma::sort_index(pvec, "descend");
//   
//   // generate the q-vector, the vector of cumulative sums
//   
//   arma::vec qvec = arma::cumsum(arma::sort(pvec, "descend"));
//   
//   // draw randomly from (0,1)
//   
//   double u = arma::randu();
//   
//   // find interval into which u falls
//   
//   for (arma::uword k = 0; k < pvec.n_elem; ++k) {
//     
//     if (u <= qvec(k))
//       return indx(k);
//     
//   }
// }

// [[Rcpp::export]]

void BupC(mat &B, mat &B2inv, vec &clsno, vec &clsmem, vec &beta, vec &beta1, const vec wl, const double sig, const double kappa, const mat &Xmu, const int random_scan_n) {
  unsigned int n = B.n_cols;
  unsigned int p = Xmu.n_cols;
  uword Totalcls=clsmem.size();
  mat Bminusi(n, n-1, fill::zeros);
  mat B2invii(n-1, n-1, fill::zeros);
  mat B2i(n-1,n-1, fill::zeros);
  uvec indexi = regspace<uvec>(1,  1,  n-1);
  vec indexf = regspace<vec>(0,  1,  n);
  double temp;
  uword i=0;
  
  
  uvec random_scan_idx = randperm( n, random_scan_n);
  for(int ri=0;ri<random_scan_idx.n_elem;ri++){
    int i = random_scan_idx(ri);
    // cout<<i<<endl;
    // for(i=0; i < n; i++){
    vec Bi = B.col(i);
    if(i > 0){indexi(i-1) = 0;}
    Bminusi  = B.cols(indexi);
    B2invii = B2inv.submat(indexi, indexi);
    colvec Biminusi = Bminusi.t() * Bi;
    rowvec B2invi = B2inv.row(i);
    B2invi = B2invi.cols(indexi);
    B2i  = B2invii - B2invi.t() * B2invi/B2inv(i,i);
    //Calculation of beta to identify the partitions
    beta1 = B.col(i)- Bminusi * B2i * Biminusi;
    //subtract the mean beta = beta - ;
    beta1 = round(beta1*10000000);
    double ubeta = mean(unique(beta1));
    beta = beta1 - ubeta; 
    vec indP = indexf.elem( find(beta > 0));//== beta(0)));//== ubeta(0)) ); //
    vec indN = indexf.elem( find(beta < 0));//!= beta(0)));//== ubeta(1)) ); 
    
    //beta.elem( find(beta > 0) ).ones();
    //beta.elem( find(beta < 0) ).zeros();
    
    //uvec indN = find(beta < 0);
    //uvec indP = find(beta > 0);
    
    //cube bigmat(indN.size(), indP.size(), Totalcls, fill::zeros); 
    
    //Identify the two node associated with i-th edge
    uvec indNZ = find(Bi != 0);
    
    int flagN=0; int flagP=0;
    uword vbbar=0;
    
    //To identify the partition with root node
    if(prod(indN)==0){flagN=1;}
    if(prod(indP)==0){flagP=1;}
    
    uword other = -1;
    if(flagN==1){
      uword indred = indNZ(1);
      //uvec clsred = find(clsno == clsno(indred));
      uword indexch = clsno(indred-1)-1;
      clsmem(indexch) = clsmem(indexch)-indP.size(); //updating the number of members
      vbbar = indP.size(); //cardinality of vb
    }
    
    if(flagP==1){
      uword indred = indNZ(1);
      //uvec clsred = find(clsno == clsno(indred)); //length n
      uword indexch = clsno(indred-1)-1;
      clsmem(indexch) = clsmem(indexch)-indN.size(); //length Totalcls
      vbbar = indN.size();
    }
    
    
    double gumgen = arma::randu();
    
    uword j=0;uword k=0;
    uword selectj=0; uword selectk=0; uword selectl=0;
    double maxvalue = -99999999;//bigmat(0, 0, 0);
    double payoff = 0;
    
    uword l=0;
    for(j =0 ; j<indN.size(); j++){
      
      Rcpp::checkUserInterrupt();
      
      for(k =0 ; k<indP.size(); k++){
        for(l =0 ; l < Totalcls; l++){
          payoff = maxvalue;
          if(indN(j)*indP(k) > 0){ //To check whether either of them is a root node
            gumgen = arma::randu();
            temp = - accu(pow(Xmu.row(indN(j)) - Xmu.row(indP(k)),2))/(2*sig) - p*log(2*3.14*sig)/2;//normal density part
            if(flagN==1){
              uword tindex = indN(j)-1;
              other = clsno(tindex)-1;}
            if(flagP==1){
              uword tindex = indP(k)-1;
              other = clsno(tindex)-1;}
            if(other == l){
              temp = temp+vbbar*log(wl(l))-log(-log(gumgen));
              temp = temp-(clsmem(other)+vbbar-2)*log(clsmem(other)+vbbar)+(clsmem(other)-2)*log(clsmem(other));
              //multi = 0;
              if(clsmem(other)+vbbar-2 < 0){
                temp = temp + (clsmem(other)+vbbar-2)*log(clsmem(other)+vbbar);
              }
              if(clsmem(other)-2 < 0){
                temp = temp - (clsmem(other)-2)*log(clsmem(other)); 
              }
              //bigmat(j, k, l) = temp;
              //clsmem is the vector of length K with number of classes.
              payoff=temp;
            }
          }
          if(indN(j)*indP(k) == 0){ //To check whether either of them is a root node
            gumgen = arma::randu();
            temp = - accu(pow(Xmu.row(indN(j)) - Xmu.row(indP(k)),2))/(2*kappa*sig) - p*log(2*3.14*kappa*sig)/2;//normal density part
            if(clsmem(l)==0){
              temp = temp-log(-log(gumgen))+vbbar*log(wl(l));
              if(vbbar>2){
                temp = temp-(vbbar-2)*log(vbbar);
              }
              //bigmat(j, k, l) = temp;
              payoff=temp;
            }
          }
          if(maxvalue < payoff){
            maxvalue = payoff;
            selectj = j;
            selectk = k;
            selectl = l;}
          
          if(j+k==0){
            if(indN(j)*indP(k) > 0){
              if(other == l){
                maxvalue = payoff;
              }
            }
            if(indN(j)*indP(k) == 0){
              if(maxvalue == -99999999){
                maxvalue = payoff;
              } 
            }
          }
        }
      }
      //sum1 = sum1 + bigmat(j, k);
    }
    
    
    ////vec u = randu<vec>(5);//-R::rexp(1.0);
    //u = exp(u);
    //temp= 0;
    
    //double indN.elem(k); double indP.elem(k);
    
    //// vec bigvec = vectorise(bigmat);
    // 
    // uword ind = samp_temp(bigvec);
    // 
    // int kk=0;int jj=0; uword count = -1;
    ////uword indP.elem(k) = indP.elem(k);
    ////uword indN.elem(j) = indN.elem(j);
    
    
    colvec veci = zeros(n+1);
    
    //int kk = 1; int jj=1;
    
    veci(indP(selectk)) = 1;
    veci(indN(selectj)) = -1;
    
    B.col(i) = veci;
    clsmem(selectl) = clsmem(selectl) + vbbar; //Update with selected class
    // beta1=arma::conv_to<arma::vec>::from(indP);
    uword tindex;
    if(flagN==1){
      for(uword ll =0 ; ll < indP.size(); ll++){
        tindex = indP(ll)-1;
        clsno(tindex)=selectl+1;
      }
    }
    if(flagP==1){
      for(uword ll = 0 ; ll < indN.size(); ll++){
        tindex = indN(ll)-1;
        clsno(tindex)=selectl+1;
      }
    }
    
    //B(indP.elem(k), i) = 1;
    //B(indN.elem(j), i) = -1;
    
    //Updating B2inv
    Biminusi = Bminusi.t() * veci;
    double st1 = accu(veci%veci);
    vec    st21 = B2i*Biminusi;
    double st2 = accu(Biminusi%st21);
    colvec st4 = zeros(n);
    st4(i) = 1/(st1-st2);
    colvec st3 = - (B2i*Biminusi);
    vec temp1 = st3*st4(i);
    st4(indexi) = temp1;
    B2inv.row(i) = st4.t();
    B2inv.col(i) = st4;
    B2inv(indexi,indexi) = B2i+temp1*temp1.t()/B2inv(i,i);
    
    
    if(i > 0){indexi(i-1) = i;}
  }
}