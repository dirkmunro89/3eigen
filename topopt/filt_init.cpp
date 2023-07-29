#include "main.h"
#include <set>
//
int filt_init(int n_e, MatrixXi& els, MatrixXd& nds, SpMat& H, int lvl, int nze)
{
//
  std::vector<T> coefs;
  coefs.reserve(nze*n_e);
//
  ArrayXi con_e;
  MatrixXi noc_e;
//
//make inverse of connect matrix, matrix of centroid coordinates, and get length scale
//
  MatrixXi noc = -1*MatrixXi::Ones(nds.rows(), 9);
  ArrayXi cnt = VectorXi::Zero(nds.rows());
  MatrixXd cen = MatrixXd::Zero(n_e,3);
  VectorXd cmp;
  double len = 1e8;
  VectorXd cmp_0;
  VectorXd Hsum = VectorXd::Zero(n_e);
//
//this is not needed with my filtering:
//
  for(int e = 0; e < n_e; e++){
    con_e = els(e,Eigen::all); 
    cmp = VectorXd::Zero(3);
    for(int i = 0; i < con_e.size(); i++){
      noc(con_e(i),cnt(con_e(i))) = e;
      cmp = cmp + nds(con_e(i),Eigen::all).transpose()/con_e.size();
    }
    cen(e,Eigen::all) = cmp;
    if(e==0){
      cmp_0 = cmp;
    }else{
      len = min(len, (cmp-cmp_0).norm());
    }
    cnt(con_e) = cnt(con_e) + 1;
    if((cnt(con_e) == 9).any()){
      std::cout << "ERROR" << endl;
    }
  }
//
//for each element
//
  double krn;
  double sum_max=0.;
  int c_max=0;
  for(int e = 0; e < n_e; e++){
//
      int c_e=0;
      double sum_e=0.;
      std::set<int> q0; 
      q0.insert(e);
      krn = (float) lvl+1.0;
      coefs.push_back(T(e, e, krn));
      c_e=c_e+1;
      sum_e=sum_e+krn;
//
//     for each level
        for(int l = 0; l < lvl; l++){ 
//       make an empty set for entries in this level
          std::set<int> qn;
//       loop over entries in q0, add neighbours to qn
          for(auto it = q0.begin(); it!=q0.end();it++) 
          {
             con_e = els(*it,Eigen::all); 
             noc_e = noc(con_e,Eigen::all); 
             std::set<int> tmp{noc_e.data(), noc_e.data() + noc_e.size()  }; 
             qn.insert(tmp.begin(),tmp.end());
          }
//       erase -1 placeholder
          qn.erase(-1); 
//       erase what was in q0
          for(auto it = q0.begin(); it!=q0.end();it++)
          {
            qn.erase(*it);
          }
//       make filter coefs
          for(auto it = qn.begin(); it!=qn.end();it++)
          {
            krn = (float)lvl - l;
            coefs.push_back(T(e, *it, krn));
            sum_e=sum_e+krn;
            c_e=c_e+1;
          }
//       add everything to q0 for next level
          q0.insert(qn.begin(),qn.end());
        }
      sum_max=fmax(sum_e,sum_max);
      c_max=max(c_e,c_max);
    }
//
//  cout <<"Max. number of non-zeros associated with an element: "<<c_max<< endl;
//
    for(int e = 0; e < n_e; e++){
        Hsum(e)=sum_max;
    }
    H.setFromTriplets(coefs.begin(), coefs.end());
    H = H * Hsum.asDiagonal().inverse();
//
return 0;
}
