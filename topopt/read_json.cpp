//
#include "main.h"
//
void read_json(string nam, MatrixXd& nds, MatrixXi& els, int *gss, MatrixXd& mat, string& mat_aux, MatrixXd& bds, MatrixXd& lds, VectorXd& dex)
{
//
  std::ifstream f(nam);
  json data = json::parse(f);
//
  int nds_r  =  data["nodes"].size();
  int nds_c  =  data["nodes"][0].size();
  int els_r  =  data["elements"].size();
  int els_c  =  data["elements"][0].size();
  int mat_e  =  data["material"].size()-1;
  int bds_r  =  data["boundary"].size();
  int bds_c  =  data["boundary"][0].size();
  int lds_r = data["load"].size();
  int lds_c = data["load"][0].size();
//
  nds.resize(nds_r,nds_c);
  for(int j = 0; j < nds_r; j++){for(int i = 0; i < nds_c; i++){ nds(j,i) = data["nodes"][j][i];}}
  els.resize(els_r,els_c);
  for(int j = 0; j < els_r; j++){for(int i = 0; i < els_c; i++){ els(j,i) = data["elements"][j][i];}}
  mat.resize(1,mat_e);
  for(int j = 0; j < mat_e; j++){mat(0,j) = data["material"][j];}
  mat_aux = data["material"][mat_e];
  bds.resize(bds_r,bds_c);
  for(int j = 0; j < bds_r; j++){for(int i = 0; i < bds_c; i++){ bds(j,i) = data["boundary"][j][i];}}
  lds.resize(lds_r,lds_c);
  for(int j = 0; j < lds_r; j++){for(int i = 0; i < lds_c; i++){ lds(j,i) = data["load"][j][i];}}
  
//
  dex.resize(els_r);
  for(int j = 0; j < els_r; j++){dex(j) = data["design"][j];}
//
  *gss = data["gauss"];
//
}
//
