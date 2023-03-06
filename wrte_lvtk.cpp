#include "main.h"
//
void wrte_lvtk(string nam, MatrixXd& nds, MatrixXi& els, int gss, MatrixXd& mat, string& mat_aux, VectorXd& dfs_sol, VectorXi& dfs_pre, VectorXd& dfs_lds)
{
  ofstream file;
  string tmp=nam.substr(0,nam.find_last_of('.'))+".vtk";
  file.open(tmp);
  file<<"# vtk DataFile Version 2.0\n";
  file<<"Results for "<< nam<< endl;
  file<<"ASCII\n";
  file<<"DATASET UNSTRUCTURED_GRID\n";
  file<<"POINTS " << nds.rows() <<" double\n";
  nds.conservativeResizeLike(MatrixXd::Zero(nds.rows(),3));
  file<< nds << endl ;
  file<<"CELLS " << els.rows() << " " << els.rows()*(els.cols()+1) << "\n";
  for(int j = 0; j < els.rows(); j++){
    file<< els.cols() <<" "<< els(j,Eigen::all) << endl;
  }
  file<<"CELL_TYPES " << els.rows() << "\n";
  int ndfs = els.cols(); int typ;
  if(ndfs==2){
    typ=3;
  }else if(ndfs==3){
    typ=5;
  }else if(ndfs==4){
    typ=9;
  }else if(ndfs==8){
    typ=12;
  }else{
    std::cout << "ERROR" << endl;
  }
  for(int j = 0; j < els.rows(); j++){file << typ << endl;}
  file << "POINT_DATA " << nds.rows() << "\n";
  file << "VECTORS displacement double" << "\n";
  file << dfs_sol.reshaped(nds.cols(),nds.rows()).transpose() << endl ;
  file << "VECTORS prescribed integer" << "\n";
  file << dfs_pre.reshaped(nds.cols(),nds.rows()).transpose() << endl ;
  file << "VECTORS loads double" << "\n";
  file << dfs_lds.reshaped(nds.cols(),nds.rows()).transpose() << endl ;
  file.close();
}
