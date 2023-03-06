//
#include <iostream>
#include <fstream>
using namespace std;
//
//Eigen header files
#include <Eigen/Dense>
#include <Eigen/Sparse>
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::ArrayXi;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::seq;
//
//json stuff from https://json.nlohmann.me/integration/  
// and/or https://github.com/nlohmann/json/releases
#include <json.hpp>
using json = nlohmann::json;
//
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat; //!!!!
typedef Eigen::Triplet<double> T;
//
void read_json(string nam, MatrixXd& nds, MatrixXi& els, int *gss, MatrixXd& mat, string& mat_aux, MatrixXd& bds, MatrixXd& lds);
//
void wrte_lvtk(string nam, MatrixXd& nds, MatrixXi& els, int gss, MatrixXd& mat, string& mat_aux, VectorXd& dfs_sol, VectorXi& dfs_pre, VectorXd& dfs_lds);
//
int smat(MatrixXd& S, ArrayXd& jac);
//
int dmat(MatrixXd& c, MatrixXd& D, ArrayXd& jac);
//
//int scg_saad(int n, SpMat A, VectorXd& b, VectorXd& x);
int scg_hond(int n, SpMat& A, VectorXd& b, VectorXd& x);
int scal (int n, SpMat& A, VectorXd& b, VectorXd& d);
//
int precon(int n, int nnz, ArrayXd& C, ArrayXi& ia, ArrayXi& iz, double alpha);
//
int pcg_hond(int n, SpMat& A, VectorXd& b, VectorXd& x);
//
int Mrhor(ArrayXd& C, int n, ArrayXi& ia, ArrayXi& iz, VectorXd& r, VectorXd& rho);
