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
void read_json(string nam, MatrixXd& nds, MatrixXi& els, int *gss, MatrixXd& mat, string& mat_aux, MatrixXd& bds, MatrixXd& lds, VectorXd& dex);
//
void wrte_lvtk(string nam, int itr, MatrixXd& nds, MatrixXi& els, int gss, MatrixXd& mat, string& mat_aux, VectorXd& dfs_sol, VectorXi& dfs_pre, VectorXd& dfs_lds, VectorXd& opt_ex, VectorXd& opt_ro);
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
//
int filt_init(int n_e, MatrixXi& els, MatrixXd& nds, SpMat& H);
//
int assy(int n_e, int n_f, int nnz, MatrixXi& els, MatrixXd& nds, VectorXi& dfs_pre, VectorXi& map_num, VectorXd& opt_ro, SpMat& ffs_K);
//
int solv(int sol, int n_f, SpMat ffs_K, VectorXd& ffs_lds, VectorXd& ffs_sol);
//
int sens(int n_d, MatrixXi& els, MatrixXd& nds, VectorXd& dfs_sol, VectorXd& opt_ro, VectorXd& opt_sns, double *opt_obj);
//
int dqpsub(VectorXd& x_k, VectorXd& g, MatrixXd& dg, MatrixXd& dc, VectorXd& dx_l, VectorXd& dx_u, VectorXi& c_s, VectorXd& x_d, VectorXd& x);
