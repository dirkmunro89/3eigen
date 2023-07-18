#include "main.h"
#include <chrono>
//
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;
//
int main(int argc, char *argv[])
{
//
// some initializations
//
    int e, i, i1, i2, j, j1, j2, k, c, err, gss, etp;
    int nnz, nnz_c, nnz_s, n_d, n_e, n_n, n_p, n_f;
    int gind1, gind2, lind1, lind2;
    MatrixXd nds; MatrixXi els; 
    MatrixXd mat; string mat_aux;
    MatrixXd bds; MatrixXd lds;
//  ArrayXi con_e; MatrixXd nds_e; MatrixXd kay_e;
    VectorXd dfs_lds;
    VectorXi dfs_pre; 
    VectorXd dfs_sol; 
//  VectorXd dfs_bfs;
    VectorXd ffs_lds;                   
//  VectorXd ffs_bfs;                   
    VectorXd ffs_sol; 
    VectorXi map_num; 
//  std::vector<T> coefs;
//
//  for opt
//
    VectorXd opt_ex;
    VectorXd opt_la;
    VectorXd opt_ro;
    VectorXd opt_ef;
    MatrixXd opt_df;
    MatrixXd opt_cf;
    VectorXd opt_xl;
    VectorXd opt_xu;
    VectorXd opt_mv;
    VectorXi opt_cs;
//
    auto t0 = high_resolution_clock::now();
//
    if (argc > 4) {
        std::cout << "Input: " << argv[1] << " " << argv[2] << " " << argv[3] << argv[4] << endl;
    }else{
        return 1;
    }
    int thr; std::stringstream str1(argv[1]); str1 >> thr;
    int sol; std::stringstream str2(argv[2]); str2 >> sol;
    int fdc; std::stringstream str3(argv[4]); str3 >> fdc;
//
    Eigen::setNbThreads(thr);
//
// read input deck (json currently)
//
    auto t1 = high_resolution_clock::now();
    read_json(argv[3], nds, els, &gss, mat, mat_aux, bds, lds, opt_ex);
    auto t2 = high_resolution_clock::now(); 
    auto dur_inp = duration_cast<milliseconds>(t2 - t1).count();
    std::cout << "Input data parsed (" << dur_inp << ")\n";
//
// reserve estimate of number of non-zeros
//  ... and define some other things ...
//
    nnz = pow(els.cols()*nds.cols(),2)*els.rows();
//  coefs.reserve(nnz);
    n_e = els.rows(); // number of elements
    n_n = nds.rows(); // number of nodes
    n_d = nds.cols()*nds.rows(); // number of dofs
    n_p = bds.rows(); // number of prescribed dofs 
    n_f = n_d-n_p; // free dofs
    std::cout << "- Elements : " << n_e << "\n";
    std::cout << "- Nodes : " << n_n << "\n";
    std::cout << "- Dofs : " << n_d << "\n";
    std::cout << "- Free : " << n_f << "\n";
//
//  for opt
//
    opt_ro = VectorXd::Zero(n_e);
    opt_ef = VectorXd::Zero(2);
    opt_df = MatrixXd::Zero(2,n_e);
    opt_cf = MatrixXd::Zero(2, n_e);
    opt_xl = VectorXd::Ones(n_e)*1e-3;
    opt_xu = VectorXd::Ones(n_e);
    opt_mv = VectorXd::Ones(n_e)*0.1;
    opt_cs = VectorXi::Ones(1);
    opt_la = VectorXd::Zero(1);
    SpMat opt_H(n_e,n_e);
//
//  filter initialisation
//
    err=filt_init(n_e, els, nds, opt_H);
//
// DOF vectors containing 
//  - loads
//  - prescribed dofs, and 
//  - (eventually) the solution vector
//
    SpMat ffs_K(n_f,n_f);
    dfs_lds.setZero(n_d); 
    dfs_pre.setZero(n_d); 
    dfs_sol.setZero(n_d);
//  dfs_bfs.setZero(n_d);
    dfs_lds(lds.col(0)) = lds.col(1);
    dfs_pre(bds.col(0)) = Eigen::VectorXi::Ones(n_p); 
//
// Free DOF vectors for loads and solution vector
    ffs_lds.setZero(n_f); 
//  ffs_bfs.setZero(n_f); 
    ffs_sol.setZero(n_f);
//
// mapping dof numbering vector with prescribed dofs removed
//
    c=0;
    map_num.setZero(n_d);
    for(i  = 0; i < n_d ; i++){
        if(dfs_pre(i)){
            map_num(i) = -1;
        }else{
            map_num(i) = c;
            c=c+1;
        }
    }
//
//  opt loop
//
    double opt_scl;
    double opt_obj;
    VectorXd opt_tmp; 
    VectorXd opt_sns; 
//
//  finite differences
//
    int fdi = 1;
    if(fdc){
        opt_ex = 0.5*VectorXd::Ones(n_e);
    }
//
    for(int itr=0; itr<100; itr++){
//
//  filter the design variables
//
        opt_ro = opt_H*opt_ex;
//
//  free dof load vector
//
        for(i  = 0; i < n_d ; i++){
            if(dfs_pre(i)){
            }else{
                ffs_lds(map_num(i)) = dfs_lds(i);
            }
        }
//
//  assemble with filtered density in the loop
//
        err=assy(n_e,n_f,nnz,els,nds,dfs_pre,map_num,opt_ro,ffs_K,ffs_lds);
//
//  solve
//
        err=solv(sol, n_f, ffs_K, ffs_lds, ffs_sol);
//
//  pack solution vector into complte dof vector
//
        c = 0;
        for(i  = 0; i < n_d ; i++){
            if(dfs_pre(i)){
                dfs_sol(i)=0.;
//              dfs_bfs(i)=0.;
            }else{
                dfs_sol(i) = ffs_sol(c);
//              dfs_bfs(i) = ffs_bfs(c);
                c=c+1;
            }
        }
//
//  compute compliance objective and sensitivities (similar to assembly)
//
        opt_obj = 0.;
        opt_tmp = VectorXd::Zero(n_e);
        opt_sns = VectorXd::Zero(n_e);
        err=sens(n_d,els,nds,dfs_sol,dfs_pre,opt_ro,opt_tmp,&opt_obj);
//
//  filter the sensitivities (need to check if transpose multiplication is fast...)
//
        opt_sns = opt_H.transpose()*opt_tmp;
//
        if(itr==0){
            opt_scl = opt_obj;
        }
        opt_obj=opt_obj/opt_scl;
        opt_sns=opt_sns/opt_scl;
//
//  construct subproblem
//
        double opt_vol = opt_ex.sum()/n_e/1. - 1.;
//
        if(fdc){
            std::cout << "Objective " << std::setprecision(15) << opt_obj << endl;
        }else{
            std::cout << "===================================================" << "\n";
            std::cout << "Volume fraction " << opt_ex.sum()/n_e << endl;
            std::cout << "Compliance " << opt_obj << endl;
            std::cout << "===================================================" << "\n";
        }
//
        opt_ef(0) = opt_obj; 
        opt_ef(1) = opt_vol; 
        opt_df(0,Eigen::all) = opt_sns;
        opt_df(1,Eigen::all) = VectorXd::Ones(n_e)/n_e/1.;
        opt_cf(0,Eigen::all) = -2.0*opt_sns.cwiseProduct(opt_ex.cwiseInverse());
        opt_cf(1,Eigen::all) = VectorXd::Zero(n_e);
        VectorXd opt_sxl = VectorXd::Zero(n_e);
        VectorXd opt_sxu = VectorXd::Zero(n_e);
        opt_cs = VectorXi::Ones(1);
//      opt_la = VectorXd::Zero(1);
        opt_sxl = opt_xl.cwiseMax(opt_ex - 0.1*(opt_xu-opt_xl));
        opt_sxu = opt_xu.cwiseMin(opt_ex + 0.1*(opt_xu-opt_xl));
//
//  solve subproblem and update design variables
//
        VectorXd opt_nx = VectorXd::Zero(n_e);
        if(fdc){
            opt_ex(fdi) = opt_ex(fdi) + 1e-6;
        }else{
            err=dqpsub(opt_ex, opt_ef, opt_df, opt_cf, opt_sxl, opt_sxu, opt_cs, opt_la, opt_nx);
            opt_ex = opt_nx;
        }
//
//  write output
//
        t1 = high_resolution_clock::now();
        wrte_lvtk(argv[3],itr, nds, els, gss, mat, mat_aux, dfs_sol, dfs_pre, dfs_lds, opt_ex, opt_ro);
        t2 = high_resolution_clock::now();
        auto dur_out = duration_cast<milliseconds>(t2 - t1).count();
        auto dur_tot = duration_cast<milliseconds>(t2 - t0).count();
        if(fdc){

        }else{
            std::cout << "Output written (" << dur_out << ")\n";
            std::cout << "Total time (iteration) : " << dur_tot << "\n";
            std::cout << "===================================================" << "\n";
        }
//
//  possibly terminate
//

//
    }
//
//  write complete time
//
}
