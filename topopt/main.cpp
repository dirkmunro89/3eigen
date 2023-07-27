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
    VectorXd dfs_lds;
    VectorXi dfs_pre; 
    VectorXd dfs_sol; 
    VectorXd ffs_lds;                   
    VectorXd ffs_sol; 
    VectorXi map_num; 
//
//  for opt
//
    VectorXd opt_ex;
    VectorXd opt_exh;
    VectorXd opt_la;
    VectorXd opt_ro;
    VectorXd opt_ef;
    VectorXd opt_efh;
    MatrixXd opt_df;
    MatrixXd opt_cf;
    VectorXd opt_xl;
    VectorXd opt_xu;
    VectorXd opt_mv;
    VectorXi opt_cs;
//
//  for fd
//
    VectorXd fdc_ef0;
    fdc_ef0 = VectorXd::Zero(2);
//
    auto t0 = high_resolution_clock::now();
//
    if (argc > 5) {
        std::cout<<"Input: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<" "<<argv[6]<<endl;
    }else{
        std::cout<<"Input: "<<"thr"<<" "<<"sol"<<" "<<"input.json"<<" "<<"fdc"<<"vol"<<"nca"<<endl;
        return 1;
    }
    int thr; std::stringstream str1(argv[1]); str1 >> thr;
    int sol; std::stringstream str2(argv[2]); str2 >> sol;
    int fdc; std::stringstream str3(argv[4]); str3 >> fdc;
    double vol; std::stringstream str4(argv[5]); str4 >> vol;
    double opt_vf=vol;
    int nca; std::stringstream str5(argv[6]); str5 >> nca;
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
    n_e = els.rows(); // number of elements
    n_n = nds.rows(); // number of nodes
    n_d = nds.cols()*nds.rows(); // number of dofs
    n_p = bds.rows(); // number of prescribed dofs 
    n_f = n_d-n_p; // free dofs
    std::cout << "- Elements : " << n_e << "\n";
    std::cout << "- Nodes : " << n_n << "\n";
    std::cout << "- Dofs : " << n_d << "\n";
    std::cout << "- Free : " << n_f << "\n";
    if(nca){
    	std::cout << "Non-convex approximation selected" << "\n";
    }
//
//  for opt
//
    opt_ef = VectorXd::Zero(2);
    opt_efh = 1e8*VectorXd::Ones(2);
    opt_df = MatrixXd::Zero(2,n_e);
    opt_cf = MatrixXd::Zero(2, n_e);
    opt_xl = VectorXd::Zero(n_e);
    opt_xu = VectorXd::Ones(n_e);
    opt_mv = VectorXd::Ones(n_e)*0.1;
    opt_cs = VectorXi::Ones(1);
    opt_la = VectorXd::Zero(1);
    VectorXd opt_sxl = VectorXd::Zero(n_e);
    VectorXd opt_sxu = VectorXd::Zero(n_e);
    opt_cs = VectorXi::Ones(1);
//
//  for subproblem storage
//
    VectorXd acc_ex = opt_ex;
    VectorXd acc_ef = opt_efh;
    MatrixXd acc_df = opt_df;
    MatrixXd acc_cf = opt_cf;
    VectorXd acc_sxl = opt_sxl;
    VectorXd acc_sxu = opt_sxu;
    VectorXd acc_la = opt_la;
    double a = 1.0;
//
//  filter initialisation
//
    opt_ro = VectorXd::Zero(n_e);
    SpMat opt_H(n_e,n_e);
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
    dfs_lds(lds.col(0)) = lds.col(1);
    dfs_pre(bds.col(0)) = Eigen::VectorXi::Ones(n_p); 
//
// Free DOF vectors for loads and solution vector
    ffs_lds.setZero(n_f); 
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
    double opt_def = 1.0;
    double opt_scl;
    double opt_obj;
    VectorXd opt_tmp; 
    VectorXd opt_sns1; 
    VectorXd opt_sns2; 
//
//  if finite differences
//
    double fdx = 0;
    int fdi = (int) (n_e / 2);
    if(fdc){
        opt_ex = 0.5*VectorXd::Ones(n_e);
    }
//
    if(fdc==0){
        printf("======================================================================\n");
        printf("  k         f0         f1       dx       df       bw       tk       tt\n");
    }else{
        printf("=======================================\n");
    }
//
    for(int itr=0; itr<400; itr++){
//
        auto ti = high_resolution_clock::now();
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
            }else{
                dfs_sol(i) = ffs_sol(c);
                c=c+1;
            }
        }
//
//  compute compliance objective and sensitivities (similar to assembly)
//
        opt_obj = 0.;
        opt_tmp = VectorXd::Zero(n_e);
        opt_sns1 = VectorXd::Zero(n_e);
        err=sens(n_d,els,nds,dfs_sol,dfs_pre,opt_ro,opt_tmp,&opt_obj);
//
//  filter the sensitivities (need to check if transpose multiplication is fast...)
//
        opt_sns1 = opt_H.transpose()*opt_tmp;
        opt_tmp = -VectorXd::Ones(n_e)/n_e/opt_vf;
        opt_sns2 = opt_H.transpose()*opt_tmp;
//
//  scaling objective and its sensitivities
//
        if(itr==0){
            opt_scl = opt_obj;
        }
        opt_obj=opt_obj/opt_scl;
        opt_sns1=opt_sns1/opt_scl;
//
        opt_ef(0) = opt_obj;
        opt_ef(1) = -opt_ro.sum()/n_e/opt_vf + 1.;
//
//  calculate delta f
//
        opt_def = opt_ef(0) - acc_ef(0);
//
//  check if feasible descent step
//
        if(opt_def > 0. && fdc==0){
//
//      retrieve the last acceptable iterate and objective value / do not update subproblem
//
            opt_ex = acc_ex;
            opt_ef = acc_ef; // just for output
            a=2.0*a;
//
        }else{
//
//      construct subproblem
//
//
            if(fdc){
                if(itr == 0){
                    fdc_ef0=opt_ef;
                    opt_df(0,Eigen::all) = opt_sns1;
                    opt_df(1,Eigen::all) = opt_sns2;
                }else{
                    VectorXd opt_err = VectorXd::Zero(2);
                    printf("%8.1e ", fdx);
                    for(int j = 0; j<2; j++){
                        opt_err(j) = ((opt_ef(j)-fdc_ef0(j))/fdx - opt_df(j,fdi))/opt_df(j,fdi);
                        printf("%14.7e ", opt_err(j));
                    }
                    printf("\n");
                }
            }else{
//
//              spherical curvature
//
                double tmp; double nrm;
                opt_df(0,Eigen::all) = opt_sns1;
                opt_df(1,Eigen::all) = opt_sns2;
                opt_cf(0,Eigen::all) = 1e-6*VectorXd::Ones(n_e);
                if(itr>0){
                    nrm = (opt_ex-opt_exh).squaredNorm();
		    if(nca){
                    	tmp=2.*fabs(opt_efh(0)-opt_ef(0)-opt_df(0,Eigen::all).dot(opt_exh-opt_ex))/nrm;
		    }else{
                    	tmp=2.*(opt_efh(0)-opt_ef(0)-opt_df(0,Eigen::all).dot(opt_exh-opt_ex))/nrm;
	            }
                    opt_cf(0,Eigen::all) = (tmp*VectorXd::Ones(n_e)).cwiseMax(1e-6);
//                  opt_cf(0,Eigen::all) = -2.0*opt_sns.cwiseProduct(opt_ex.cwiseInverse());
                }
                opt_cf(1,Eigen::all) = VectorXd::Zero(n_e);
//
//              subproblem bounds
//
                opt_sxl = opt_xl.cwiseMax(opt_ex - 0.1*(opt_xu-opt_xl));
                opt_sxu = opt_xu.cwiseMin(opt_ex + 0.1*(opt_xu-opt_xl));
//
//              save this acceptable iterate
//
                acc_ex = opt_ex; acc_df = opt_df; acc_cf = opt_cf;
                acc_sxl = opt_sxl; acc_sxu = opt_sxu; acc_la = opt_la;
                opt_exh=opt_ex; acc_ef = opt_ef; a=1.0;
            }
        }
//
//  solve subproblem and update design variables
//
        double opt_dx = 0.;
        VectorXd opt_nx = VectorXd::Zero(n_e);
        if(fdc){
            fdx = 1e-16*pow(10.0,itr);
            opt_ex = 0.5*VectorXd::Ones(n_e);
            opt_ex(fdi) = opt_ex(fdi) + fdx;
        }else{
            MatrixXd tmp_cf = acc_cf*a;
            err=dqpsub(acc_ex, acc_ef, acc_df, tmp_cf, acc_sxl, acc_sxu, opt_cs, acc_la, opt_nx);
            opt_dx = (opt_nx - opt_ex).array().abs().maxCoeff();
            opt_ex = opt_nx; opt_efh = opt_ef;
        }
//
//  write output
//
        t1 = high_resolution_clock::now();
        wrte_lvtk(argv[3],itr,nds,els,gss,mat,mat_aux,dfs_sol,dfs_pre,dfs_lds,opt_ex,opt_ro);
        t2 = high_resolution_clock::now();
        auto dur_out = duration_cast<milliseconds>(t2 - t1).count();
        auto dur_itr = duration_cast<milliseconds>(t2 - ti).count();
        auto dur_tot = duration_cast<milliseconds>(t2 - t0).count();
        if(fdc){
            if(itr==18){
                break;
            }
        }else{
            double opt_baw = 0.0;
            for(int j = 0; j<n_e;j++){
                if( opt_ex(j) < opt_xl(j) + 1e-3 ){
                    opt_baw = opt_baw + 1.0/n_e;
                }else if ( opt_ex(j) > opt_xu(j) - 1e-3) {
                    opt_baw = opt_baw + 1.0/n_e;
                }
            }
            printf("%3d %10.3e %10.3e %8.1e %8.1e %8.1e %8.1e %8.1e\n", itr, opt_ef(0)*opt_scl, 
                opt_ef(1), opt_dx, opt_def, opt_baw, (double)dur_itr/1000, (double)dur_tot/1000.);
        }
        fflush(stdout);
//
//  possibly terminate
//
        if(fdc==0 && opt_dx < 1e-6 && a < 1.5 && itr > 100 ){
            break;
        }
//
    }
//
}
