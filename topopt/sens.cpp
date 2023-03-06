#include "main.h"
//
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;
//
int sens(int n_d, MatrixXi& els, MatrixXd& nds, VectorXd& dfs_sol, VectorXd& opt_ro, VectorXd& opt_sns, double *opt_obj){
//
    int e, i, i1, i2, j, j1, j2, k, c, err, gss, etp;
    int gind1, gind2, lind1, lind2;
    ArrayXi con_e; MatrixXd nds_e; MatrixXd kay_e;
//
    auto t1 = high_resolution_clock::now();
    auto t11 = high_resolution_clock::now();
    auto t22 = high_resolution_clock::now();
    double dur_elm = 0; double dur_trp = 0;
//
// element matrices (just init here, each is computed individually; keep this way, for unstructured)
//
    MatrixXd D; //= MatrixXd::Zero(48,24); //!!
    MatrixXd S; //= MatrixXd::Zero(48, 48); //!!
    ArrayXd jac; //= ArrayXd::Zero(8); //!!
//
// loop over elements
//
    for(e = 0; e < els.rows(); e++){
//
        t11 = high_resolution_clock::now();
//
        D.setZero(48,24);
        S.setZero(48,48);
        jac.setZero(8);
//
//      coordinates of elements' nodes
        con_e = els(e,Eigen::all);
        nds_e = nds(con_e,Eigen::all);
//
//      get D matrix
        err=dmat(nds_e, D, jac);
//
//      get S matrix
        err=smat(S, jac);
//
//      make element k matrix 
        kay_e = D.transpose() * S * D ;
//
        t22 = high_resolution_clock::now();
        duration<double, std::milli> dur_tmp1 = t22 - t11;
        dur_elm = dur_elm + dur_tmp1.count();
//
//      Can/should definitely be done nicer with Eigen structures
//
        t11 = high_resolution_clock::now();
        for(i1 = 0; i1 < con_e.size(); i1++){
            for(j1 = 0; j1 < nds_e.cols(); j1++){
                gind1 = nds_e.cols()*con_e[i1] + j1;
                lind1 = i1*nds_e.cols() + j1;
                for(i2 = 0; i2 < con_e.size(); i2++){
                    for(j2 = 0; j2 < nds_e.cols(); j2++){
                        gind2 = nds_e.cols()*con_e[i2] + j2;
                        lind2 = i2*nds_e.cols() + j2;
                        *opt_obj= *opt_obj +  pow(opt_ro(e),3.)*kay_e(lind1,lind2)*dfs_sol(gind1)*dfs_sol(gind2);
                        opt_sns(e) = opt_sns(e) - 3.*pow(opt_ro(e),2.)*kay_e(lind1,lind2)*dfs_sol(gind1)*dfs_sol(gind2);
                    }
                }
            }
        }
//
        t22 = high_resolution_clock::now();
        duration<double, std::milli> dur_tmp2 = t22 - t11;
        dur_trp = dur_trp + dur_tmp2.count();
//
    }
//
    auto t2 = high_resolution_clock::now();
    auto dur_asy = duration_cast<milliseconds>(t2 - t1).count();
    std::cout << "Compliance sensitivity (" << dur_asy << ")\n";
    std::cout << "- Element level matrix operations : " << dur_elm << "\n";
    std::cout << "- Local-global dof mult. : " << dur_trp << "\n";
//
return 0;
}

