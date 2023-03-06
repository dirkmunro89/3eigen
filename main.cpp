#include "main.h"
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
    ArrayXi con_e; MatrixXd nds_e; MatrixXd kay_e;
    VectorXd dfs_lds; VectorXi dfs_pre; VectorXd dfs_sol; 
    VectorXd ffs_lds;                   VectorXd ffs_sol; 
    VectorXi map_num; 
    std::vector<T> coefs;
//
//
    auto t0 = high_resolution_clock::now();
//
    if (argc > 1) {
        std::cout << "Input: " << argv[1] << " " << argv[2] << " " << argv[3] << endl;
    }else{
        std::cout << "Usage: " << argv[0] << " filename" << std::endl;
        return 1;
    }
    int thr;
    std::stringstream str1(argv[1]);
    str1 >> thr;
    int sol;
    std::stringstream str2(argv[2]);
    str2 >> sol;
//
    Eigen::setNbThreads(thr);
//
// read input deck (json currently)
//
    auto t1 = high_resolution_clock::now();
    read_json(argv[3], nds, els, &gss, mat, mat_aux, bds, lds);
    auto t2 = high_resolution_clock::now(); 
    auto dur_inp = duration_cast<milliseconds>(t2 - t1).count();
    std::cout << "Input data parsed (" << dur_inp << ")\n";
//
// reserve estimate of number of non-zeros
//  ... and define some other things ...
//
    nnz = pow(els.cols()*nds.cols(),2)*els.rows();
    coefs.reserve(nnz);
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
// DOF vectors containing 
//  - loads
//  - prescribed dofs, and 
//  - (eventually) the solution vector
//
//    SpMat ffs_K(n_d,n_d);
    SpMat ffs_K(n_f,n_f);
    dfs_lds.setZero(n_d); dfs_pre.setZero(n_d); dfs_sol.setZero(n_d);
    dfs_lds(lds.col(0)) = lds.col(1);
    dfs_pre(bds.col(0)) = Eigen::VectorXi::Ones(n_p); 
//
// Free DOF vectors for loads and solution vector
    ffs_lds.setZero(n_f); ffs_sol.setZero(n_f);
//
// mapping dof numbering vector with prescribed dofs removed
// and the corresponding load vector
//
    c=0;
    map_num.setZero(n_d);
    for(i  = 0; i < n_d ; i++){
        if(dfs_pre(i)){
            map_num(i) = -1;
        }else{
            map_num(i) = c;
            ffs_lds(map_num(i)) = dfs_lds(i);
            c=c+1;
        }
    }
//
    t1 = high_resolution_clock::now();
    auto t11 = high_resolution_clock::now();
    auto t22 = high_resolution_clock::now(); 
    double dur_elm = 0; double dur_trp = 0; 
//
    nnz_c = 0;
    nnz_s = 0;
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
        kay_e = D.transpose() * S * D;
//
        t22 = high_resolution_clock::now();
        duration<double, std::milli> dur_tmp1 = t22 - t11;
        dur_elm = dur_elm + dur_tmp1.count();
//
//      add little k to global k triplets (should be able to do nicer)
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
                        if( dfs_pre(gind2)){

                        }else if( dfs_pre(gind1) ){

                        }else{
                            coefs.push_back(T(map_num(gind1),map_num(gind2),kay_e(lind1,lind2)));
                        }
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
    t2 = high_resolution_clock::now();
    auto dur_asy = duration_cast<milliseconds>(t2 - t1).count();
    std::cout << "K triplets made (" << dur_asy << ")\n";
    std::cout << "- Element level matrix operations : " << dur_elm << "\n";
    std::cout << "- Making triplets : " << dur_trp << "\n";
//
    t1 = high_resolution_clock::now();
    ffs_K.setFromTriplets(coefs.begin(), coefs.end());
    std::cout << "Estimated nnz : " << nnz << " / Actual nnz : " << nnz_c << endl;
    nnz_c = ffs_K.nonZeros();
    t2 = high_resolution_clock::now();
    auto dur_spr = duration_cast<milliseconds>(t2 - t1).count();
    std::cout << "K made from triplets (" << dur_spr << ")\n";
    std::cout << "Estimated nnz : " << nnz << " / Actual nnz : " << nnz_c << endl;
    std::cout << "Filled : " << (double) nnz_c/n_f/n_f << endl;
//
    t1 = high_resolution_clock::now();
    if (sol==0){
        Eigen::SimplicialLDLT<SpMat> solver(ffs_K); 
        ffs_sol = solver.solve(ffs_lds); 
    }else if (sol==1){
        err=scg_hond(n_f, ffs_K, ffs_lds, ffs_sol);
    }else if (sol==2){
        err=pcg_hond(n_f, ffs_K, ffs_lds, ffs_sol);
    }else{
        Eigen::ConjugateGradient<SpMat,
            Eigen::Lower|Eigen::Upper,Eigen::DiagonalPreconditioner<double>> solver;
        solver.setMaxIterations(1000000);
        solver.setTolerance(1e-6);
        solver.compute(ffs_K);
        ffs_sol = solver.solve(ffs_lds);
        std::cout << "- #iterations:     " << solver.iterations() << std::endl;
        std::cout << "- estimated error: " << solver.error()      << std::endl;
    }
//
    t2 = high_resolution_clock::now();
    auto dur_sol = duration_cast<milliseconds>(t2 - t1).count();
    std::cout << "Solved (" << dur_sol << ")\n";
//
//  pack solution vector into complte dof vector
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
  t1 = high_resolution_clock::now();
  wrte_lvtk(argv[3], nds, els, gss, mat, mat_aux, dfs_sol, dfs_pre, dfs_lds);
  t2 = high_resolution_clock::now();
  auto dur_out = duration_cast<milliseconds>(t2 - t1).count();
  auto dur_tot = duration_cast<milliseconds>(t2 - t0).count();
  std::cout << "Output written (" << dur_out << ")\n";
  std::cout << "Total time : " << dur_tot << "\n";
//
}
