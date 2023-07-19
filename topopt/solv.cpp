#include "main.h"
#include <chrono>
//
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;
//
int solv(int sol, int n_f, SpMat ffs_K, VectorXd& ffs_lds, VectorXd& ffs_sol){
//
    int err;
//
    auto t1 = high_resolution_clock::now();
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
        solver.setTolerance(1e-12);
        solver.compute(ffs_K);
        ffs_sol = solver.solve(ffs_lds);
//      std::cout << "- #iterations:     " << solver.iterations() << std::endl;
//      std::cout << "- estimated error: " << solver.error()      << std::endl;
    }
//
    auto t2 = high_resolution_clock::now();
    auto dur_sol = duration_cast<milliseconds>(t2 - t1).count();
//  std::cout << "Solved (" << dur_sol << ")\n";
//

//
return 0;
}
