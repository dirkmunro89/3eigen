//
#include "main.h"
//
//typedef Eigen::Matrix<double, 12, 12> MatrixS;
//
int smat(MatrixXd& S, ArrayXd& jac)
{
//
//  S = MatrixXd::Zero(6*8,6*8);
//
    double E = 1.;
    double nu = 0.3;
//
    double lambda = E*nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double mu     = E*1.0 / (2.0 * (1.0 + nu));
//
    for(int i = 0; i < 8; i++){
        lambda = E*nu / ((1.0 + nu) * (1.0 - 2.0 * nu)) * jac(i);
        mu     = E*1.0 / (2.0 * (1.0 + nu)) * jac(i);
        S(0+6*i,0+6*i) = lambda + 2.0 * mu; 
        S(0+6*i,1+6*i) = lambda; 
        S(0+6*i,2+6*i) = lambda;
        S(1+6*i,0+6*i) = lambda; 
        S(1+6*i,1+6*i) = lambda + 2.0 * mu; 
        S(1+6*i,2+6*i) = lambda;
        S(2+6*i,0+6*i) = lambda; 
        S(2+6*i,1+6*i) = lambda; 
        S(2+6*i,2+6*i) = lambda + 2.0 * mu;
        S(3+6*i,3+6*i) = mu;
        S(4+6*i,4+6*i) = mu;
        S(5+6*i,5+6*i) = mu;
    }
//
return 0;
//
}
//
