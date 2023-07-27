//
#include <Eigen/Core>
#include <iostream>
#include <LBFGSB.h>  // Note the different header file
//
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;
using Eigen::all;
using Eigen::seq;
using Eigen::last;
using namespace LBFGSpp;
//
class Dual
{
private:
    VectorXd x_k;
    VectorXd g;
    MatrixXd dg;
    VectorXd dx_l;
    VectorXd dx_u;
    MatrixXd dc;
public:
    Dual(VectorXd& x_k_, VectorXd& g_, MatrixXd& dg_, VectorXd& dx_l_, VectorXd& dx_u_, MatrixXd& dc_) : x_k(x_k_), g(g_), dg(dg_), dx_l(dx_l_), dx_u(dx_u_), dc(dc_) {}
    double operator()(const VectorXd& x_d, VectorXd& grad)
    {
//
//      update primal variables
//
        VectorXd tmp1;
        VectorXd tmp2;
        VectorXd x;
        tmp1 = (dc(0,all) + x_d.transpose()*dc(seq(1,last),all));
        tmp2 = dg(0,all) + x_d.transpose()*dg(seq(1,last),all);
        x = ((x_k - tmp2.cwiseProduct(tmp1.cwiseInverse())).cwiseMax(dx_l)).cwiseMin(dx_u);
//
//      compute dual function value
//
        double W = -g(0);
        W = W - dg(0,all).dot(x-x_k) - 1./2.*tmp1.dot((x-x_k).cwiseProduct(x-x_k));
        W = W - x_d.transpose()*( g(seq(1,last)) + dg(seq(1,last),all)*(x-x_k) );
//
//      compute gradient of dual 
//
        grad = -g(seq(1,last));
        grad = grad - dg(seq(1,last),all)*(x-x_k);
        grad = grad - 1./2.*dc(seq(1,last),all)*((x-x_k).cwiseProduct(x-x_k));

        return W;
    }
};

int dqpsub(VectorXd& x_k, VectorXd& g, MatrixXd& dg, MatrixXd& dc, VectorXd& dx_l, VectorXd& dx_u, VectorXi& c_s, VectorXd& x_d, VectorXd& x)
{
    // set up parameters
    LBFGSBParam<double> param;
    param.epsilon = 1e-12;
    param.ftol = 1e-12;
    param.max_linesearch = 1000000;
    param.m = 20;
//
    // create solver and function object
    LBFGSBSolver<double> solver(param);  
    Dual fun(x_k, g, dg, dx_l, dx_u, dc);
//
    // dual variable bounds
    VectorXd ub = VectorXd::Ones(c_s.size())*1e8;
    VectorXd lb = -1e8*VectorXd::Ones(c_s.size());
    for(int i=0; i<c_s.size(); i++){
        if(c_s(i))
        {
            lb(i) = 0.0;
        }
    }
//
    double fx;
    int niter = solver.minimize(fun, x_d, fx, lb, ub);
//
//  std::cout << niter << " iterations" << std::endl;
//  std::cout << "x_d = \n" << x_d.transpose() << std::endl;
//  std::cout << "===================================================" << "\n";
//  std::cout << "f(x) = " << fx << std::endl;
//
//  update primal variables
    VectorXd tmp1; VectorXd tmp2;
    tmp1 = (dc(0,all) + x_d.transpose()*dc(seq(1,last),all));
    tmp2 = dg(0,all) + x_d.transpose()*dg(seq(1,last),all);
    x = ((x_k - tmp2.cwiseProduct(tmp1.cwiseInverse())).cwiseMax(dx_l)).cwiseMin(dx_u);
//
    return 0;
}

