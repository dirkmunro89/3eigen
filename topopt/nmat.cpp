//
#include "main.h"
//
int ShapFunc(double xi, double eta, double zeta, VectorXd& N);
int aux2(MatrixXd& c, MatrixXd& B, ArrayXd& jac, int redInt);
//
int nmat(MatrixXd& c, MatrixXd& bN, ArrayXd& jac)
{
//
    int i;
    double gx, gy;
    double x1,y1,x2,y2,x3,y3,x4,y4;
    double jaci;
//
    int redInt = 0;
//
    int err = aux2(c,bN,jac,redInt);
//
return 0;
//
}
//
int aux2(MatrixXd& c, MatrixXd& bN, ArrayXd& jac, int redInt){
//
    VectorXd X(8); VectorXd Y(8); VectorXd Z(8);
//
    X = c(Eigen::all,0);
    Y = c(Eigen::all,1);
    Z = c(Eigen::all,2);
//
    ArrayXd GP(2);
    if(redInt){
        GP(0) = 0.;
    }else{
        GP(0) = -0.577350269189626;
        GP(1) = 0.577350269189626;
    }
//
    VectorXd  dNdxi(8);
    VectorXd  dNdeta(8);
    VectorXd  dNdzeta(8);
    MatrixXd  J(3,3);
    VectorXd N(8);
//
    int count = 0;
    for(int ii = 0; ii < 2 - redInt; ii++)
    {
        for(int jj = 0; jj < 2 - redInt; jj++)
        {
            for(int kk = 0; kk < 2 - redInt; kk++)
            {
                double xi   = GP(ii); double eta  = GP(jj); double zeta = GP(kk);

                int err = ShapFunc(xi, eta, zeta, N);
                err = DiffShapFunc(xi, eta, zeta, dNdxi, dNdeta, dNdzeta);

                J(0,0)  =  dNdxi.dot(X); J(0,1)  =  dNdxi.dot(Y); J(0,2)  =  dNdxi.dot(Z);
                J(1,0)  =  dNdeta.dot(X); J(1,1)  =  dNdeta.dot(Y); J(1,2)  =  dNdeta.dot(Z);
                J(2,0)  =  dNdzeta.dot(X); J(2,1)  =  dNdzeta.dot(Y); J(2,2)  =  dNdzeta.dot(Z);

                double detJ = J(0,0) * (J(1,1) * J(2,2) - J(2,1) * J(1,2)) -
                       J(0,1) * (J(1,0) * J(2,2) - J(2,0) * J(1,2)) +
                       J(0,2) * (J(1,0) * J(2,1) - J(2,0) * J(1,1));
                jac( count) = detJ;


                // Add contributions to strain-displacement matrix
                for (int j = 0; j < 8; j++) {
                    bN(0 + 3*count,3*j) = bN(0 + 3*count,3*j) +  N(j) * jac(count);
                    bN(1 + 3*count,3*j+1) = bN(1 + 3*count,3*j+1) + N(j) * jac(count);
                    bN(2 + 3*count,3*j+2) = bN(2 + 3*count,3*j+2) + N(j) * jac(count);
                }
                count=count+1;
            }
        }
    }
//
return 0;
}
//
int ShapFunc(double xi, double eta, double zeta, VectorXd& N){
//
    N(0) = 0.125 * (1.0 - xi) * (1.0 - zeta) * (1.0 - eta);
    N(1) = 0.125 * (1.0 + xi) * (1.0 - zeta) * (1.0 - eta);
    N(2) = 0.125 * (1.0 + xi) * (1.0 - zeta) * (1.0 + eta);
    N(3) = 0.125 * (1.0 - xi) * (1.0 - zeta) * (1.0 + eta);
    N(4) = 0.125 * (1.0 - xi) * (1.0 + zeta) * (1.0 - eta);
    N(5) = 0.125 * (1.0 + xi) * (1.0 + zeta) * (1.0 - eta);
    N(6) = 0.125 * (1.0 + xi) * (1.0 + zeta) * (1.0 + eta);
    N(7) = 0.125 * (1.0 - xi) * (1.0 + zeta) * (1.0 + eta);
//
return 0;
//
}

