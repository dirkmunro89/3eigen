//
#include "main.h"
//
int DiffShapFunc(double xi, double eta, double zeta, VectorXd& dNdxi, VectorXd& dNdeta, VectorXd& dNdzeta);
int aux(MatrixXd& c, MatrixXd& B, ArrayXd& jac, int redInt);
//
int dmat(MatrixXd& c, MatrixXd& D, ArrayXd& jac)
{
//
//    D.setZero(6,24);
//
    int i;
    double gx, gy;
    double x1,y1,x2,y2,x3,y3,x4,y4;
    double jaci;
//
//    MatrixXd B =;
//
    int redInt = 0;
//
    int err = aux(c,D,jac,redInt);
//
 //   D = B;
//
return 0;
//
}
//
int aux(MatrixXd& c, MatrixXd& B, ArrayXd& jac, int redInt){
//
    VectorXd X(8); VectorXd Y(8); VectorXd Z(8);
//
    X = c(Eigen::all,0);
    Y = c(Eigen::all,1);
    Z = c(Eigen::all,2);
//
    double dx = 1;
    double dy = 1;
    double dz = 1;
//
//  X(0)= 0.; X(1)= dx; X(2)= dx; X(3)= 0.; X(4)= 0.; X(5)= dx; X(6)= dx; X(7)= 0.;
//  Y(0)= 0.; Y(1)= 0.; Y(2)= dy; Y(3)= dy; Y(4)= 0.; Y(5)= 0.; Y(6)= dy; Y(7)= dy;
//  Z(0)= 0.; Z(1)= 0.; Z(2)= 0.; Z(3)= 0.; Z(4)= dz; Z(5)= dz; Z(6)= dz; Z(7)= dz;
//
    ArrayXd GP(2);
    if(redInt){
        GP(0) = 0.;
    }else{
        GP(0) = -0.577350269189626;
        GP(1) = 0.577350269189626;
    }
//
    MatrixXd alpha1 = MatrixXd::Zero(6,3);
    MatrixXd alpha2 = MatrixXd::Zero(6,3);
    MatrixXd alpha3 = MatrixXd::Zero(6,3);
//
    alpha1(0,0) = 1.0; alpha1(3,1) = 1.0; alpha1(5,2) = 1.0;
    alpha2(1,1) = 1.0; alpha2(3,0) = 1.0; alpha2(4,2) = 1.0;
    alpha3(2,2) = 1.0; alpha3(4,1) = 1.0; alpha3(5,0) = 1.0;
//
    VectorXd  dNdxi(8);
    VectorXd  dNdeta(8);
    VectorXd  dNdzeta(8);
    MatrixXd  Jay(3,3);
    MatrixXd  invJ(3,3);
    MatrixXd  beta(6,3);
//    B = MatrixXd::Zero(6*8,24); // Note: Small enough to be allocated on stack
//    B.setZero(48,24);
    MatrixXd  J(3,3);
    ArrayXd dN;
//
    int count = 0;
    for(int ii = 0; ii < 2 - redInt; ii++)
    {
        for(int jj = 0; jj < 2 - redInt; jj++)
        {
            for(int kk = 0; kk < 2 - redInt; kk++)
            {
                double xi   = GP(ii); double eta  = GP(jj); double zeta = GP(kk);

                int err = DiffShapFunc(xi, eta, zeta, dNdxi, dNdeta, dNdzeta);

                J(0,0)  =  dNdxi.dot(X); J(0,1)  =  dNdxi.dot(Y); J(0,2)  =  dNdxi.dot(Z);
                J(1,0)  =  dNdeta.dot(X); J(1,1)  =  dNdeta.dot(Y); J(1,2)  =  dNdeta.dot(Z);
                J(2,0)  =  dNdzeta.dot(X); J(2,1)  =  dNdzeta.dot(Y); J(2,2)  =  dNdzeta.dot(Z);

                double detJ = J(0,0) * (J(1,1) * J(2,2) - J(2,1) * J(1,2)) -
                       J(0,1) * (J(1,0) * J(2,2) - J(2,0) * J(1,2)) +
                       J(0,2) * (J(1,0) * J(2,1) - J(2,0) * J(1,1));
                jac( count) = detJ;
                invJ(0,0) = (J(1,1) * J(2,2) - J(2,1) * J(1,2)) / detJ;
                invJ(0,1) = -(J(0,1) * J(2,2) - J(0,2) * J(2,1)) / detJ;
                invJ(0,2) = (J(0,1) * J(1,2) - J(0,2) * J(1,1)) / detJ;
                invJ(1,0) = -(J(1,0) * J(2,2) - J(1,2) * J(2,0)) / detJ;
                invJ(1,1) = (J(0,0) * J(2,2) - J(0,2) * J(2,0)) / detJ;
                invJ(1,2) = -(J(0,0) * J(1,2) - J(0,2) * J(1,0)) / detJ;
                invJ(2,0) = (J(1,0) * J(2,1) - J(1,1) * J(2,0)) / detJ;
                invJ(2,1) = -(J(0,0) * J(2,1) - J(0,1) * J(2,0)) / detJ;
                invJ(2,2) = (J(0,0) * J(1,1) - J(1,0) * J(0,1)) / detJ;


                for (int ll = 0; ll < 3; ll++) {
                    // Add contributions from the different derivatives
                    if (ll == 0) {
                        dN = dNdxi;
                    }
                    if (ll == 1) {
                        dN = dNdeta;
                    }
                    if (ll == 2) {
                        dN = dNdzeta;
                    }
                    beta.setZero(6,3);
                    // Assemble strain operator 
                    for (int i = 0; i < 6; i++) {
                        for (int j = 0; j < 3; j++) {
                        beta(i,j)=invJ(0,ll)*alpha1(i,j)+invJ(1,ll)*alpha2(i,j)+invJ(2,ll)*alpha3(i,j);
                        }
                    }
                    // Add contributions to strain-displacement matrix
                    for (int i = 0; i < 6; i++) {
                        for (int j = 0; j < 24; j++) {
                            B(i + 6*count,j) = B(i + 6*count,j) + beta(i,j % 3) * dN(j / 3);
                        }
                    }
                }
                count=count+1;
            }
        }
    }
//
return 0;
}
//
int DiffShapFunc(double xi, double eta, double zeta, VectorXd& dNdxi, VectorXd& dNdeta, VectorXd& dNdzeta){
//
    dNdxi(0) = -0.125 * (1.0 - eta) * (1.0 - zeta);
    dNdxi(1) = 0.125 * (1.0 - eta) * (1.0 - zeta);
    dNdxi(2) = 0.125 * (1.0 + eta) * (1.0 - zeta);
    dNdxi(3) = -0.125 * (1.0 + eta) * (1.0 - zeta);
    dNdxi(4) = -0.125 * (1.0 - eta) * (1.0 + zeta);
    dNdxi(5) = 0.125 * (1.0 - eta) * (1.0 + zeta);
    dNdxi(6) = 0.125 * (1.0 + eta) * (1.0 + zeta);
    dNdxi(7) = -0.125 * (1.0 + eta) * (1.0 + zeta);
    // With respect to eta:
    dNdeta(0) = -0.125 * (1.0 - xi) * (1.0 - zeta);
    dNdeta(1) = -0.125 * (1.0 + xi) * (1.0 - zeta);
    dNdeta(2) = 0.125 * (1.0 + xi) * (1.0 - zeta);
    dNdeta(3) = 0.125 * (1.0 - xi) * (1.0 - zeta);
    dNdeta(4) = -0.125 * (1.0 - xi) * (1.0 + zeta);
    dNdeta(5) = -0.125 * (1.0 + xi) * (1.0 + zeta);
    dNdeta(6) = 0.125 * (1.0 + xi) * (1.0 + zeta);
    dNdeta(7) = 0.125 * (1.0 - xi) * (1.0 + zeta);
    // With respect to zeta:
    dNdzeta(0) = -0.125 * (1.0 - xi) * (1.0 - eta);
    dNdzeta(1) = -0.125 * (1.0 + xi) * (1.0 - eta);
    dNdzeta(2) = -0.125 * (1.0 + xi) * (1.0 + eta);
    dNdzeta(3) = -0.125 * (1.0 - xi) * (1.0 + eta);
    dNdzeta(4) = 0.125 * (1.0 - xi) * (1.0 - eta);
    dNdzeta(5) = 0.125 * (1.0 + xi) * (1.0 - eta);
    dNdzeta(6) = 0.125 * (1.0 + xi) * (1.0 + eta);
    dNdzeta(7) = 0.125 * (1.0 - xi) * (1.0 + eta);
//
return 0;
//
}

