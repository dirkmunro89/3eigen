#include "main.h"
#include <chrono>
//
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;
//
int scal (int n, SpMat& A, VectorXd& b, VectorXd& d)
{
//
//  Copy of Hond; can be done nicer in Eigen frame
//
    VectorXd r;
    VectorXd tmp;
    SpMat mmp;
//
    tmp = A.diagonal().array().sqrt();
    for(int i = 0; i<n; i++)
    {
        d(i)=1.0/tmp(i);
    }
    b = b.cwiseProduct(-d);
//
    for (int k=0; k<A.outerSize(); ++k) // for each row
    {
        for (SpMat::InnerIterator it(A,k); it; ++it) // over all column entries
        {
            it.valueRef() = it.value()*d(it.row())*d(it.col());
        }
    }
//
    return 0;
//
}
int precon(int nnz, int n, ArrayXd& tmp, ArrayXi& ia, ArrayXi& iz, double alpa)
{
//
//  Copy of Hond; can be done nicer in Eigen frame
//
    double factor; 
    int i, j, k, l, m, jlo, jup, klo, kup, llo, lup, id;
//
    factor = 1.0/(1.0 + alpa);
    ArrayXd C = ArrayXd::Zero(tmp.size());
//
//  scale all off diagonal terms by factor
//
    C(0) = tmp(0);
    for(i=1; i<n; i++)
    {
        jlo=iz(i-1)+1;
        jup=iz(i);
        C(jup) = tmp(jup);
        for(j=jlo; j<jup; j++){
            C(j) = tmp(j)*factor;
        }
    }
//
//  can not see how to do this in terms of sparse mat vec operations
//  will be amazing if one could... people do do it on GPUs etc., so it is possible, i think
//
    for(i=1; i<n; i++)
    {
        jlo=iz(i-1)+1;
        jup=iz(i);
        for(j=jlo; j<jup; j++)
        {
            C(j) /= C(iz(ia(j))); 
            klo = j+1;
            kup = jup;
            for(k=klo; k<=kup; k++)
            {
                m = ia(k); 
                llo = iz(m-1)+1;
                lup = iz(m);
                for(l=llo; l<=lup; l++)
                {
                    if(ia(l)>ia(j)) break; 
                    if(ia(l)<ia(j)) continue; 
                    C(k) -= C(j)*C(l);
                    break;
                }
            }
        }
        id = iz(i);
        if (C(id)<1.0e-6){
            return 0;
        }
        C(id) = sqrt(C(id));
    }
//
    tmp=C;
    return 1;
//
}
int Mrhor(ArrayXd& C, int n, ArrayXi& ia, ArrayXi& iz, VectorXd& r, VectorXd& rho)
{
//
//  Copy of Hond; can be done nicer in Eigen frame
//
    int             i=0, j=0, jlo=0, jup=0;
    double          s=0.0;
//
//  rho = CT C r
//
    rho(0) = r(0);
    for (i=1; i<n; i++)
    {
        s = 0.0;
        jlo = iz(i-1)+1;            /*..first non-zero element in current row...... */
        jup = iz(i);                /*..diagonal element in current row............ */
        for (j=jlo; j<jup; j++)     /*..all non-zero off-diagonal element.......... */
            s = s + C(j)*rho(ia(j));
        rho[i] = (r(i)-s)/C(jup);
    }
    for (i=n-1; i>0; i--)
    {
        rho(i) = rho(i) / C(iz(i));
        jlo = iz(i-1)+1;            /*..first non-zero element in current row...... */
        jup = iz(i)-1;              /*..diagonal element in current row............ */
        for (j=jlo; j<=jup; j++)    /*..all non-zero off-diagonal element.......... */
            rho(ia(j)) = rho(ia(j)) - C(j)*rho(i);
    }
    return 0;
}
int pcg_hond(int n, SpMat& A, VectorXd& b, VectorXd& x)
{
//
    VectorXd r;
    VectorXd p;
    VectorXd z;
    VectorXd g;
    VectorXd rho;
    VectorXd tmp;
    double  gz, qk, ekm1, rr, rrho, rrho1;
    int iam=0;
    double qam=0;
    double ram=0;
    double err=0.;
    double c1=0.005;
    int i,j,k,c;
//
//  scaling with diag
//
    VectorXd b0=b; 
    VectorXd d(n);
    err=scal(n,A,b,d);
//
//  hond arrays
//
    int nnz = A.nonZeros();
    ArrayXi iz(n);
    ArrayXi ia((nnz-n)/2+n); // pointer to diagonal element
    ArrayXd C((nnz-n)/2+n);
//
//  we pack the upper triangular row wise into honds arrays,
//  that is equivalent to his lower triangular stacked column-wise
//
//
    c = 0;
    int cc = 0;
    for (int k=0; k<A.outerSize(); ++k) // for each row
    {
        for (SpMat::InnerIterator it(A,k); it; ++it) // over all column entries
        {
            if (it.row() > it.col()){
                ia(cc) = it.col();
                C(cc) = it.value();
                cc=cc+1;
            }    
            if (it.row() == it.col()){
                iz(k) = cc;
                ia(cc) = it.col();
                C(cc) = it.value();
                cc=cc+1;
            }
        }
        c=c+1;
    }
//
//  std::cout << "Check sym nnz: " << cc << " " << nnz/2 + n/2 << endl;
//
//  init sol and residual vectors
//
    x = VectorXd::Zero(n); // can of course do init start here
    r = A*x+b;
    rho = VectorXd::Zero(n);
    for(i = 0; i<n; i++){
        err = abs(r(i));
        if(err>1e-20){
            qam = qam+err;
            iam=iam+1;
        }
    }
    if(iam>0){
        qam=qam/iam;
    }
    qam=qam*1e-2;
//
    auto t1 = high_resolution_clock::now();
    c= 0;
    double alpa = 0.;
//  std::cout << "alpha = "<< alpa << endl;
    err=precon(nnz,n,C,ia,iz,alpa);
    while(err == 0){
        if(alpa<=0.){ 
            alpa = 0.005;
        }
        alpa=alpa+alpa;
//      std::cout << err << " alpha = "<< alpa << endl;
        err=precon(nnz,n,C,ia,iz,alpa);
        c=c+1;
        if(c==1000){
            return 1;
        }
    }
    auto t2 = high_resolution_clock::now();
    auto dur_chl = duration_cast<milliseconds>(t2 - t1).count();
//  std::cout << "Preconditioned (" << dur_chl << ")\n";
//
//
//  try making sparse rep...
//
/*    std::vector<T> coefs;
    coefs.reserve(nnz);///2+n/2);
    SpMat Cs(n,n);
//
    int jlo;
    int jup;
    coefs.push_back(T(0,0,C(iz(0))));
    for (i=1; i<n; i++){
        jlo = iz(i-1)+1;
        jup = iz(i);
        coefs.push_back(T(i,i,C(iz(i))));
        for (j=jlo; j<jup; j++){
            coefs.push_back(T(i,ia(j),C(j)));
        }
    }
    Cs.setFromTriplets(coefs.begin(), coefs.end());*/
//
    t1 = high_resolution_clock::now();
//
//  main loop
//
    for(k=1;k<1000000;k++){
//
//      M rho = r, M=C CT
//      rho = CT C r
//
//        tmp=r;
//        Cs.triangularView<Eigen::Lower>().solveInPlace(tmp);
//        Cs.triangularView<Eigen::Lower>().adjoint().solveInPlace(tmp);
//        rho=tmp;
        err=Mrhor(C,n,ia,iz,r,rho);
        rrho = r.dot(rho);
//
        rr = r.dot(r);
        if (k!=1 && ram<=c1*qam){
            break;
        }
        if(k!=1){
            ekm1=rrho/rrho1;
            g=ekm1*g-rho;
        }else{
            g = -rho;
        }
        z = A*g;
        gz = g.dot(z);
        qk = rrho/gz;
        ram = 0.;
        x = x+qk*g;
        r = r+qk*z;
        ram = r.lpNorm<Eigen::Infinity>();
        rrho1=rrho;
//
    }
//  std::cout << "iteration = " << k << " error = " << ram << " limit = "<< c1*qam << endl;
    t2 = high_resolution_clock::now();
    auto dur_itr = duration_cast<milliseconds>(t2 - t1).count();
//  std::cout << " ("<< dur_itr << ")" << endl;
//
//  back-scale solution vector and load vector
//
    x = x.cwiseProduct(d);
    b = b0;//.cwiseProduct(-d.pow(-1.));
//
    return 0;
//
}
int scg_hond(int n, SpMat& A, VectorXd& b, VectorXd& x)
{
//
    VectorXd r;
    VectorXd p;
    VectorXd z;
    VectorXd tmp;
    std::vector<T> coefs;
    coefs.reserve(n);
    double pz, qk, ekm1, rr, rro;
    int iam=0;
    double qam=0;
    double ram=0;
    double err=0.;
    double c1=0.005;
    int i,j,k;
//
//  scaling of diagonal
//
    VectorXd b0=b; 
    VectorXd d(n); 
    err=scal(n,A,b,d);
//
//  init sol, search, and residual vectors
//
    x = VectorXd::Zero(n);
    r = A*x+b;
    p = -r;
    for(i = 0; i<n; i++){
        err = abs(r(i));
        if(err>1e-20){
            qam = qam+err;
            iam=iam+1;
        }
    }
    if(iam>0){
        qam=qam/n;
    }
//  qam=1e-9; for finite differences, set this
//
//  main loop
//
    for(k=1;k<1000000;k++){
//
        rr = r.dot(r);
        if (k!=1 && ram<=c1*qam){
            break;
        }
        if(k!=1){
            ekm1=rr/rro;
            p=ekm1*p-r;
        }
        z = A*p;
        pz = p.dot(z);
        qk = rr/pz;
        ram = 0.;
        x = x+qk*p;
        r = r+qk*z;
        ram = r.lpNorm<Eigen::Infinity>();
        rro=rr;
//
    }
//
//  std::cout << "iteration = " << k << " error = " << ram << " limit = "<< c1*qam << endl;
//
//  back-scale solution vector
//
    x = x.cwiseProduct(d);
    b = b0;//.cwiseProduct(-d.pow(-1));
//
  return 0;
}
