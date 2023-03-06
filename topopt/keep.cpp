
int scg_saad(int n, SpMat A, VectorXd& b, VectorXd& x)
{
//
    double tmp;
    VectorXd r;
    VectorXd r1;
    VectorXd p;
    double alph = 0.;
    double beta = 0.;
    std::vector<T> coefs;
    coefs.reserve(n);
//
    SpMat I(n,n);
    for(int i =0; i<n; i++){
        coefs.push_back(T(i,i,1.));
    }
    I.setFromTriplets(coefs.begin(), coefs.end());
//
    r = b - A*x;
    p = r;
//
    double err = 1.;
//
    for(int k=0;k<1000000;k++){
//
//       tmp = b.dot(I*x); // we do it this way for sake of (efficient) generality (benefits from MT)
                         // b.dot(x) is faster in general.
//
        alph = r.dot(I*r)/(p.dot(A*p)); //!! bad way to do it
        x = x + alph*p;
        r1 = r - alph*A*p;
        beta = r1.dot(I*r1)/(r.dot(I*r)); //!! bad way to do it
        r=r1;
        p=r+beta*p;
//
        err = r.lpNorm<Eigen::Infinity>();
        if (err < 1e-6){
            std::cout << k << " " << err << endl;
            return 0;
        }
//
    }
//
  return 1;
}
//
