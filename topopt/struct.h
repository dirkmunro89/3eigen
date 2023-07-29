//
// structures
//
struct { // for optimization data
    VectorXd _ex;
    VectorXd _exh;
    VectorXd _la;
    VectorXd _ef;
    VectorXd _efh;
    MatrixXd _df;
    MatrixXd _cf;
    VectorXd _xl;
    VectorXd _xu;
    VectorXd _dxl;
    VectorXd _dxu;
    VectorXd _mv;
    VectorXi _cs;
} opt;
struct { // for last acceptable subproblem data
    VectorXd _ex;
    VectorXd _ef;
    MatrixXd _df;
    MatrixXd _cf;
    VectorXd _dxl;
    VectorXd _dxu;
    VectorXd _la;
} acc;
struct { // for density filtering data
    VectorXd _ro;
    SpMat _H;
} flt;
