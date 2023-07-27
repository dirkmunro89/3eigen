#include "main.h"
#include <set>
//
int filt_init(int n_e, MatrixXi& els, MatrixXd& nds, SpMat& H)
{
//
    std::vector<T> coefs;
    ArrayXi con_e;
    MatrixXi noc_e;
//
    int lvl;
//
    lvl=1;
//
//  make inverse of connect matrix
//
    MatrixXi noc = -1*MatrixXi::Ones(nds.rows(), 9);
    ArrayXi cnt = VectorXi::Zero(nds.rows());
    MatrixXi one; 
//
    coefs.reserve(27*n_e);
//
    for(int e = 0; e < n_e; e++){
        con_e = els(e,Eigen::all); 
        for(int i = 0; i < con_e.size(); i++){
            noc(con_e(i),cnt(con_e(i))) = e;
        }
        cnt(con_e) = cnt(con_e) + 1;
        if ((cnt(con_e) == 9).any()){
            std::cout << "ERROR" << endl;
        }
    }
//
    for(int e = 0; e < n_e; e++){
//
        con_e = els(e,Eigen::all); // nodes of element
        for(int l = 0; l < lvl; l++){ //!!!! would need mod for larger radii; ideally recursive
                                      // function
            noc_e = noc(con_e,Eigen::all); // elements connected to each node
            std::set<int> q{noc_e.data(), noc_e.data() + noc_e.size()  }; // get unique element ids
            for (auto it = q.begin(); it!=q.end();it++) // loop over neighbours
            {
                if (*it != -1){
                    coefs.push_back(T(e, *it, 1./(27.))); //removing -1 and adding factor 
//                  coefs.push_back(T(e, *it, 1./(q.size()))); //removing -1 and adding factor 
                                                               //2 for current element cancels here
                }
            }
            coefs.push_back(T(e, e, 1./27.));
        }
    }
//
    H.setFromTriplets(coefs.begin(), coefs.end());
//
return 0;
}
