#include "main.h"
#include <set>
//
int filt_init(int n_e, MatrixXi& els, MatrixXd& nds, SpMat& H, int lvl)
{
//
    std::vector<T> coefs;
//  coefs.reserve(200*n_e); // 27
    ArrayXi con_e;
    MatrixXi noc_e;
//
//  int lvl;
//
    cout << lvl << endl;
//
//  lvl=2;
//
// make inverse of connect matrix, matrix of centroid coordinates, and get length scale
//
    MatrixXi noc = -1*MatrixXi::Ones(nds.rows(), 9);
    ArrayXi cnt = VectorXi::Zero(nds.rows());
    MatrixXd cen = MatrixXd::Zero(n_e,3);
    VectorXd cmp;
    double len = 1e8;
    VectorXd cmp_0;
    VectorXd Hsum = VectorXd::Zero(n_e);
//
    for(int e = 0; e < n_e; e++){
        con_e = els(e,Eigen::all); 
        cmp = VectorXd::Zero(3);
        for(int i = 0; i < con_e.size(); i++){
            noc(con_e(i),cnt(con_e(i))) = e;
            cmp = cmp + nds(con_e(i),Eigen::all).transpose()/con_e.size();
        }
        cen(e,Eigen::all) = cmp;
        if(e==0){
            cmp_0 = cmp;
        }else{
            len = min(len, (cmp-cmp_0).norm());
        }
        cnt(con_e) = cnt(con_e) + 1;
        if ((cnt(con_e) == 9).any()){
            std::cout << "ERROR" << endl;
        }
    }
//
// for each element
//
    for(int e = 0; e < n_e; e++){
//
//     get nodes connected to element e
        con_e = els(e,Eigen::all); 
//     elements connected to each node (immediate neighbours)
        noc_e = noc(con_e,Eigen::all); 
//     make set of unique entries
        std::set<int> q0{noc_e.data(), noc_e.data() + noc_e.size()  }; 
//     erase -1 place-holder
        q0.erase(-1); 
//
//     for each level
        for(int l = 0; l < lvl; l++){ 
//       make an empty set for entries in this level
            std::set<int> qn;
//       loop over entries in q0, add entries to qn
            for(auto it = q0.begin(); it!=q0.end();it++) 
            {
                    con_e = els(*it,Eigen::all); 
                    noc_e = noc(con_e,Eigen::all); 
                    std::set<int> tmp{noc_e.data(), noc_e.data() + noc_e.size()  }; 
                    qn.insert(tmp.begin(),tmp.end());
            }
//        erase -1 placeholder
            qn.erase(-1); 
//        add to q0
            q0.insert(qn.begin(),qn.end());
        }
//
//     make filter coefs
        double krn;
        for(auto it = q0.begin(); it!=q0.end();it++)
        {
            krn = lvl*len*1.01 - (cen(e,Eigen::all)-cen(*it,Eigen::all)).norm();
            if(krn > 0.){
                coefs.push_back(T(e, *it, krn));
                Hsum(e) = Hsum(e) + krn;
            }
        }
//
    }
//
    H.setFromTriplets(coefs.begin(), coefs.end());
    H = H * Hsum.asDiagonal().inverse();
    cout << H.nonZeros() << endl;
//  H.conservativeResize(n_e,n_e);
//  cout << H.nonZeros() << endl;
//
return 0;
}
