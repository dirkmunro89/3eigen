////
#include "main.h"
#include "struct.h"
////
//
//main
//
////
int main(int argc, char *argv[])
{
//
//time
//
  auto t0 = high_resolution_clock::now();
  auto t1 = high_resolution_clock::now();
  auto t2 = high_resolution_clock::now();
  auto t21 = duration_cast<milliseconds>(t2-t1).count();
//
//check input
//
  printf("======================================================================\n");
  if (argc == 9) {
    printf("Input\n");
    printf("--------------------------------------------------------------------\n");
    cout<<" 1 - number of threads: "<<argv[1]<<endl<<
    " 2 - solver: "<<argv[2]<<endl<<
    " 3 - mesh file: "<<argv[3]<<endl<<
    " 4 - finite diff.: "<<argv[4]<<endl<<
    " 5 - volume: "<<argv[5]<<endl<<
    " 6 - non-convex approx.:  "<<argv[6]<<endl<<
    " 7 - filter level: "<<argv[7]<<endl<<
    " 8 - strict convex. param.:  "<<argv[8]<<endl;
  }else{
    cout<<"Input: "<<endl<<
    " 1 - number of threads (int) "<<endl<<
    " 2 - solver (0,1,2)"<<endl<<
    " 3 - mesh file (.json) "<<endl<<
    " 4 - finite diff. flag (0,1)"<<endl<<
    " 5 - volume (float)"<<endl<<
    " 6 - non-convex approx. flag (0,1)"<<endl<<
    " 7 - filter level (0,1,2,3..)"<<endl<<
    " 8 - strict convexity param. (float) "<<endl;
    return 1;
  }
//
//get input
//
  int thr; std::stringstream str1(argv[1]); str1 >> thr;
  int sol; std::stringstream str2(argv[2]); str2 >> sol;
  int fdc; std::stringstream str3(argv[4]); str3 >> fdc;
  double opt_vf; std::stringstream str4(argv[5]); str4 >> opt_vf;
  int nca; std::stringstream str5(argv[6]); str5 >> nca;
  flt._lvl; std::stringstream str6(argv[7]); str6 >> flt._lvl;
  double eps; std::stringstream str7(argv[8]); str7 >> eps;
//
//set threads
//
  Eigen::setNbThreads(thr);
//
//generic initializations
//
  int e, i, i1, i2, j, j1, j2, k, c, err, gss;
  int nnz, n_d, n_e, n_n, n_p, n_f;
//
  MatrixXd nds; MatrixXi els; 
  MatrixXd mat; string mat_aux;
  MatrixXd bds; MatrixXd lds;
  VectorXd dfs_lds;
  VectorXi dfs_pre; 
  VectorXd dfs_sol; 
  VectorXd ffs_lds;                   
  VectorXd ffs_sol; 
  VectorXi map_num; 
//
//for fd
//
  VectorXd fdc_ef0;
  fdc_ef0 = VectorXd::Zero(2);
//
//read input deck (json currently)
//
  printf("======================================================================\n");
  printf("Mesh\n"); 
  printf("--------------------------------------------------------------------\n");
//
  printf(" file: %s\n",argv[3]);
//
  t1 = high_resolution_clock::now();
  read_json(argv[3], nds, els, &gss, mat, mat_aux, bds, lds, opt._ex);
  t2 = high_resolution_clock::now(); t21 = duration_cast<milliseconds>(t2-t1).count();
//
  printf(" read in %7.1e seconds.\n",(double)t21/1000.);
  nnz = pow(els.cols()*nds.cols(),2)*els.rows();
  n_e = els.rows(); 
  n_n = nds.rows(); 
  n_d = nds.cols()*nds.rows();
  n_p = bds.rows(); 
  n_f = n_d-n_p;
  printf(" elements: %d\n",n_e);
  printf(" nodes: %d\n",n_n);
  printf(" dofs: %d\n",n_d);
  printf(" prescribed: %d\n",n_p);
  printf(" free: %d\n",n_f);
  printf(" reserved nnz: %d\n",nnz);
//
//filter init.
//
  printf("======================================================================\n");
  printf("Filtering\n");
  printf("--------------------------------------------------------------------\n");
//
  printf(" levels: %d\n",flt._lvl);
//
  flt._ro = VectorXd::Zero(n_e);
  flt._H = SpMat(n_e,n_e);
  flt._nze= (2*flt._lvl+1)*(2*flt._lvl+1)*(2*flt._lvl+1); //NON_ZERO ESTIMATE FOR STRUCTURED GRID!
//
  t1 = high_resolution_clock::now();
  err = filt_init(n_e, els, nds, flt._H, flt._lvl, flt._nze);
  t2 = high_resolution_clock::now(); t21 = duration_cast<milliseconds>(t2-t1).count();
//
  printf(" filter matrix made in %7.1e seconds\n",(double)t21/1000.);
  printf(" estimated nnz per element: %d\n",flt._nze);
  printf(" reserved nnz: %d\n",flt._nze*n_e);
  printf(" actual nnz: %ld\n",flt._H.nonZeros()); 
  printf(" filled: %4.2f\n",(double)flt._H.nonZeros()/n_e/n_e);
//
//opt. init.
//
  opt._ef = VectorXd::Zero(2);
  opt._efh = 1e8*VectorXd::Ones(2);
  opt._df = MatrixXd::Zero(2,n_e);
  opt._cf = MatrixXd::Zero(2, n_e);
  opt._xl = VectorXd::Zero(n_e);
  opt._xu = VectorXd::Ones(n_e);
  opt._mv = VectorXd::Ones(n_e)*0.1;
  opt._cs = VectorXi::Ones(1);
  opt._la = VectorXd::Zero(1);
  opt._dxl = VectorXd::Zero(n_e);
  opt._dxu = VectorXd::Zero(n_e);
  opt._cs = VectorXi::Ones(1);
//
//acceptable subproblem init.
//
  acc._ex = opt._ex;
  acc._ef = opt._efh;
  acc._df = opt._df;
  acc._cf = opt._cf;
  acc._dxl = opt._dxl;
  acc._dxu = opt._dxu;
  acc._la = opt._la;
  double a = 1.0;
//
//linear system init.
//
//dof vectors containing 
// - loads
// - prescribed dofs, and 
// - (eventually) the solution vector
//
  SpMat ffs_K(n_f,n_f);
  dfs_lds.setZero(n_d); 
  dfs_pre.setZero(n_d); 
  dfs_sol.setZero(n_d);
  dfs_lds(lds.col(0)) = lds.col(1);
  dfs_pre(bds.col(0)) = VectorXi::Ones(n_p); 
//
//free DOF vectors for loads and solution vector
  ffs_lds.setZero(n_f); 
  ffs_sol.setZero(n_f);
//
//mapping dof numbering vector with prescribed dofs removed
//
  c=0;
  map_num.setZero(n_d);
  for(i  = 0; i < n_d ; i++){
    if(dfs_pre(i)){
      map_num(i) = -1;
    }else{
      map_num(i) = c;
      c=c+1;
    }
  }
//
//make aux and/or aux struct with nice generic names (tmp1, vec1)
//
  double opt_def = 1.0;
  double opt_scl;
  double opt_obj;
  VectorXd opt_tmp; 
  VectorXd opt_sns1; 
  VectorXd opt_sns2; 
//
//if finite differences set x to 0.5
//
  double fdx = 0;
  int fdi = (int) (n_e / 2);
  if(fdc){
    opt._ex = 0.5*VectorXd::Ones(n_e);
  }
//
//print headers 
//
  if(fdc==0){
    printf("======================================================================\n");
    printf("  k         f0         f1       dx       df       bw       tk       tt\n");
  }else{
    printf("=======================================\n");
  }
//
//opt loop
//
  for(int itr=0; itr<400; itr++){
//
//  time
//
    auto ti = high_resolution_clock::now();
//
//  filter the design variables
//
    flt._ro = flt._H*opt._ex;
//
//  free dof load vector
//
    for(i  = 0; i < n_d ; i++){
      if(dfs_pre(i)){
      }else{
        ffs_lds(map_num(i)) = dfs_lds(i);
      }
    }
//
//  assemble with filtered density in the loop
//
    err=assy(n_e,n_f,nnz,els,nds,dfs_pre,map_num,flt._ro,ffs_K,ffs_lds);
//
//  solve
//
    err=solv(sol, n_f, ffs_K, ffs_lds, ffs_sol);
//
//  pack solution vector into complte dof vector
//
    c = 0;
    for(i  = 0; i < n_d ; i++){
        if(dfs_pre(i)){
        dfs_sol(i)=0.;
      }else{
        dfs_sol(i) = ffs_sol(c);
        c=c+1;
      }
    }
//
//  compute compliance objective and sensitivities (similar to assembly)
//
    opt_obj = 0.;
    opt_tmp = VectorXd::Zero(n_e);
    opt_sns1 = VectorXd::Zero(n_e);
    err=sens(n_d,els,nds,dfs_sol,dfs_pre,flt._ro,opt_tmp,&opt_obj);
//
//  filter the sensitivities (need to check if transpose multiplication is fast...)
//
    opt_sns1 = flt._H.transpose()*opt_tmp;
    opt_tmp = -VectorXd::Ones(n_e)/n_e/opt_vf;
    opt_sns2 = flt._H.transpose()*opt_tmp;
//
//  scaling objective and its sensitivities
//
    if(itr==0){
      opt_scl = opt_obj;
    }
    opt_obj=opt_obj/opt_scl;
    opt_sns1=opt_sns1/opt_scl;
//
    opt._ef(0) = opt_obj;
    opt._ef(1) = -flt._ro.sum()/n_e/opt_vf + 1.;
//
//  calculate delta f
//
    opt_def = opt._ef(0) - acc._ef(0);
//
//  check if feasible descent step
//
    if(opt_def > 1e9 && fdc==0){
//
//    retrieve the last acceptable iterate and objective value / do not update subproblem
//
      opt._ex = acc._ex;
      opt._ef = acc._ef; // just for output
      a=2.0*a;
//
    }else{
//
//    construct subproblem
//
      if(fdc){
        if(itr == 0){
          fdc_ef0=opt._ef;
          opt._df(0,all) = opt_sns1;
          opt._df(1,all) = opt_sns2;
        }else{
          VectorXd opt_err = VectorXd::Zero(2);
          printf("%8.1e ", fdx);
          for(int j = 0; j<2; j++){
            opt_err(j) = ((opt._ef(j)-fdc_ef0(j))/fdx - opt._df(j,fdi))/opt._df(j,fdi);
            printf("%14.7e ", opt_err(j));
          }
          printf("\n");
        }
      }else{
//
//      spherical curvature
//
        double tmp; double nrm;
        opt._df(0,all) = opt_sns1;
        opt._df(1,all) = opt_sns2;
        opt._cf(0,all) = eps*VectorXd::Ones(n_e);
        if(itr>0){
          nrm = (opt._ex-opt._exh).squaredNorm();
		  if(nca){
            tmp=2.*fabs(opt._efh(0)-opt._ef(0)-opt._df(0,all).dot(opt._exh-opt._ex))/nrm;
		  }else{
            tmp=2.*(opt._efh(0)-opt._ef(0)-opt._df(0,all).dot(opt._exh-opt._ex))/nrm;
	      }
          opt._cf(0,all) = (tmp*VectorXd::Ones(n_e)).cwiseMax(eps);
//        opt._cf(0,all) = -2.0*opt_sns.cwiseProduct(opt._ex.cwiseInverse());
        }
        opt._cf(1,all) = VectorXd::Zero(n_e);
//
//      subproblem bounds
//
        opt._dxl = opt._xl.cwiseMax(opt._ex - 0.1*(opt._xu-opt._xl));
        opt._dxu = opt._xu.cwiseMin(opt._ex + 0.1*(opt._xu-opt._xl));
//
//      save this acceptable iterate
//
        acc._ex = opt._ex; acc._df = opt._df; acc._cf = opt._cf;
        acc._dxl = opt._dxl; acc._dxu = opt._dxu; acc._la = opt._la;
        opt._exh=opt._ex; acc._ef = opt._ef; a=1.0;
      }
    }
//
//  solve subproblem and update design variables
//
    double opt_dx = 0.;
    VectorXd opt_nx = VectorXd::Zero(n_e);
    if(fdc){
      fdx = 1e-16*pow(10.0,itr);
      opt._ex = 0.5*VectorXd::Ones(n_e);
      opt._ex(fdi) = opt._ex(fdi) + fdx;
    }else{
      MatrixXd tmp_cf = acc._cf*a;
      err=dqpsub(acc._ex, acc._ef, acc._df, tmp_cf, acc._dxl, acc._dxu, opt._cs, acc._la, opt_nx);
      opt_dx = (opt_nx - opt._ex).array().abs().maxCoeff();
      opt._ex = opt_nx; opt._efh = opt._ef;
    }
//
//  write output
//
    wrte_lvtk(argv[3],itr,nds,els,gss,mat,mat_aux,dfs_sol,dfs_pre,dfs_lds,opt._ex,flt._ro);
//
//  time, black and white, and screen output
//
    t2 = high_resolution_clock::now();
    auto dur_itr = duration_cast<milliseconds>(t2 - ti).count();
    auto dur_tot = duration_cast<milliseconds>(t2 - t0).count();
    if(fdc){
      if(itr==18){
        break;
      }
    }else{
      double opt_baw = 0.0;
      for(int j = 0; j<n_e;j++){
        if( opt._ex(j) < opt._xl(j) + 1e-3 ){
          opt_baw = opt_baw + 1.0/n_e;
        }else if ( opt._ex(j) > opt._xu(j) - 1e-3) {
          opt_baw = opt_baw + 1.0/n_e;
        }
      }
      printf("%3d %10.3e %10.3e %8.1e %8.1e %8.1e %8.1e %8.1e\n", itr, opt._ef(0)*opt_scl, 
        opt._ef(1), opt_dx, opt_def, opt_baw, (double)dur_itr/1000, (double)dur_tot/1000.);
    }
    fflush(stdout);
//
//  possibly terminate
//
    if(fdc==0 && opt_dx < 1e-2 && a < 1.5 && fabs(opt_def) < 1e-6 ){
      break;
    }
//
  }
//
}
