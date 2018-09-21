
#include <vector>

#include "main.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// estimate f by maximum likelihood
// [[Rcpp::export]]
Rcpp::List estimate_f_cpp(Rcpp::List args) {
  
  // extract inputs
  vector<vector<int>> x = rcpp_to_mat_int(args["x"]);
  vector<double> p = rcpp_to_vector_double(args["p"]);
  bool report_progress = rcpp_to_bool(args["report_progress"]);
  int f_breaks = rcpp_to_int(args["f_breaks"]);
  
  // get/set dimensions
  int n = x.size();
  int L = x[0].size();
  
  // create q = 1-p
  vector<double> q(L);
  for (int j=0; j<L; ++j) {
    q[j] = 1.0 - p[j];
  }
  
  // create lookup tables
  vector<vector<double>> lookup_homo1(f_breaks, vector<double>(L));
  vector<vector<double>> lookup_homo2(f_breaks, vector<double>(L));
  vector<vector<double>> lookup_het(f_breaks, vector<double>(L));
  for (int k=0; k<f_breaks; ++k) {
    double f = k/double(f_breaks-1);
    for (int j=0; j<L; ++j) {
      lookup_homo1[k][j] = log((1-f)*p[j]*p[j] + f*p[j]);
      lookup_homo2[k][j] = log((1-f)*q[j]*q[j] + f*q[j]);
      lookup_het[k][j] = log((1-f)*2*p[j]*q[j]);
    }
  }
  
  // create objects for storing results
  vector<double> loglike_vec(f_breaks);
  vector<vector<double>> ret(n, vector<double>(n));
  
  // loop through all pairwise samples
  for (int i1=0; i1<(n-1); ++i1) {
    
    // report progress
    if (report_progress) {
      print("sample", i1, "of", n);
    }
    
    for (int i2=(i1+1); i2<n; ++i2) {
      
      // calculate loglike for every value of f
      for (int k=0; k<f_breaks; ++k) {
        double loglike = 0;
        for (int j=0; j<L; ++j) {
          if (x[i1][j] == -1 || x[i2][j] == -1) {
            continue;
          }
          if (x[i1][j] == 1 && x[i2][j] == 1) {
            loglike += lookup_homo1[k][j];
          } else if (x[i1][j] == 0 && x[i2][j] == 0) {
            loglike += lookup_homo2[k][j];
          } else {
            loglike += lookup_het[k][j];
          }
        }
        loglike_vec[k] = loglike;
      }
      
      // store maximum likelihood f
      double best_loglike = loglike_vec[0];
      for (int k=1; k<f_breaks; ++k) {
        if (loglike_vec[k] > best_loglike) {
          best_loglike = loglike_vec[k];
          ret[i1][i2] = k/double(f_breaks-1);
        }
      }
      
    }  // end i2 loop
  }  // end i1 loop
  
  // return as Rcpp object
  return Rcpp::List::create(Rcpp::Named("ret") = ret);
}

//------------------------------------------------
// calculate proportion IBS
// [[Rcpp::export]]
Rcpp::List calculate_IBS_cpp(Rcpp::List args) {
  
  // extract inputs
  vector<vector<int>> x = rcpp_to_mat_int(args["x"]);
  bool report_progress = rcpp_to_bool(args["report_progress"]);
  
  // get dimensions
  int n = x.size();
  int L = x[0].size();
  
  // create objects for storing results
  vector<vector<double>> ret(n, vector<double>(n));
  
  // loop through all pairwise samples
  int IBS_num = 0;
  int IBS_denom = 0;
  for (int i1=0; i1<(n-1); ++i1) {
    
    // report progress
    if (report_progress) {
      print("sample", i1, "of", n);
    }
    
    for (int i2=(i1+1); i2<n; ++i2) {
      
      // calculate proportion identical in state
      IBS_num = 0;
      IBS_denom = 0;
      for (int j=0; j<L; ++j) {
        if (x[i1][j] == -1 || x[i2][j] == -1) {
          continue;
        }
        if (x[i1][j] == x[i2][j]) {
          IBS_num++;
        }
        IBS_denom++;
      }
      ret[i1][i2] = IBS_num/double(IBS_denom);
      
    }  // end i2 loop
  }  // end i1 loop
  
  // return as Rcpp object
  return Rcpp::List::create(Rcpp::Named("ret") = ret);
}
