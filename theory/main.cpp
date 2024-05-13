#include <iostream>
#include <random>
#include <map>
#include <vector>
#include <set>
#include <cmath>
#include <numbers>
using namespace std;

double outdegree_plus_one(double x){
  double g;
  double outDegree = pow(x, 10);
  g = x * outDegree;
  return g;
}

double bernoulli(double x, double p0){
  double f = 0;
  f += p0;
  f += (1 - p0) * outdegree_plus_one(x);
  return f;
}
//Needed to get fbar at each time point
double dsdt(double si, double fi, double pi, double mu, double fbar, double z){
  double dsdt1 = 0;
  dsdt1 += (1 - mu)*z*(1 - si); 
  dsdt1 += mu*pi*z*(1-si); 
  dsdt1 += -1*mu*(1 - pi)*(1 + z*si);
  dsdt1 *= (si * fi);

  double dsdt2 = 0;
  dsdt2 += -1*(1- mu)*z*si;
  dsdt2 += mu*pi*(1 + z*(1-si));
  dsdt2 += -1*mu*(1 - pi)*z*si;
  dsdt2 *= (fbar - (si*fi));

  return dsdt1 + dsdt2;
}

double dqdt(double qi, double ai, double p0i){
  double dqdt = 0;
  dqdt += bernoulli(qi, p0i);
  dqdt -= qi;
  dqdt *= ai;
  return dqdt;
}
//I want to have a function where I can get just one si curve, or all of them

//double polylog(complex<double> x, int kMin){
//  double pLog = 0;
//  for(int i = kMin; i != int(1e4); ++i){
//    pLog += pow(x, i)/pow(i, kappa);
//  }
//  return pLog;
//}



int main(){
  vector<double> fitness = {0.9, 1.1};
  vector<vector<double>> Q = {{0},{0}};
  vector<vector<double>> S = {{0.5},{0.5}};
  vector<double> p0 = {0,0};
  vector<double> a = {0,0};
  vector<double> FQ = {0,0};
  double maxTime = 10;
  double time = 0;
  double dt = 0.01;
  double fbar;
  double mu = 0;
  double z = 10;
  
  while(time < maxTime){
    fbar = 0;
    for(int i  = 0; i != fitness.size(); ++i){
      fbar += fitness[i] * S[i].back();
    }
    for(int i  = 0; i != fitness.size(); ++i){
      S[i].push_back(S[i].back() + (dt * dsdt(S[i].back(), fitness[i], 0.5, mu, fbar, z)));
      a[i] = z*fbar + fitness[i];
      p0[i] = (z*fbar + fitness[i]) * mu/a[i];
      FQ[i] = bernoulli(Q[i].back(), p0[i]);
      Q[i].push_back(Q[i].back() + (dt * dqdt(Q[i].back(), a[i], p0[i])));
    }
  }

  for(int i = 0; i != fitness.size(); ++i){
    for(int j = 0; j != Q[0].size(); ++j){
      cout << 0.01 * j << ',' << Q[i][j] << ", q" << i << endl;
      cout << 0.01 * j << ',' << S[i][j] << ", s" << i << endl;
    }
  }

}
