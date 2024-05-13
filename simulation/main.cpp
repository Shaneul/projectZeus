#include <iostream>
#include <random>
#include <map>
#include <vector>
#include <set>
#include <cmath>
#include <numbers>
using namespace std;



static double riemann_zeta(double gamma){
  double zeta = 0;
  for (int i = 1; i < 10000000; ++i){
    zeta += pow(i, -gamma); 
  }
  return zeta;
}

//function to generate a powerlaw cdf
vector<double> powerlawCdf(int kMax, double gamma){
  vector<double> cdf;
  map<int, double> pdf;
  double zeta = riemann_zeta(gamma);
  for (int i = 1; i < kMax; ++i){
    pdf[i] = pow(i, -gamma)/zeta;
    cdf.push_back(0);
  }
  for (int itr = 0; itr != cdf.size(); ++itr){
    for (int i = 1; i < kMax; ++i){
      if (i <= itr){
        cdf[itr] += pdf[i];
      }
    }
  }
  return cdf;
}


// uses inverse transform method to generate a degree sequence from a given cdf with specified kmin and kmax
void Graph::powerlawDegreeList(int kMin, int kMax){
  vector<double> cdf = powerlawCdf(kMax, kappa);
  mt19937 gen(random_device{}()); 
  uniform_real_distribution<double> unif(cdf[kMin-1], cdf[kMax]);

  while (degreeList.size() < N){
    double p = unif(gen);
    for (int i = 0; i != kMax; ++i){
      if(p > cdf[i]){
        if(p <= cdf[i+1]){
          degreeList.push_back(i+1);
        }
      }
    }
  }
  //return degreeList;
void Graph::generateGraph(){
  
  /**
  Function to take the desired in and out degree distributions for a Graph and generate the neighbor sets
  */
  for(int node = 0; node != N; ++node){
    neighbors[node] = set<int> {};
  }
  set<string> allowedOutDistributions = {"regular", "powerlaw"};
  set<string> allowedInDistributions = {"regular", "poisson"};
  
  if(allowedOutDistributions.find(outDegreeDistribution) == allowedOutDistributions.end()){
    cerr << "Out degree distribution not defined" << endl;
  }else if(allowedInDistributions.find(inDegreeDistribution) == allowedInDistributions.end()){
    cerr << "In degree distribution not defined" << endl;
  }else{
    //degree distribution choices are allowed, genrate the graph
    mt19937 gen(random_device{}());
    uniform_int_distribution<int> unif(0, N-1);
    if(outDegreeDistribution == "regular" && inDegreeDistribution == "regular"){
      vector<int> stubsIn;
      vector<int> stubsOut;
      for(int i = 0; i != N; ++i){
        for(int j = 0; j != kMean; ++j){
          stubsIn.push_back(i);
          stubsOut.push_back(i);
        }
      }
      shuffle(stubsIn.begin(), stubsIn.end(), gen);
      shuffle(stubsOut.begin(), stubsOut.end(), gen);
      for(int i = 0; i != stubsIn.size(); ++i){
        int v1 = stubsOut[i];
        int v2 = stubsIn[i];
        neighbors[v1].insert(v2);
      }

    }else if(outDegreeDistribution == "powerlaw" && inDegreeDistribution == "poisson"){
      powerlawDegreeList(4, N);
      for(int i = 0; i != N; ++i){
        while(neighbors[i].size() < degreeList[i]){
          int neighbor = unif(gen);
          if(neighbor != i){
            neighbors[i].insert(neighbor);
          }
        }
      }      
    }else{
      powerlawDegreeList(4, N);
      vector<int> stubsOut = degreeList;
      vector<int> stubsIn;
      int degreeSum = 0;
      for(int i = 0; i != N; ++i){
        degreeSum += degreeList[i];
      }
      int nodeIndex = 0;
      while(degreeSum > kMean * N){
        for(int j = 0; j != kMean + 1; ++j){
          stubsIn.push_back(nodeIndex);
        }
        nodeIndex += 1;
        degreeSum -= kMean + 1;
      }
      while(degreeSum % kMean != 0){
        for(int j = 0; j != kMean + 1; ++j){
          stubsIn.push_back(nodeIndex);
        }
        nodeIndex += 1;
        degreeSum -= kMean + 1;
      }
      for(int i = nodeIndex; i != N; ++i){
        for(int j = 0; j != kMean; ++j){
          stubsIn.push_back(i);
        }
      }
      shuffle(stubsIn.begin(), stubsIn.end(), gen);
      shuffle(stubsOut.begin(), stubsOut.end(), gen);
      for(int i = 0; i != stubsIn.size(); ++i){
        int v1 = stubsOut[i];
        int v2 = stubsIn[i];
        neighbors[v1].insert(v2);
      }
    }
  }
}
double Graph::dsdt(double si, double fi, double pi){
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


//I want to have a function where I can get just one si curve, or all of them

void Graph::eulersMethod(){
  double t0 = 0;
  double dt = 0.001;
  //calculate initial actual s values 
  for(int i = 0; i != inputMemeFitness.size(); ++i){
    double si = 0;
    for(int node = 0; node != N; ++node){
      if(memeFitness[meme[node]] == inputMemeFitness[i]){
        si += 1/double(N);
      } 
    }
    sOverTime.push_back(si);
  }
  while(t0 < maxTime){
    fbar = 0;
    for(int i = 0; i != sOverTime.size(); ++i){
      fbar += sOverTime[i] * inputMemeFitness[i];
    }
    for(int i = 0; i != sOverTime.size(); ++i){
      double ds = dt * dsdt(sOverTime[i], inputMemeFitness[i], inputFitnessPdf[i]);
      sOverTime[i] += ds;
      cout << t0 + dt << ',' << sOverTime[i] << ',' << "mf_s" << i << endl;
    }
    t0 += dt;
  }
}
void Graph::set_ri(complex<double> r){
  ri[0] = r;
}

complex<double> Graph::regular_pgf(complex<double> x){
  complex<double> f;
  f = pow(x, z);
  return f;
}

complex<double> Graph::bernoulli_pgf(complex<double> x){
  return bernoulli_pgf(x, 0);
}

complex<double> Graph::bernoulli_pgf(complex<double> x, int index){
  complex<double> unity = 1;
  complex<double> binPgf = unity - ri[index] + (x * ri[index]);
  return binPgf;
}

complex<double> Graph::polylog(complex<double> x, int kMin){
  complex<double> pLog = 0;
  for(int i = kMin; i != int(1e4); ++i){
    pLog += pow(x, i)/pow(i, kappa);
  }
  return pLog;
}

complex<double> Graph::powerlaw_pgf(complex<double> x){
  complex<double> numerator = polylog(x, kappa);
  return numerator/zeta;
}

complex<double> Graph::outdegree_pgf(complex<double> x){
  complex<double> f;
  if(outDegreeDistribution == "regular"){
    f = regular_pgf(x);
  }
  if(outDegreeDistribution == "powerlaw"){
    f = powerlaw_pgf(x);
  }
  return f;
}

complex<double> Graph::outdegree_plus_one(complex<double> x){
  complex<double> g;
  complex<double> outDegree = outdegree_pgf(x);
  g = x * outDegree;
  return g;
}

complex<double> Graph::semelparity(complex<double> x){
  return semelparity(x, 0);
}

complex<double> Graph::semelparity(complex<double> x, int index){
  complex<double> pgf = bernoulli_pgf(regular_pgf(x), index);
  return pgf;
}

complex<double> Graph::recursive_semelparity(complex<double> x, int g){
  return recursive_semelparity(x, g, 0);
}

complex<double> Graph::recursive_semelparity(complex<double> x, int g, int index){
  vector<complex<double>> Ts;
  Ts.push_back(x);
  for(int i = 0; i != g; ++i){
    Ts.push_back(x * semelparity(Ts.back(), index));
  }
  return Ts.back();
}

complex<double> Graph::popularity(complex<double> x, int g){
  return popularity(x, g, 0);
}

complex<double> Graph::popularity(complex<double> x, int g, int index){
  vector<complex<double>> omegas;
  omegas.push_back(x);
  for(int i = 0; i < g; ++i){
    omegas.push_back(x * outdegree_plus_one(bernoulli_pgf(omegas.back(), index)));
  }
  return omegas.back();
}
