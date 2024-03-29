#pragma once
#include "qlmps/qlmps.h"
using qlmps::CaseParamsParserBasic;


struct CaseParams : public CaseParamsParserBasic {
  CaseParams(const char *pf) : CaseParamsParserBasic(pf) {
    Lx = ParseInt("Lx");
    Ly = ParseInt("Ly");
    ts = ParseDouble("ts");
    td = ParseDouble("td");
    tsd_xy = ParseDouble("tsd_xy");
    tsd_nn = ParseDouble("tsd_nn");
    Uss = ParseDouble("Uss");
    Usd = ParseDouble("Usd");
    Udd = ParseDouble("Udd");
    noise = ParseDoubleVec("noise");
    Numhole = ParseInt("Numhole");
    Sweeps = ParseInt("Sweeps");
    Dmin = ParseInt("Dmin");
    Dmax = ParseInt("Dmax");
    CutOff = ParseDouble("CutOff");
    LanczErr = ParseDouble("LanczErr");
    MaxLanczIter = ParseInt("MaxLanczIter");
    TotalThreads = ParseInt("TotalThreads");
    Perturbation = ParseBool("Perturbation");
    if (Perturbation) {
      PA = ParseDouble("PerturbationAmplitude");
      PerturbationPeriod = ParseInt("PerturbationPeriod");
    } else {
      PA = 0.0;
      PerturbationPeriod = 1;
    }
  }

  size_t Lx;
  size_t Ly;
  double ts;
  double td;
  double tsd_xy;//diag
  double tsd_nn;//hopping between s-orbital and d_x^2-y^2 orbital
  double Uss;
  double Udd;
  double Usd;
  size_t Sweeps;
  size_t Dmin;
  size_t Dmax;
  int Numhole;
  double CutOff;
  double LanczErr;
  size_t MaxLanczIter;
  size_t TotalThreads;
  std::vector<double> noise;
  bool Perturbation;
  double PA;
  size_t PerturbationPeriod;
};
