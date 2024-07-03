/*
 * File Name: measure_nk.cpp
 * Description: measure static single particle correlation function for every two points
 * Created by Hao-Xin on 2024/07/02.
 *
 * usage:
 *     mpirun -n 5 ./measure_nk params.json
 * Processor number only determine speed
 */


#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <ctime>
#include "gqdouble.h"
#include "operators.h"
#include "params_case.h"
#include "myutil.h"
#include "my_measure.h"
#include "qlten/utility/timer.h"
#include "boost/serialization/complex.hpp"

using namespace qlmps;
using namespace qlten;
using namespace std;

int main(int argc, char *argv[]) {
  namespace mpi = boost::mpi;
  mpi::environment env;
  mpi::communicator world;

  CaseParams params(argv[1]);
  size_t Lx = params.Lx, Ly = params.Ly;
  size_t N = 2 * Lx * Ly;
  if (GetNumofMps() != N) {
    std::cout << "The number of mps files are inconsistent with mps size!" << std::endl;
    exit(1);
  }

  clock_t startTime, endTime;
  startTime = clock();

  qlten::hp_numeric::SetTensorManipulationThreads(params.TotalThreads);


  OperatorInitial();
  const SiteVec<TenElemT, U1U1QN> sites = SiteVec<TenElemT, U1U1QN>(N, pb_out);
  using FiniteMPST = qlmps::FiniteMPS<TenElemT, U1U1QN>;
  FiniteMPST mps(sites);
//  mps.Load();
//  cout << "mps loaded" << endl;
//  cout << "bond dimension of middle mps = ";
//  cout << mps[N / 2].GetShape()[0] << endl;


  /******** define the two site sets of two point function ********/
  vector<vector<size_t>> two_point_sites_setF;
  size_t beginx = 0;
  size_t endx = N;
  two_point_sites_setF.reserve((endx - beginx) * (endx - beginx));
  for (size_t i = 0; i < N; i++) {
    for (size_t j = i + 1; j < N; j++) {
      two_point_sites_setF.push_back({i, j});
    }
  }


  /******** define the measure_tasks ********/
  std::vector<MeasureGroupTask> measure_tasks;
  measure_tasks.reserve(N);

  for (size_t i = 0; i < N; i++) {
    const size_t site1 = i;
    std::vector<size_t> site2;
    site2.reserve(N - i);
    for (size_t j = i + 1; j < N; j++) {
      site2.push_back(j);
    }
    measure_tasks.push_back(MeasureGroupTask(site1, site2));
  }
  MeasureTwoSiteOp(mps, kMpsPath, bupcF, bupa, measure_tasks, "single_particle_correlation0", world, f);
  MeasureTwoSiteOp(mps, kMpsPath, bdnc, Fbdna, measure_tasks, "single_particle_correlation1", world, f);
  MeasureTwoSiteOp(mps, kMpsPath, bupaF, bupc, measure_tasks, "single_particle_correlation2", world, f);
  MeasureTwoSiteOp(mps, kMpsPath, bdna, Fbdnc, measure_tasks, "single_particle_correlation3", world, f);

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}