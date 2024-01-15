//
// Created by haoxinwang on 15th, Jan 2024.
//

/*
    measure2.cpp
    for measure spin, charge, on-site pair, single-particle correlation function.
    memory optimized and parallel version.
    usage:
        mpirun -n 2*Ly ./measure2
    Optional arguments:
      --start=
      --end=
    Which are set as start=Lx/4, end = 3*Lx/4+2 by default
*/

#include "gqmps2/gqmps2.h"
#include "gqten/gqten.h"
#include <ctime>
#include "gqdouble.h"
#include "operators.h"
#include "params_case.h"
#include "myutil.h"
#include "my_measure.h"
#include "gqten/utility/timer.h"

#include "boost/mpi.hpp"

using std::cout;
using std::endl;
using std::vector;
using QNT = U1U1QN;
using FiniteMPST = gqmps2::FiniteMPS<TenElemT, QNT>;
using gqmps2::SiteVec;
using gqmps2::MeasureTwoSiteOp;
using gqten::Timer;
using gqmps2::MeasureGroupTask;
using gqmps2::kMpsPath;

int main(int argc, char *argv[]) {
  namespace mpi = boost::mpi;
  mpi::environment env;
  mpi::communicator world;
  if (argc == 1) {
    std::cout << "Usage: \n mpirun -np 2*Ly ./measure2 <params file> --start=<start x coor> --end=<end x coor>\n";
    return 0;
  } else if (argc == 2) {
    std::cout << "No start and end parameters. set start as Lx/4 and end as 3*Lx/4" << std::endl;
    std::cout
        << "The complete usage is: \n mpirun -np 2*Ly ./measure2 <params file> --start=<start x coor> --end=<end x coor>\n";
  }

  clock_t startTime, endTime;
  startTime = clock();

  size_t beginx;
  size_t endx;
  bool start_argument_has = ParserMeasureSite(argc, argv, beginx, endx);

  CaseParams params(argv[1]);

  size_t Lx = params.Lx, Ly = params.Ly;
  size_t N = 2 * Lx * Ly;
  if (GetNumofMps() != N) {
    std::cout << "The number of mps files are inconsistent with mps size!" << std::endl;
    exit(1);
  }

  if (!start_argument_has) {
    beginx = Lx / 4;
    endx = beginx + Lx / 2 + 2;
  }

  OperatorInitial();

  const SiteVec<TenElemT, QNT> sites = SiteVec<TenElemT, QNT>(N, pb_out);
  FiniteMPST mps(sites);
  gqten::hp_numeric::SetTensorTransposeNumThreads(params.TotalThreads);
  gqten::hp_numeric::SetTensorManipulationThreads(params.TotalThreads);

  Timer two_site_timer("measure two site operators");

  std::vector<MeasureGroupTask> measure_tasks;
  measure_tasks.reserve(N);
  //for structure factor.
//  for (size_t site1 = beginx * Ly; site1 < endx * Ly; site1++) {
//    std::vector<size_t> site2;
//    site2.reserve(N - site1);
//    for (size_t j = site1 + 1; j < endx * Ly; j++) {
//      site2.push_back(j);
//    }
//    measure_tasks.push_back(MeasureGroupTask(site1, site2));
//  }

  for (size_t y = 0; y < (2 * Ly); ++y) {
    auto site1 = beginx * (2 * Ly) + y;
    std::vector<size_t> site2;
    site2.reserve((2 * Ly) * (endx - beginx));
    for (size_t j = (beginx + 2) * (2 * Ly); j < endx * (2 * Ly); j++) {
      site2.push_back(j);
    }
    measure_tasks.push_back(MeasureGroupTask(site1, site2));
  }

  Timer two_site_measure_timer("measure spin_ structure factors");
  MeasureTwoSiteOp(mps, kMpsPath, sz, sz, measure_tasks, "zzsf", world);
  if (world.rank() == 0) {
    std::cout << "measured sz sz correlation." << std::endl;
  }
  world.barrier();
  MeasureTwoSiteOp(mps, kMpsPath, sp, sm, measure_tasks, "pmsf", world);
  MeasureTwoSiteOp(mps, kMpsPath, sm, sp, measure_tasks, "mpsf", world);
  MeasureTwoSiteOp(mps, kMpsPath, nf, nf, measure_tasks, "nfnf", world);
  MeasureTwoSiteOp(mps, kMpsPath, bupc, bupa, measure_tasks, "bupcbupa",
                   world, f);  //directly equal to the single-particle correlatin function
  MeasureTwoSiteOp(mps, kMpsPath, bdnc, bdna, measure_tasks, "bupcbupa",
                   world, f);  //directly equal to the single-particle correlatin function
  MeasureTwoSiteOp(mps,
                   kMpsPath,
                   cupccdnc,
                   cdnacupa,
                   measure_tasks,
                   "onsitepair",
                   world);//directly equal on-site pair correlation
  two_site_measure_timer.PrintElapsed();

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

  return 0;
}


