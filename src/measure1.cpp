// SPDX-License-Identifier: LGPL-3.0-only
/*
* Author: Hao-Xin Wang <wanghaoxin1996@gmail.com>
* Creation Date: 15th, Jan, 2024
*
* Description: Measurement of the one-point function
* Usage: mpirun -n 1 ./measure1 params.json
 * or ./measure1 params.json
*/

///<  TODO: optimize the CPU cost.


#include "gqten/gqten.h"
#include <ctime>
#include "gqdouble.h"
#include "operators.h"
#include "params_case.h"
#include "myutil.h"
#include "my_measure.h"


using namespace gqmps2;
using namespace gqten;
using namespace std;

using FiniteMPST = gqmps2::FiniteMPS<TenElemT, U1U1QN>;

int main(int argc, char *argv[]) {
//  namespace mpi = boost::mpi;
//  mpi::environment env;
//  mpi::communicator world;

  CaseParams params(argv[1]);
  size_t Lx = params.Lx;
  size_t N = 2 * Lx * params.Ly;

  clock_t startTime, endTime;
  startTime = clock();

  gqten::hp_numeric::SetTensorTransposeNumThreads(params.TotalThreads);
  gqten::hp_numeric::SetTensorManipulationThreads(params.TotalThreads);

  OperatorInitial();

  const SiteVec<TenElemT, U1U1QN> sites = SiteVec<TenElemT, U1U1QN>(N, pb_out);

  FiniteMPST mps(sites);

  Timer one_site_timer("measure one site operators");
  MeasureOneSiteOp(mps, kMpsPath, {nf, sz}, {"nf", "sz"});
  cout << "measured one point function.<====" << endl;
  one_site_timer.PrintElapsed();

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}
