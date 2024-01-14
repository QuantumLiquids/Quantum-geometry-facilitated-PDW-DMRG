/**
 * DMRG for Huang PDW model
 *
 */
#include "gqmps2/gqmps2.h"
#include "./gqdouble.h"
#include "./operators.h"
#include "./myutil.h"
#include "./params_case.h"

using FiniteMPST = gqmps2::FiniteMPS<TenElemT, U1U1QN>;

int main(int argc, char *argv[]) {
  using namespace std;
  using namespace gqmps2;
  using namespace gqten;
  namespace mpi = boost::mpi;
  mpi::environment env;
  mpi::communicator world;
  CaseParams params(argv[1]);

  if (world.rank() == 0 && world.size() > 1 && params.TotalThreads > 2) {
    gqten::hp_numeric::SetTensorTransposeNumThreads(params.TotalThreads - 2);
    gqten::hp_numeric::SetTensorManipulationThreads(params.TotalThreads - 2);
  } else {
    gqten::hp_numeric::SetTensorTransposeNumThreads(params.TotalThreads);
    gqten::hp_numeric::SetTensorManipulationThreads(params.TotalThreads);
  }

  /******** Model parameter ********/
  size_t Lx = params.Lx, Ly = params.Ly;
  size_t N = 2 * Lx * Ly;//two orbital
  double ts = params.ts, td = params.td, tsd_nn = params.tsd_nn, tsd_xy = params.tsd_xy;
  double Uss = params.Uss, Udd = params.Udd, Usd = params.Usd;
  cout << "System size = (" << Lx << "," << Ly << ")" << endl;
  cout << "The number of electron sites =" << Lx * Ly << endl;
  cout << "Model parameter: ts :" << ts << ", Uss :" << Uss
       << endl;

  /****** DMRG parameter *******/
  gqmps2::FiniteVMPSSweepParams sweep_params(
      params.Sweeps,
      params.Dmin, params.Dmax, params.CutOff,
      gqmps2::LanczosParams(params.LanczErr, params.MaxLanczIter),
      params.noise
  );

  clock_t startTime, endTime;
  startTime = clock();
  OperatorInitial();
  const SiteVec<TenElemT, U1U1QN> sites = SiteVec<TenElemT, U1U1QN>(N, pb_out);

  /****** Initialize MPS ******/
  FiniteMPST mps(sites);
  if (world.rank() == 0) {
    if (!IsPathExist(kMpsPath) || !(N == GetNumofMps())) {
      cout << "Initial mps as direct product state." << endl;
      std::vector<size_t> stat_labs(N, 0);
      size_t sitenumber_perhole;
      if (params.Numhole > 0) {
        sitenumber_perhole = N / params.Numhole;
      } else {
        sitenumber_perhole = 4 * N;
      }

      size_t half_filled_qn_label = 1;
      for (size_t i = 0; i < N; i++) {
        if (i % sitenumber_perhole == sitenumber_perhole / 2) {
          stat_labs[i] = 3;   //punch a hole
        } else {
          stat_labs[i] = half_filled_qn_label;
          half_filled_qn_label = 3 - half_filled_qn_label;
        }
      }
      gqmps2::DirectStateInitMps(mps, stat_labs);
      mps.Dump(sweep_params.mps_path, true);
    }
  }


  /*******  Calculation Subroutine *******/
  // creation the mpo/mro by the newest phonon displacement data
  gqmps2::MPOGenerator<TenElemT, U1U1QN> mpo_gen(sites, qn0);
  for (size_t i = 0; i < N; i = i + 2) {
    mpo_gen.AddTerm(-Uss, Uterm, i);
    mpo_gen.AddTerm(-Udd, Uterm, i + 1);
    mpo_gen.AddTerm(Usd, cupccdnc, i, cdnacupa, i + 1);
    mpo_gen.AddTerm(Usd, cdnacupa, i, cupccdnc, i + 1);
  }

  //horizontal hopping inner band
  for (size_t i = 0; i < N - 2 * Ly; i = i + 2) {
    size_t site1 = i, site2 = i + 2 * Ly;
    size_t x = i / (2 * Ly);
    mpo_gen.AddTerm(-ts, bupcF, site1, bupa, site2, f);
    mpo_gen.AddTerm(-ts, bdnc, site1, Fbdna, site2, f);
    mpo_gen.AddTerm(ts, bupaF, site1, bupc, site2, f);
    mpo_gen.AddTerm(ts, bdna, site1, Fbdnc, site2, f);

    site1++;
    site2++;
    mpo_gen.AddTerm(-td, bupcF, site1, bupa, site2, f);
    mpo_gen.AddTerm(-td, bdnc, site1, Fbdna, site2, f);
    mpo_gen.AddTerm(td, bupaF, site1, bupc, site2, f);
    mpo_gen.AddTerm(td, bdna, site1, Fbdnc, site2, f);
  }

  //vertical hopping inner band
  for (size_t i = 0; i < N; i += 2) {
    size_t y = i % (2 * Ly), x = i / (2 * Ly);
    if (y < 2 * Ly - 2) {
      size_t site1 = i, site2 = i + 2;
      mpo_gen.AddTerm(-ts, bupcF, site1, bupa, site2);
      mpo_gen.AddTerm(-ts, bdnc, site1, Fbdna, site2);
      mpo_gen.AddTerm(ts, bupaF, site1, bupc, site2);
      mpo_gen.AddTerm(ts, bdna, site1, Fbdnc, site2);
      cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;

      site1++;
      site2++;
      mpo_gen.AddTerm(-td, bupcF, site1, bupa, site2);
      mpo_gen.AddTerm(-td, bdnc, site1, Fbdna, site2);
      mpo_gen.AddTerm(td, bupaF, site1, bupc, site2);
      mpo_gen.AddTerm(td, bdna, site1, Fbdnc, site2);
      cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;
    } else if (Ly > 2) {
      size_t site1 = i - 2 * Ly + 2, site2 = i;
      mpo_gen.AddTerm(-ts, bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(-ts, bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(ts, bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(ts, bdna, site1, Fbdnc, site2, f);
      cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;

      site1++;
      site2++;
      mpo_gen.AddTerm(-ts, bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(-ts, bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(ts, bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(ts, bdna, site1, Fbdnc, site2, f);
      cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;

    }
  }

  //intra band hopping
  for (size_t i = 0; i < N - (2 * Ly); i += 2) {
    size_t y = i % (2 * Ly), x = i / (2 * Ly);
    if (y < 2 * Ly - 2) {
      size_t site1 = i, site2 = i;
      mpo_gen.AddTerm(-ts, bupcF, site1, bupa, site2);
      mpo_gen.AddTerm(-ts, bdnc, site1, Fbdna, site2);
      mpo_gen.AddTerm(ts, bupaF, site1, bupc, site2);
      mpo_gen.AddTerm(ts, bdna, site1, Fbdnc, site2);
      cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;

      site1++;
      site2++;
      mpo_gen.AddTerm(-td, bupcF, site1, bupa, site2);
      mpo_gen.AddTerm(-td, bdnc, site1, Fbdna, site2);
      mpo_gen.AddTerm(td, bupaF, site1, bupc, site2);
      mpo_gen.AddTerm(td, bdna, site1, Fbdnc, site2);
      cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;
    } else if (Ly > 2) {

    }

  }

  auto mro = mpo_gen.GenMatReprMPO();
  cout << "MRO generated." << endl;

  // dmrg
  double e0 = gqmps2::FiniteDMRG(mps, mro, sweep_params, world);

  return 0;
}
