/**
 * DMRG for Huang PDW model
 *
 */
#include "qlmps/qlmps.h"
#include "./ql_hubbard_fermion_qn_only_spin_conservation.h"
#include "./operators.h"
#include "./myutil.h"
#include "./params_case.h"

using FiniteMPST = qlmps::FiniteMPS<TenElemT, U1QN>;

int main(int argc, char *argv[]) {
  using namespace std;
  using namespace qlmps;
  using namespace qlten;
  namespace mpi = boost::mpi;
  mpi::environment env;
  mpi::communicator world;

  if (argc == 1) {
    if (world.rank() == 0) {
      std::cout
          << "Usage: \n mpirun -np <num_proc> ./dmrg <params file> --D=<list of bond dimension, connected by comma>\n";
    }
    return 0;
  } else if (argc == 2) {
    if (world.rank() == 0)
      std::cout
          << "The complete usage can be: Usage: \n mpirun -np <num_proc> ./dmrg <params file> --D=<list of bond dimension, connected by comma>"
          << std::endl;
  }

  CaseParams params(argv[1], 1);

  if (world.rank() == 0 && world.size() > 1 && params.TotalThreads > 2) {
    qlten::hp_numeric::SetTensorManipulationThreads(params.TotalThreads - 2);
  } else {
    qlten::hp_numeric::SetTensorManipulationThreads(params.TotalThreads);
  }

  /******** Model parameter ********/
  size_t Lx = params.Lx, Ly = params.Ly;
  size_t N = 2 * Lx * Ly;//two orbital
  double ts = params.ts, td = params.td, tsd_nn = params.tsd_nn, tsd_xy = params.tsd_xy;
  double Uss = params.Uss, Udd = params.Udd, Usd = params.Usd;
  if (world.rank() == 0) {
    cout << "System size = (" << Lx << "," << Ly << ")" << endl;
    cout << "The number of electron sites =" << Lx * Ly << endl;
    cout << "Model parameter: ts :" << ts << ", Uss :" << Uss
         << endl;
  }
  /****** DMRG parameter *******/
  qlmps::FiniteVMPSSweepParams sweep_params(
      params.Sweeps,
      params.Dmin, params.Dmax, params.CutOff,
      qlmps::LanczosParams(params.LanczErr, params.MaxLanczIter),
      params.noise
  );

  clock_t startTime, endTime;
  startTime = clock();
  OperatorInitial();
  const SiteVec<TenElemT, U1QN> sites = SiteVec<TenElemT, U1QN>(N, pb_out);

  std::vector<size_t> input_D_set;
  bool has_bond_dimension_parameter = ParserBondDimension(
      argc, argv,
      input_D_set);
  size_t DMRG_time = input_D_set.size();
  std::vector<size_t> MaxLanczIterSet(DMRG_time);
  if (has_bond_dimension_parameter) {
    MaxLanczIterSet.back() = params.MaxLanczIter;
    if (DMRG_time > 1) {
      size_t MaxLanczIterSetSpace;
      MaxLanczIterSet[0] = 3;
      MaxLanczIterSetSpace = (params.MaxLanczIter - 3) / (DMRG_time - 1);
      if (world.rank() == 0)
        std::cout << "Setting MaxLanczIter as : [" << MaxLanczIterSet[0] << ", ";
      for (size_t i = 1; i < DMRG_time - 1; i++) {
        MaxLanczIterSet[i] = MaxLanczIterSet[i - 1] + MaxLanczIterSetSpace;
        if (world.rank() == 0)
          std::cout << MaxLanczIterSet[i] << ", ";
      }
      if (world.rank() == 0)
        std::cout << MaxLanczIterSet.back() << "]" << std::endl;
    } else {
      if (world.rank() == 0)
        std::cout << "Setting MaxLanczIter as : [" << MaxLanczIterSet[0] << "]" << std::endl;
    }
  }

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

      size_t half_filled_qn_label = 0;
      for (size_t i = 0; i < N; i++) {
        if (i % sitenumber_perhole == sitenumber_perhole / 2) {
          stat_labs[i] = 3;   //punch a hole
        } else {
          stat_labs[i] = half_filled_qn_label;
          half_filled_qn_label = 1 - half_filled_qn_label;
        }
      }
      qlmps::DirectStateInitMps(mps, stat_labs);
      mps.Dump(sweep_params.mps_path, true);
    }
  }


  /*******  Creation MPO/MRO *******/
  qlmps::MPOGenerator<TenElemT, U1QN> mpo_gen(sites, qn0);
  for (size_t i = 0; i < N; i = i + 2) {
    mpo_gen.AddTerm(-Uss, Uterm, i);
    mpo_gen.AddTerm(-params.mu_s, nf, i);
    mpo_gen.AddTerm(-Udd, Uterm, i + 1);
    mpo_gen.AddTerm(-params.mu_d, nf, i + 1);
    mpo_gen.AddTerm(Usd, cupccdnc, i, cdnacupa, i + 1);
    mpo_gen.AddTerm(Usd, cdnacupa, i, cupccdnc, i + 1);
  }

  //horizontal hopping inner band
  for (size_t i = 0; i < N - 2 * Ly; i = i + 2) {
    size_t site1 = i, site2 = i + 2 * Ly;
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
      mpo_gen.AddTerm(-ts, bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(-ts, bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(ts, bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(ts, bdna, site1, Fbdnc, site2, f);
      if (world.rank() == 0)
        cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;

      site1++;
      site2++;
      mpo_gen.AddTerm(-td, bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(-td, bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(td, bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(td, bdna, site1, Fbdnc, site2, f);
      if (world.rank() == 0)
        cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;
    } else if (Ly > 2) {
      size_t site1 = i - 2 * Ly + 2, site2 = i;
      mpo_gen.AddTerm(-ts, bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(-ts, bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(ts, bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(ts, bdna, site1, Fbdnc, site2, f);
      if (world.rank() == 0)
        cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;

      site1++;
      site2++;
      mpo_gen.AddTerm(-td, bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(-td, bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(td, bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(td, bdna, site1, Fbdnc, site2, f);
      if (world.rank() == 0)
        cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;
    }
  }

  //intra band diagonal hopping, tsd_xy
  for (size_t i = 0; i < N - (2 * Ly); i += 2) {
    size_t y = i % (2 * Ly), x = i / (2 * Ly);
    if (y < 2 * Ly - 2) {
      // left-upper-s <====> right-lower-d_xy
      size_t site1 = i, site2 = i + 3 + (2 * Ly);
      mpo_gen.AddTerm(-tsd_xy, bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(-tsd_xy, bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(tsd_xy, bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(tsd_xy, bdna, site1, Fbdnc, site2, f);
      if (world.rank() == 0)
        cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;

      // left-upper-d_xy <====> right-lower-s
      site1 = i + 1;
      site2 = i + 2 + (2 * Ly);
      mpo_gen.AddTerm(-tsd_xy, bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(-tsd_xy, bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(tsd_xy, bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(tsd_xy, bdna, site1, Fbdnc, site2, f);
      if (world.rank() == 0)
        cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;

    } else if (Ly > 2) {
      //winding term, incomplete
      // left-upper-s <====> right-lower-d_xy
      size_t site1 = i, site2 = i + 3;
      mpo_gen.AddTerm(-tsd_xy, bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(-tsd_xy, bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(tsd_xy, bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(tsd_xy, bdna, site1, Fbdnc, site2, f);
      if (world.rank() == 0)
        cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;
      // left-upper-d_xy <====> right-lower-s
      site1 += 1;
      site2 = site2 - 1;
      mpo_gen.AddTerm(-tsd_xy, bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(-tsd_xy, bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(tsd_xy, bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(tsd_xy, bdna, site1, Fbdnc, site2, f);
      if (world.rank() == 0)
        cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;
    }

    if (y > 0) {
      // left-lower-s <====> right-upper-d_xy
      size_t site1 = i, site2 = i + (2 * Ly) - 1;
      mpo_gen.AddTerm(tsd_xy, bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(tsd_xy, bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(-tsd_xy, bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(-tsd_xy, bdna, site1, Fbdnc, site2, f);
      if (world.rank() == 0)
        cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;

      // left-lower-d_xy <====> right-upper-s
      site1 = i + 1;
      site2 = site2 - 1;
      mpo_gen.AddTerm(tsd_xy, bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(tsd_xy, bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(-tsd_xy, bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(-tsd_xy, bdna, site1, Fbdnc, site2, f);
      if (world.rank() == 0)
        cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;
    } else if (Ly > 2) {
      //winding term
      // left-lower-s <====> right-upper-d_xy
      size_t site1 = i, site2 = i + (4 * Ly) - 1;
      mpo_gen.AddTerm(tsd_xy, bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(tsd_xy, bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(-tsd_xy, bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(-tsd_xy, bdna, site1, Fbdnc, site2, f);
      if (world.rank() == 0)
        cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;

      // left-lower-d_xy <====> right-upper-s
      site1 = i + 1;
      site2 = site2 - 1;
      mpo_gen.AddTerm(tsd_xy, bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(tsd_xy, bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(-tsd_xy, bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(-tsd_xy, bdna, site1, Fbdnc, site2, f);
      if (world.rank() == 0)
        cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;
    }
  }

  //intra band nn hopping, tsd_x^2-y^2
  for (size_t i = 0; i < N; i += 2) {
    size_t y = i % (2 * Ly), x = i / (2 * Ly);
    //horizontal
    if (i < N - (2 * Ly)) {
      size_t site1 = i, site2 = i + (2 * Ly) + 1;
      mpo_gen.AddTerm(-tsd_nn, bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(-tsd_nn, bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(tsd_nn, bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(tsd_nn, bdna, site1, Fbdnc, site2, f);
      if (world.rank() == 0)
        cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;
    }
    //vertical
    if (y < 2 * Ly - 2) {
      size_t site1 = i, site2 = i + 3;
      mpo_gen.AddTerm(tsd_nn, bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(tsd_nn, bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(-tsd_nn, bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(-tsd_nn, bdna, site1, Fbdnc, site2, f);
      if (world.rank() == 0)
        cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;
    } else if (Ly > 2) {
      //winding term
      size_t site1 = i - (2 * Ly) + 3, site2 = i;
      mpo_gen.AddTerm(tsd_nn, bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(tsd_nn, bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(-tsd_nn, bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(-tsd_nn, bdna, site1, Fbdnc, site2, f);
      if (world.rank() == 0)
        cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;
    }
  }

  auto mro = mpo_gen.GenMatReprMPO(false);
  if (world.rank() == 0)
    cout << "MRO generated." << endl;

  // dmrg
  double e0;
  if (!has_bond_dimension_parameter) {
    e0 = qlmps::FiniteDMRG(mps, mro, sweep_params, world);
  } else {
    for (size_t i = 0; i < DMRG_time; i++) {
      size_t D = input_D_set[i];
      if (world.rank() == 0) {
        std::cout << "D_max = " << D << std::endl;
      }
      qlmps::FiniteVMPSSweepParams sweep_params(
          params.Sweeps,
          D, D, params.CutOff,
          qlmps::LanczosParams(params.LanczErr, MaxLanczIterSet[i]),
          params.noise
      );
      e0 = qlmps::FiniteDMRG(mps, mro, sweep_params, world);
    }

  }
  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}
