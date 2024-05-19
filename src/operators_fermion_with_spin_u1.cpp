#include "./ql_hubbard_fermion_qn_only_spin_conservation.h"

Tensor sz = Tensor({pb_in, pb_out});
Tensor sp = Tensor({pb_in, pb_out});
Tensor sm = Tensor({pb_in, pb_out});
Tensor id = Tensor({pb_in, pb_out});
Tensor sx = Tensor({pb_in, pb_out});
Tensor isy =
    Tensor({pb_in, pb_out}); // 1i*sy to write as a real value, note to change the coefficient,i.e. sy*sy=- isy*isy
//auto sy = Tensor({pb_in, pb_out});

Tensor f = Tensor({pb_in, pb_out}); //fermion's insertion operator

auto bupc = Tensor({pb_in, pb_out}); //hardcore boson, b_up^creation, used for JW transformation
auto bupa = Tensor({pb_in, pb_out}); //hardcore boson, b_up^annihilation
auto bdnc = Tensor({pb_in, pb_out}); //hardcore boson, b_down^creation
auto bdna = Tensor({pb_in, pb_out}); //hardcore boson, b_down^annihilation


auto bupcF = Tensor({pb_in, pb_out}); // matrix product of bupc * f
auto bupaF = Tensor({pb_in, pb_out});
auto Fbdnc = Tensor({pb_in, pb_out});
auto Fbdna = Tensor({pb_in, pb_out});

auto cupccdnc = Tensor({pb_in, pb_out}); // c_up^creation * c_down^creation=b_up^creation*b_down^creation*F

auto cdnacupa = Tensor({pb_in, pb_out}); // onsite pair, usually c_down*c_up


auto Uterm = Tensor({pb_in, pb_out}); // Hubbard Uterm, nup*ndown

auto nf = Tensor({pb_in, pb_out}); // nup+ndown, fermion number

auto nfsquare = Tensor({pb_in, pb_out}); // nf^2
auto nup = Tensor({pb_in, pb_out}); // fermion number of spin up
auto ndn = Tensor({pb_in, pb_out}); // ndown


void OperatorInitial() {
  static bool initialized = false;
  if (!initialized) {
    sz({0, 0}) = 0.5;
    sz({1, 1}) = -0.5;
    sp({0, 1}) = 1.0;
    sm({1, 0}) = 1.0;
    id({0, 0}) = 1;
    id({1, 1}) = 1;
    id({2, 2}) = 1;
    id({3, 3}) = 1;
//    sx({0, 1}) = 0.5;
//    sx({1, 0}) = 0.5;
//    isy({0, 1}) = 0.5;
//    isy({1, 0}) = -0.5;
    //sy({1,2}) = -0.5i;
    //sy({2,1}) = 0.5i;

    f({0, 0}) = -1;
    f({1, 1}) = -1;
    f({2, 2}) = 1;
    f({3, 3}) = 1;

    bupc({2, 1}) = 1;
    bupc({0, 3}) = 1;
    bdnc({2, 0}) = 1;
    bdnc({1, 3}) = 1;
    bupa({1, 2}) = 1;
    bupa({3, 0}) = 1;
    bdna({0, 2}) = 1;
    bdna({3, 1}) = 1;

    bupcF({2, 1}) = -1;
    bupcF({0, 3}) = 1;
    Fbdnc({2, 0}) = 1;
    Fbdnc({1, 3}) = -1;
    bupaF({1, 2}) = 1;
    bupaF({3, 0}) = -1;
    Fbdna({0, 2}) = -1;
    Fbdna({3, 1}) = 1;

    cupccdnc({2, 3}) = 1;
    cdnacupa({3, 2}) = 1;

    Uterm({2, 2}) = 1;

    nf({0, 0}) = 1;
    nf({1, 1}) = 1;
    nf({2, 2}) = 2;

    nfsquare({0, 0}) = 1;
    nfsquare({1, 1}) = 1;
    nfsquare({2, 2}) = 4;

    nup({2, 2}) = 1;
    nup({0, 0}) = 1;
    ndn({2, 2}) = 1;
    ndn({1, 1}) = 1;

    initialized = true;
  }
}
