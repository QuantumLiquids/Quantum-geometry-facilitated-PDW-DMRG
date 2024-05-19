#ifndef QL_HUBBARD_FERMION_QN_ONLY_SPIN_CONSERVATION_H
#define QL_HUBBARD_FERMION_QN_ONLY_SPIN_CONSERVATION_H
#include "qlten/qlten.h"

using qlten::QNCard;
using qlten::U1QNVal;
using qlten::TenIndexDirType;

using TenElemT = qlten::QLTEN_Double;
using QNT = qlten::special_qn::U1QN;
using U1QN = QNT;
using Tensor = qlten::QLTensor<TenElemT, QNT>;

using QNSctT = qlten::QNSector<QNT>;
using IndexT = qlten::Index<QNT>;

const auto qn0 = QNT(
    {QNCard("Sz", U1QNVal(0))}
);

const IndexT pb_out = IndexT({QNSctT(QNT({QNCard("Sz", U1QNVal(1))}), 1), // spin up
                              QNSctT(QNT({QNCard("Sz", U1QNVal(-1))}), 1), // spin down
                              QNSctT(QNT({QNCard("Sz", U1QNVal(0))}), 2)}, // double occupancy, empty
                             TenIndexDirType::OUT
);
const auto pb_in = qlten::InverseIndex(pb_out);

void OperatorInitial();

#endif  //QL_HUBBARD_FERMION_QN_ONLY_SPIN_CONSERVATION_H