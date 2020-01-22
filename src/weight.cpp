#include "weight.h"
#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/vector.h"
#include <array>
#include <iostream>
#include <string>

using namespace diag;
using namespace std;

void weight::ReadDiagrams() {
  int ID = 0;
  for (auto &name : Para.GroupName) {
    // construct filename based on format string and group id
    string FileName = fmt::format(Para.DiagFileFormat, name);
    ifstream DiagFile(FileName);
    ASSERT_ALLWAYS(DiagFile.is_open(),
                   "Unable to find the file " << FileName << endl);
    // group Group;
    LOG_INFO("Find " << FileName << "\n");
    // vector<green> GList;
    istream &DiagFileStream = DiagFile;
    Groups.push_back(ReadOneGroup(DiagFileStream));
    Groups.back().Name = name;
    Groups.back().ID = ID;
    ID++;
  }

  // cout << "After read" << endl;
  // cout << ToString(*(GroupList[0].DiagList[0].GIndex[0])) << endl;
  Initialization();
}

void weight::Initialization() {

  LOG_INFO("Initializating diagram states ...")
  for (auto &group : Groups) {
    group.ReWeight = 1.0;
  }

  // vector<dse::channel> Chan = {dse::T, dse::U, dse::S};
  dse::channel Chan[4] = {dse::I, dse::T, dse::U, dse::S};
  for (int c = 0; c < 4; c++)
    for (int order = 1; order <= Para.Order; order++) {
      vector<dse::channel> chan = {Chan[c]};
      Ver4Root[order][c] =
          VerDiag.Build(Var.LoopMom, order, chan, dse::caltype::PARQUET);
      LOG_INFO(VerDiag.ToString(Ver4Root[order][c]));
    }

  LOG_INFO("Initializating MC variables ...")
  // initialize momentum variables
  for (auto &mom : Var.LoopMom)
    for (int i = 0; i < D; i++)
      mom[i] = Random.urn() * Para.Kf / sqrt(D);

  for (int i = 0; i < MaxTauNum; i++) {
    Var.Tau[i] = Random.urn() * Para.Beta;
  }

  // initialize spin variables
  // for (auto &sp : Var.LoopSpin)
  //   sp = (spin)(Random.irn(0, 1));

  Var.CurrExtMomBin = 0;

  // Var.LoopMom[0].fill(0.0);
  Var.LoopMom[0] = Para.ExtMomTable[Var.CurrExtMomBin];

  for (int i = 1; i < D; i++) {
    Var.LoopMom[1][i] = 0.0;
    Var.LoopMom[2][i] = 0.0;
  }
  Var.LoopMom[1][0] = Para.Kf;
  Var.LoopMom[2][0] = Para.Kf;

  // initialize external tau
  // Var.Tau[0] = 0.0;
  // Var.Tau[1] = 1.0e-10; // do not make Tau[1]==Tau[0], otherwise the Green's
  // function is not well-defined

  Var.CurrTau = Var.Tau[1] - Var.Tau[0];

  // initialize group

  Var.CurrVersion = 0;

  Var.CurrGroup = &Groups[0];

  Var.CurrIRScaleBin = ScaleBinSize / 1.5;

  Var.CurrChannel = dse::T;

  // initialize RG staff
  // Var.CurrScale = ScaleBinSize - 1;
  Var.CurrScale = Para.Kf;

  LOG_INFO("Calculating the weights of all objects...")

  // ChangeGroup(*Var.CurrGroup, true);
  GetNewWeight(*Var.CurrGroup);
  AcceptChange(*Var.CurrGroup);

  LOG_INFO("Initializating variables done.")
}

double weight::GetNewWeight(group &Group) {
  Group.NewWeight = Evaluate(Group.Order, Var.CurrChannel);
  return Group.NewWeight.Sum();
}

void weight::AcceptChange(group &Group) {
  Var.CurrVersion++;
  Var.CurrGroup = &Group;
  Group.Weight = Group.NewWeight.Sum(); // accept group  newweight
}

void weight::RejectChange(group &Group) { return; }

void weight::Measure(double WeightFactor) {
  if (Para.Type == RG && Para.Vertex4Type == MOM_ANGLE) {
    // if (Var.CurrScale >= Para.ScaleTable[Var.CurrIRScaleBin])
    VerQTheta.Measure(Var.LoopMom[1], Var.LoopMom[2], Var.CurrExtMomBin,
                      Var.CurrGroup->Order,
                      Var.Tau[Var.CurrGroup->TauNum - 1] - Var.Tau[0],
                      Var.CurrChannel, Var.CurrGroup->Weight, WeightFactor);
  }
}

void weight::Save(bool Simple) {
  if (Para.Type == RG && Para.Vertex4Type == MOM_ANGLE) {
    VerQTheta.Save(Simple);
  }
}

void weight::ClearStatis() {
  if (Para.Type == RG && Para.Vertex4Type == MOM_ANGLE) {
    VerQTheta.ClearStatis();
  }
}

void weight::LoadWeight() { VerQTheta.LoadWeight(); }