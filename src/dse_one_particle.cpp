#include "dse.h"
#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/vector.h"
#include <array>
#include <iostream>
#include <string>

using namespace dse;
using namespace std;

momentum *sigmaDiag::NextMom() {
  MomNum += 1;
  ASSERT_ALLWAYS(MomNum < MaxMomNum, "Too many momentum variables! " << MomNum);
  return &(*LoopMom)[MomNum - 1];
}

sigma sigmaDiag::Build(array<momentum, MaxMomNum> &loopMom, int LoopNum,
                       caltype Type) {
  ASSERT_ALLWAYS(LoopNum > 0, "LoopNum must be larger than zero!");
  MomNum = MaxLoopNum;
  LoopMom = &loopMom;

  sigma Sigma;
  //   if (Type == PARQUET)
  Sigma.RexpandBare = false;
  Sigma.LegK = &(*LoopMom)[1];
  Sigma.TauNum = LoopNum;
  Sigma.InTidx = 0;
  Sigma.OutTidx = Sigma.TauNum - 1;

  Sigma.G[0] = gMatrix(Sigma.TauNum, Sigma.InTidx, &(*LoopMom)[2]);
  Sigma.G[1] = gMatrix(Sigma.TauNum, Sigma.InTidx, &(*LoopMom)[3]);

  //   for (int ol = 0; ol < LoopNum; ++ol) {
  //     ver4 LVer =
  //   }
}

// ver4 verDiag::Vertex(array<momentum *, 4> LegK, int InTL, int LoopNum,
//                      int LoopIndex, vector<channel> Channel, int Side,
//                      bool RenormVer4, bool RexpandBare, bool IsFullVer4) {