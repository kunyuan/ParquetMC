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

int AddToTList(vector<array<int, 4>> &TList, const array<int, 4> T) {
  // find the T array in the list, if failed, create a new array
  for (int i = 0; i < TList.size(); i++) {
    auto t = TList[i];
    ASSERT_ALLWAYS(t[INL] == T[INL],
                   "left Tin must be the same for all subvertex!");
    if (t[OUTL] == T[OUTL] && t[INR] == T[INR] && t[OUTR] == T[OUTR])
      return i;
  }
  TList.push_back(T);
  return TList.size() - 1;
}

momentum *verDiag::NextMom() {
  MomNum += 1;
  ASSERT_ALLWAYS(MomNum < MaxMomNum, "Too many momentum variables! " << MomNum);
  return &(*LoopMom)[MomNum - 1];
}

ver4 verDiag::Build(array<momentum, MaxMomNum> &loopMom, int LoopNum,
                    vector<channel> Channel, caltype Type) {
  ASSERT_ALLWAYS(LoopNum > 0, "LoopNum must be larger than zero!");
  DiagNum = 0;
  MomNum = MaxLoopNum;
  LoopMom = &loopMom;
  array<momentum *, 4> LegK;
  // if (Channel.size() == 1 && Channel[0] == S) {
  //   LegK = {&(*LoopMom)[1], &(*LoopMom)[2], NextMom(), NextMom()};
  // } else {
  //   LegK = {&(*LoopMom)[1], NextMom(), &(*LoopMom)[2], NextMom()};
  // }

  if (Channel.size() == 1) {
    LegK = {&(*LoopMom)[1], NextMom(), &(*LoopMom)[2], NextMom()};
    // LegK = {&(*LoopMom)[1], &(*LoopMom)[1], &(*LoopMom)[2], &(*LoopMom)[2]};
  }

  if (Type == PARQUET)
    return Vertex(LegK, 0, LoopNum, 3, Channel, LEFT, true, false, false);
  else if (Type == BARE)
    return Vertex(LegK, 0, LoopNum, 3, Channel, LEFT, false, false, false);
  else if (Type == RENORMALIZED)
    return Vertex(LegK, 0, LoopNum, 3, Channel, LEFT, true, true, true);
  else
    ABORT("Not implemented!");
}

ver4 verDiag::Vertex(array<momentum *, 4> LegK, int InTL, int LoopNum,
                     int LoopIndex, vector<channel> Channel, int Side,
                     bool RenormVer4, bool RexpandBare, bool IsFullVer4) {
  ver4 Ver4;
  Ver4.ID = DiagNum;
  DiagNum++;
  Ver4.LoopNum = LoopNum;
  Ver4.TauNum = LoopNum + 1;
  Ver4.LegK = LegK;
  Ver4.Side = Side;
  Ver4.InTL = InTL;
  Ver4.Loopidx = LoopIndex;
  Ver4.IsFullVer4 = IsFullVer4;
  Ver4.RenormVer4 = RenormVer4;
  Ver4.RexpandBare = RexpandBare;

  // ASSERT_ALLWAYS(
  //     RexpandBare && RenormVer4,
  //     "RenormVer4 and RexpandBare can not be true at the same time!");

  vector<channel> UST;
  vector<channel> II;
  for (auto &chan : Channel) {
    if (chan == I)
      II.push_back(chan);
    else
      UST.push_back(chan);
  }

  if (LoopNum == 0) {
    // the same for left and right vertex with loopnum=0
    Ver4 = Ver0(Ver4, InTL);
  } else {
    // normal diagram
    Ver4 = ChanI(Ver4, II, InTL, LoopNum, LoopIndex, false);
    Ver4 = ChanUST(Ver4, UST, InTL, LoopNum, LoopIndex, false);

    // counter diagrams if the vertex is on the right
    if (IsFullVer4) {
      if (Ver4.RenormVer4) {
        Ver4 = ChanI(Ver4, II, InTL, LoopNum, LoopIndex, true);
        Ver4 = ChanUST(Ver4, UST, InTL, LoopNum, LoopIndex, true);
      }
    } else {
      if (Ver4.RexpandBare) {
        // counter diagrams if the vertex is on the left
        Ver4 = ChanI(Ver4, {I}, InTL, LoopNum, LoopIndex, true);
        Ver4 = ChanUST(Ver4, {T, U, S}, InTL, LoopNum, LoopIndex, true);
      }
    }
  }

  Ver4.Weight.resize(Ver4.T.size());
  for (auto &d : Ver4.Weight)
    d.SetZero();

  return Ver4;
}

ver4 verDiag::Ver0(ver4 Ver4, int InTL) {
  ////////////// bare interaction ///////////
  // if (Ver4.Side == LEFT)
  //   // Side==left, then make sure INL Tau are the last TauIndex
  //   Ver4.T.push_back({InTL, InTL, InTL, InTL});
  // else
  //   // Side==right, then make sure INR Tau are the last TauIndex
  //   Ver4.T.push_back({InTL + 1, InTL + 1, InTL + 1, InTL + 1});

  // if (Ver4.RexpandBare == true) {
  //   //////////// dressed interaction ///////////
  //   Ver4.T.push_back({InTL, InTL, InTL + 1, InTL + 1});
  //   Ver4.T.push_back({InTL, InTL + 1, InTL + 1, InTL});
  // }

  Ver4.T.push_back({InTL, InTL, InTL, InTL});
  return Ver4;
}

vector<mapT2> CreateMapT2(ver4 &Ver4, ver4 LVer, ver4 RVer, channel Chan,
                          bool IsProjected) {
  ///////////   External and Internal Tau  ////////////////
  vector<mapT2> Map;
  array<int, 2> G0T;
  array<array<int, 2>, 4> GT;
  array<array<int, 4>, 4> LegT;
  int Tidx;

  for (int lt = 0; lt < LVer.T.size(); ++lt)
    for (int rt = 0; rt < RVer.T.size(); ++rt) {

      auto &LvT = LVer.T[lt];
      auto &RvT = RVer.T[rt];

      G0T = {LvT[OUTR], RvT[INL]};
      GT[T] = {RvT[OUTL], LvT[INR]};
      GT[U] = {RvT[OUTL], LvT[INR]};
      GT[S] = {LvT[OUTL], RvT[INR]};

      if (IsProjected == false) {
        LegT[T] = {LvT[INL], LvT[OUTL], RvT[INR], RvT[OUTR]};
        LegT[U] = {LvT[INL], RvT[OUTR], RvT[INR], LvT[OUTL]};
        LegT[S] = {LvT[INL], RvT[OUTL], LvT[INR], RvT[OUTR]};
      } else {
        // LegT[T] = {LvT[INL], LvT[INL], RvT[INR], RvT[INR]};
        // LegT[U] = {LvT[INL], RvT[INR], RvT[INR], LvT[INL]};
        LegT[S] = {LvT[INL], LvT[INL], LvT[INL], LvT[INL]};

        LegT[T] = {LvT[INL], LvT[INL], LvT[INL], LvT[INL]};
        LegT[U] = {LvT[INL], LvT[INL], LvT[INL], LvT[INL]};
      }

      // add T array into the T pool of the vertex
      Tidx = AddToTList(Ver4.T, LegT[Chan]);
      Map.push_back(mapT2{lt, rt, Tidx, G0T, GT[Chan]});
    }
  return Map;
}

ver4 verDiag::ChanUST(ver4 Ver4, vector<channel> Channel, int InTL, int LoopNum,
                      int LoopIndex, bool IsProjected) {
  bubble Bubble;
  Bubble.IsProjected = IsProjected;
  Bubble.InTL = InTL;
  Bubble.Channel = Channel;
  array<double, 4> SymFactor = {0.0, -1.0, 1.0, -0.5};
  bool HasT = bool(std::count(Channel.begin(), Channel.end(), T));
  bool HasU = bool(std::count(Channel.begin(), Channel.end(), U));
  bool HasS = bool(std::count(Channel.begin(), Channel.end(), S));
  if (HasT || HasU)
    Bubble.HasTU = true;
  else
    Bubble.HasTU = false;
  if (HasS)
    Bubble.HasS = true;
  else
    Bubble.HasS = false;

  if (IsProjected) {
    for (auto &s : SymFactor)
      s *= -1;
    // Bubble.LegK[T][INL] = NextMom();
    // Bubble.LegK[T][INR] = NextMom();
    // Bubble.LegK[T][OUTL] = Bubble.LegK[T][INL];
    // Bubble.LegK[T][OUTR] = Bubble.LegK[T][INR];

    // Bubble.LegK[U][INL] = Bubble.LegK[T][INL];
    // Bubble.LegK[U][INR] = Bubble.LegK[T][INR];
    // Bubble.LegK[U][OUTL] = Bubble.LegK[T][INL];
    // Bubble.LegK[U][OUTR] = Bubble.LegK[T][INR];

    // Bubble.LegK[S][INL] = Bubble.LegK[T][INL];
    // Bubble.LegK[S][INR] = Bubble.LegK[T][INR];
    // Bubble.LegK[S][OUTL] = Bubble.LegK[T][INL];
    // Bubble.LegK[S][OUTR] = Bubble.LegK[T][INR];
    for (auto &c : Channel)
      Bubble.LegK[c] = Ver4.LegK;
  } else
    for (auto &c : Channel)
      Bubble.LegK[c] = Ver4.LegK;

  auto &G = Bubble.G;
  auto &LegK = Bubble.LegK;

  G[0] = gMatrix(Ver4.TauNum, InTL, &(*LoopMom)[LoopIndex]);
  for (auto &c : Bubble.Channel)
    G[c] = gMatrix(Ver4.TauNum, InTL, NextMom());

  for (int ol = 0; ol < LoopNum; ol++) {
    // left and right vertex external LegK
    array<momentum *, 4> LLegK[4], RLegK[4];

    ////////////////// T channel ////////////////////////////
    LLegK[T] = {LegK[T][INL], LegK[T][OUTL], G[T].K, G[0].K};
    RLegK[T] = {G[0].K, G[T].K, LegK[T][INR], LegK[T][OUTR]};

    ////////////////// U channel ////////////////////////////
    LLegK[U] = {LegK[U][INL], LegK[U][OUTR], G[U].K, G[0].K};
    RLegK[U] = {G[0].K, G[U].K, LegK[U][INR], LegK[U][OUTL]};

    ////////////////// S channel ////////////////////////////
    LLegK[S] = {LegK[S][INL], G[S].K, LegK[S][INR], G[0].K};
    RLegK[S] = {G[0].K, LegK[S][OUTL], G[S].K, LegK[S][OUTR]};

    for (auto &c : Bubble.Channel) {
      // if (ol == 1 && c == T && LoopNum == 2)
      //   continue;

      // if (IsProjected && c == S)
      //   continue;

      pair Pair;
      Pair.Channel = c;
      Pair.SymFactor = SymFactor[c];
      ////////////////////   Right SubVer  ///////////////////
      int oR = LoopNum - 1 - ol;
      int RInTL = InTL + (ol + 1);
      int Rlopidx = LoopIndex + 1 + ol;

      if (c == U || c == T) {
        Pair.LVer = Vertex(LLegK[c], InTL, ol, LoopIndex + 1, {I, U, S}, LEFT,
                           Ver4.RenormVer4, Ver4.RexpandBare, false);
        Pair.RVer = Vertex(RLegK[c], RInTL, oR, Rlopidx, {I, T, U, S}, RIGHT,
                           Ver4.RenormVer4, Ver4.RenormVer4, true);
      } else if (c == S) {
        Pair.LVer = Vertex(LLegK[c], InTL, ol, LoopIndex + 1, {I, T, U}, LEFT,
                           Ver4.RenormVer4, Ver4.RexpandBare, false);
        Pair.RVer = Vertex(RLegK[c], RInTL, oR, Rlopidx, {I, T, U, S}, RIGHT,
                           Ver4.RenormVer4, Ver4.RenormVer4, true);
      }
      Pair.Map = CreateMapT2(Ver4, Pair.LVer, Pair.RVer, c, IsProjected);
      Bubble.Pair.push_back(Pair);
    }
  }

  Ver4.Bubble.push_back(Bubble);
  return Ver4;
}

vector<mapT4> CreateMapT4(ver4 &Ver4, ver4 LDVer, ver4 LUVer, ver4 RDVer,
                          ver4 RUVer) {
  vector<mapT4> Map;
  array<array<int, 2>, 9> GT; // G Tau pair
  array<int, 4> LegT[4], Tidx;

  for (int ldt = 0; ldt < LDVer.T.size(); ldt++)
    for (int lut = 0; lut < LUVer.T.size(); lut++)
      for (int rdt = 0; rdt < RDVer.T.size(); rdt++)
        for (int rut = 0; rut < RUVer.T.size(); rut++) {
          auto &ldT = LDVer.T[ldt];
          auto &luT = LUVer.T[lut];
          auto &rdT = RDVer.T[rdt];
          auto &ruT = RUVer.T[rut];

          // Tau Index for all possible internal G
          GT[0] = {ldT[OUTR], rdT[INL]};
          GT[1] = {ldT[OUTL], luT[INL]};
          GT[2] = {rdT[OUTL], luT[INR]};
          GT[3] = {ruT[OUTL], ldT[INR]};
          GT[4] = {luT[OUTR], ruT[INL]};
          GT[5] = {rdT[OUTR], ruT[INR]};
          GT[6] = {luT[OUTR], ruT[INL]};
          GT[7] = {ruT[OUTR], rdT[INR]};
          GT[8] = {ruT[OUTR], rdT[INR]};

          // external T for four envelope diagram
          // INL, OUTL, INR, OUTR
          LegT[0] = {ldT[INL], luT[OUTL], rdT[INR], ruT[OUTR]};
          LegT[1] = {ldT[INL], ruT[OUTR], rdT[INR], luT[OUTL]};
          LegT[2] = {ldT[INL], luT[OUTL], ruT[INR], rdT[OUTR]};
          LegT[3] = {ldT[INL], rdT[OUTR], ruT[INR], luT[OUTL]};

          for (int i = 0; i < 4; i++)
            Tidx[i] = AddToTList(Ver4.T, LegT[i]);

          Map.push_back(mapT4{ldt, lut, rdt, rut, Tidx, GT});
        }
  return Map;
}

ver4 verDiag::ChanI(ver4 Ver4, vector<channel> Channel, int InTL, int LoopNum,
                    int LoopIndex, bool IsProjected) {
  if (Channel.size() == 0 || Channel[0] != I)
    return Ver4;

  if (LoopNum != 3)
    return Ver4;

  if (IsProjected)
    return Ver4;

  envelope Env;
  Env.IsProjected = IsProjected;
  Env.InTL = InTL;
  if (IsProjected == false)
    Env.LegK = Ver4.LegK;
  auto &G = Env.G;

  int LDInTL = InTL;
  int LUInTL = InTL + 1;
  int RDInTL = InTL + 2;
  int RUInTL = InTL + 3;

  /////// Initialize G Tau and K Table  /////////
  G[0] = g2Matrix(LDInTL, RDInTL, &(*LoopMom)[LoopIndex]);
  G[1] = g2Matrix(LDInTL, LUInTL, &(*LoopMom)[LoopIndex + 1]);
  G[2] = g2Matrix(RDInTL, LUInTL, &(*LoopMom)[LoopIndex + 2]);
  G[3] = g2Matrix(RUInTL, LDInTL, NextMom());
  G[4] = g2Matrix(LUInTL, RUInTL, NextMom());
  G[5] = g2Matrix(RDInTL, RUInTL, NextMom());
  G[6] = g2Matrix(LUInTL, RUInTL, NextMom());
  G[7] = g2Matrix(RUInTL, RDInTL, NextMom());
  G[8] = g2Matrix(RUInTL, RDInTL, NextMom());

  momentum *InL = Ver4.LegK[INL];
  momentum *OutL = Ver4.LegK[OUTL];
  momentum *InR = Ver4.LegK[INR];
  momentum *OutR = Ver4.LegK[OUTR];

  //////// Initialize all sub-vertex ///////////

  array<momentum *, 4> LegK[10];
  vector<channel> ALL = {I, U, S, T};

  // LD Vertex
  LegK[0] = {InL, G[1].K, G[3].K, G[0].K};
  Env.Ver[0] = Vertex(LegK[0], LDInTL, 0, LoopIndex, ALL, LEFT, Ver4.RenormVer4,
                      Ver4.RenormVer4, true);

  // LU Vertex
  LegK[1] = {G[1].K, OutL, G[2].K, G[4].K};
  LegK[2] = {G[1].K, OutR, G[2].K, G[6].K};
  Env.Ver[1] = Vertex(LegK[1], LUInTL, 0, LoopIndex, ALL, LEFT, Ver4.RenormVer4,
                      Ver4.RenormVer4, true);
  Env.Ver[2] = Vertex(LegK[2], LUInTL, 0, LoopIndex, ALL, LEFT, Ver4.RenormVer4,
                      Ver4.RenormVer4, true);

  // RD Vertex
  LegK[3] = {G[0].K, G[2].K, InR, G[5].K};
  LegK[4] = {G[0].K, G[2].K, G[7].K, OutR};
  LegK[5] = {G[0].K, G[2].K, G[8].K, OutL};
  for (int i = 3; i <= 5; i++)
    Env.Ver[i] = Vertex(LegK[i], RDInTL, 0, LoopIndex, ALL, RIGHT,
                        Ver4.RenormVer4, Ver4.RenormVer4, true);

  // RU Vertex
  LegK[6] = {G[4].K, G[3].K, G[5].K, OutR};
  LegK[7] = {G[6].K, G[3].K, G[5].K, OutL};
  LegK[8] = {G[4].K, G[3].K, InR, G[7].K};
  LegK[9] = {G[6].K, G[3].K, InR, G[8].K};
  for (int i = 6; i <= 9; i++)
    Env.Ver[i] = Vertex(LegK[i], RUInTL, 0, LoopIndex, ALL, RIGHT,
                        Ver4.RenormVer4, Ver4.RenormVer4, true);

  //////// T map (for all four envelope diagram) //////
  // four diagrams have the same sub-vertex Tau configuration
  // so here we just use the first diagram
  Env.Map = CreateMapT4(Ver4, Env.Ver[0], Env.Ver[1], Env.Ver[3], Env.Ver[6]);

  Env.SymFactor = {-1.0, 1.0, -1.0, 1.0};
  Ver4.Envelope.push_back(Env);
  return Ver4;
}

string verDiag::ToString(const ver4 &Ver4, string indent, int Level) {
  string SideStr;
  if (Ver4.Side == LEFT)
    SideStr = "L";
  else
    SideStr = "R";
  // string Info =
  //     indent +
  //     "==============================================================\n";
  string Info =
      indent +
      fmt::format("├Level: {0}, {1}Ver4 ID: {2}, LoopNum: {3}, "
                  "RexpandBare: {4}, RenormVer4: {5}, IsFullVer4: {6}\n",
                  Level, SideStr, Ver4.ID, Ver4.LoopNum, Ver4.RexpandBare,
                  Ver4.RenormVer4, Ver4.IsFullVer4);
  Info += indent + "├─T: ";
  for (auto &t : Ver4.T)
    Info +=
        fmt::format("({0}, {1}, {2}, {3}), ", t[INL], t[OUTL], t[INR], t[OUTR]);

  Info += "\n";
  if (Ver4.Bubble.size() > 0)
    Info += indent + "├─LegK: ";
  else
    Info += indent + "└─LegK: ";

  momentum *start = &(*LoopMom)[0];
  // momentum *InL = Ver4.LegK[INL];

  Info += fmt::format("{0}, {1}, {2}, {3}", Ver4.LegK[INL] - start,
                      Ver4.LegK[OUTL] - start, Ver4.LegK[INR] - start,
                      Ver4.LegK[OUTR] - start);

  // Info += "\n" + indent;
  // Info += "---------------------------------------------------------\n";
  Info += "\n";

  for (auto &bubble : Ver4.Bubble) {
    Info +=
        indent + fmt::format("└─BUBBLE - Projected: {0}\n", bubble.IsProjected);
    // Info += indent + ". │\n";
    for (int p = 0; p < bubble.Pair.size(); p++) {
      pair pp = bubble.Pair[p];
      Info += indent + fmt::format(". ├PAIR - Channel: {0}, LVerLoopNum: {1}\n",
                                   pp.Channel, pp.LVer.LoopNum);

      // Info += indent + " . │\n";
      Info += indent + fmt::format(". ├─Map:");
      for (auto &m : pp.Map)
        Info += fmt::format("({0},{1}):{2}, ", m.LVerTidx, m.RVerTidx, m.Tidx);

      Info += "\n" + indent + ". │\n";
      Info += ToString(pp.LVer, indent + ". ", Level + 1);
      // Info +=
      //     indent +
      //     ".....................................................\n";

      Info += indent + ". │\n";
      Info += ToString(pp.RVer, indent + ". ", Level + 1);
      Info += indent + ".  \n";
      // Info += "\n";
      // Info += indent + "\n";
      // Info += "\n" + indent + " . │\n";
      // Info += "\n" + indent;
      // Info += ".....................................................\n";

      // Info += fmt::format("  G1 Internal T Map: ");
      // for (auto &m : pp.Map)
      //   Info += fmt::format("({0}, {1}): {2}-{3},
      //   ", m.LVerTidx, m.RVerTidx, m.GT[0],
      //                       m.G1T[1]);
      // Info += "\n";

      // Info += fmt::format("  G2 Internal T Map: ");
      // for (auto &m : pp.Map)
      //   Info += fmt::format("({0}, {1}): {2}-{3},
      //   ", m.LVerT, m.RVerT, m.G2T[0],
      //                       m.G2T[1]);
      // Info += indent + "__\n";

      // Info += indent + "--------------------------------------"
      //                  "----------------------\n";
    }
  }

  // Info += "=======================================================\n";
  // Info += indent +
  // "----------------------------------------------------\n";
  return Info;
}

bool verTest() {
  verDiag VerDiag;
  vector<channel> Chan = {T, U, S};
  //   ver4 *Root = VerDiag.Build(1, Chan, NORMAL);
  //   cout << VerDiag.ToString(*Root) << endl;
  return true;
}
