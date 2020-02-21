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

int AddToGList(vector<green> &GList, const array<int, 2> T) {
  // find the T array in the list, if failed, create a new array
  for (int i = 0; i < GList.size(); i++) {
    auto t = GList[i].T;
    if (t[OUT] == T[OUT] && t[IN] == T[IN])
      return i;
  }
  GList.push_back(green({T, 0.0}));
  return GList.size() - 1;
}

ver4 verDiag::Vertex(int Level, int LoopNum, int LoopIndex, int InTL,
                     const vector<channel> &Channel,
                     const array<momentum *, 4> &LegK, int Side, bool InBox) {
  ver4 Ver4;
  Ver4.Level = Level;
  Ver4.LoopNum = LoopNum;
  Ver4.TauNum = LoopNum + 1;
  Ver4.Side = Side;
  Ver4.InTL = InTL;
  Ver4.Loopidx = LoopIndex;
  Ver4.InBox = InBox;
  Ver4.LegK = LegK;

  vector<channel> UST;
  vector<channel> II;
  for (auto &chan : Channel) {
    // if one wants the bare diagrams, filter all counter diagrams!
    if (Para.Type == BARE && chan >= 4)
      continue;

    if (chan == I)
      II.push_back(chan);
    else
      UST.push_back(chan);
  }

  // check if there is any counter-term or not
  Ver4.HasCT = false;
  if (Para.Type != BARE) {
    for (auto &chan : Channel)
      if (chan >= 4)
        Ver4.HasCT = true;
  }

  if (LoopNum == 0) {
    // the same for left and right vertex with loopnum=0
    Ver4 = Ver0(Ver4);
  } else {
    // normal diagram
    Ver4 = ChanI(Ver4, II);
    Ver4 = ChanUST(Ver4, UST);
  }

  Ver4.Weight.resize(Ver4.T.size());
  return Ver4;
}

ver4 verDiag::Ver0(ver4 Ver4) {
  ////////////// bare interaction ///////////
  int InTL = Ver4.InTL;
  Ver4.T.push_back({InTL, InTL, InTL, InTL});
  return Ver4;
}

vector<indexMap> CreateIndexMap(ver4 &Ver4, const ver4 &LVer, const ver4 &RVer,
                                channel Chan) {
  ///////////   External and Internal Tau  ////////////////
  vector<indexMap> Map;
  int GT0, GT, Tidx;
  array<int, 2> GTpair;
  array<int, 4> LegT;
  ASSERT_ALLWAYS(Chan != I, "CreateIndexMap is not for I channel!");

  for (int lt = 0; lt < LVer.T.size(); ++lt)
    for (int rt = 0; rt < RVer.T.size(); ++rt) {

      auto &LvT = LVer.T[lt];
      auto &RvT = RVer.T[rt];

      GT0 = AddToGList(Ver4.G[0], {LvT[OUTR], RvT[INL]});

      if (Chan == T || Chan == TC || Chan == U || Chan == UC)
        GTpair = {RvT[OUTL], LvT[INR]};
      else if (Chan == S)
        GTpair = {LvT[OUTL], RvT[INR]};
      else
        ABORT("The channel does not exist!");

      GT = AddToGList(Ver4.G[Chan], GTpair);

      switch (Chan) {
      case T:
        LegT = {LvT[INL], LvT[OUTL], RvT[INR], RvT[OUTR]};
        break;
      case U:
        LegT = {LvT[INL], RvT[OUTR], RvT[INR], LvT[OUTL]};
        break;
      case S:
        LegT = {LvT[INL], RvT[OUTL], LvT[INR], RvT[OUTR]};
        break;
      case TC:
      case UC:
        // counterterms are equal-time
        LegT = {LvT[INL], LvT[INL], LvT[INL], LvT[INL]};
      default:
        ABORT("The channel does not exist!");
        break;
      }

      // add T array into the T pool of the vertex
      Tidx = AddToTList(Ver4.T, LegT);
      Map.push_back(indexMap{lt, rt, Tidx, GT0, GT});
    }
  return Map;
}

ver4 verDiag::ChanUST(ver4 Ver4, vector<channel> Channel) {
  int InTL = Ver4.InTL;
  auto LegK = Ver4.LegK;
  auto K = Ver4.K;

  vector<channel> FULL = {I, T, U, S, TC, UC};
  vector<channel> F = {I, U, S, TC, UC};
  vector<channel> V = {I, T, U, TC, UC};
  vector<channel> FULL_CT = {I, T, TC};
  vector<channel> F_CT = {I, TC};

  array<momentum *, 4> LLegK[4], RLegK[4];

  ////////////////// T channel ////////////////////////////
  LLegK[T] = {LegK[INL], LegK[OUTL], &K[T], &K[0]};
  RLegK[T] = {&K[0], &K[T], LegK[INR], LegK[OUTR]};

  ////////////////// U channel ////////////////////////////
  LLegK[U] = {LegK[INL], LegK[OUTR], &K[U], &K[0]};
  RLegK[U] = {&K[0], &K[U], LegK[INR], LegK[OUTL]};

  ////////////////// S channel ////////////////////////////
  LLegK[S] = {LegK[INL], &K[S], LegK[INR], &K[0]};
  RLegK[S] = {&K[0], LegK[OUTL], &K[S], LegK[OUTR]};

  for (int ol = 0; ol < Ver4.LoopNum; ol++) {
    for (auto &chan : Channel) {
      pair Pair;
      ////////////////////   Right SubVer  ///////////////////
      int Level = Ver4.Level + 1;
      int oR = Ver4.LoopNum - 1 - ol;
      int RInTL = Ver4.InTL + (ol + 1);
      int Llopidx = Ver4.Loopidx + 1;
      int Rlopidx = Ver4.Loopidx + 1 + ol;

      switch (chan) {
      case T:
      case U:
        Pair.LVer =
            Vertex(Level, ol, Llopidx, InTL, F, LLegK[chan], LEFT, Ver4.InBox);
        Pair.RVer = Vertex(Level, oR, Rlopidx, RInTL, FULL, RLegK[chan], RIGHT,
                           Ver4.InBox);
        break;
      case S:
        Pair.LVer =
            Vertex(Level, ol, Llopidx, InTL, V, LLegK[chan], LEFT, Ver4.InBox);
        Pair.RVer = Vertex(Level, oR, Rlopidx, RInTL, FULL, RLegK[chan], RIGHT,
                           Ver4.InBox);
        break;
      case TC:
      case UC:
        Pair.LVer =
            Vertex(Level, ol, Llopidx, InTL, F_CT, LLegK[chan], LEFT, true);
        Pair.RVer = Vertex(Level, oR, Rlopidx, RInTL, FULL_CT, RLegK[chan],
                           RIGHT, true);
        break;
      default:
        ABORT("The channel does not exist!");
        break;
      }

      Pair.Channel = chan;
      Pair.Map = CreateIndexMap(Ver4, Pair.LVer, Pair.RVer, chan);
      Ver4.Pair.push_back(Pair);
    }
  }
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

ver4 verDiag::ChanI(ver4 Ver4, vector<channel> Channel) {
  if (Channel.size() == 0)
    return Ver4;

  if (Ver4.LoopNum != 3)
    return Ver4;

  // envelope Env;
  // Env.IsProjected = IsProjected;
  // Env.InTL = InTL;
  // if (IsProjected == false)
  //   Env.LegK = Ver4.LegK;
  // auto &G = Env.G;

  // int LDInTL = InTL;
  // int LUInTL = InTL + 1;
  // int RDInTL = InTL + 2;
  // int RUInTL = InTL + 3;

  // /////// Initialize G Tau and K Table  /////////
  // G[0] = g2Matrix(LDInTL, RDInTL, &(*LoopMom)[LoopIndex]);
  // G[1] = g2Matrix(LDInTL, LUInTL, &(*LoopMom)[LoopIndex + 1]);
  // G[2] = g2Matrix(RDInTL, LUInTL, &(*LoopMom)[LoopIndex + 2]);
  // G[3] = g2Matrix(RUInTL, LDInTL, NextMom());
  // G[4] = g2Matrix(LUInTL, RUInTL, NextMom());
  // G[5] = g2Matrix(RDInTL, RUInTL, NextMom());
  // G[6] = g2Matrix(LUInTL, RUInTL, NextMom());
  // G[7] = g2Matrix(RUInTL, RDInTL, NextMom());
  // G[8] = g2Matrix(RUInTL, RDInTL, NextMom());

  // momentum *InL = Ver4.LegK[INL];
  // momentum *OutL = Ver4.LegK[OUTL];
  // momentum *InR = Ver4.LegK[INR];
  // momentum *OutR = Ver4.LegK[OUTR];

  // //////// Initialize all sub-vertex ///////////

  // array<momentum *, 4> LegK[10];
  // vector<channel> ALL = {I, U, S, T};

  // bool HasBeenBoxed = Ver4.HasBeenBoxed || IsProjected;

  // // LD Vertex
  // LegK[0] = {InL, G[1].K, G[3].K, G[0].K};
  // Env.Ver[0] = Vertex(LegK[0], LDInTL, 0, LoopIndex, ALL, LEFT,
  // Ver4.RenormVer4,
  //                     Ver4.RenormVer4, true, HasBeenBoxed);

  // // LU Vertex
  // LegK[1] = {G[1].K, OutL, G[2].K, G[4].K};
  // LegK[2] = {G[1].K, OutR, G[2].K, G[6].K};
  // Env.Ver[1] = Vertex(LegK[1], LUInTL, 0, LoopIndex, ALL, LEFT,
  // Ver4.RenormVer4,
  //                     Ver4.RenormVer4, true, Ver4.HasBeenBoxed);
  // Env.Ver[2] = Vertex(LegK[2], LUInTL, 0, LoopIndex, ALL, LEFT,
  // Ver4.RenormVer4,
  //                     Ver4.RenormVer4, true, HasBeenBoxed);

  // // RD Vertex
  // LegK[3] = {G[0].K, G[2].K, InR, G[5].K};
  // LegK[4] = {G[0].K, G[2].K, G[7].K, OutR};
  // LegK[5] = {G[0].K, G[2].K, G[8].K, OutL};
  // for (int i = 3; i <= 5; i++)
  //   Env.Ver[i] = Vertex(LegK[i], RDInTL, 0, LoopIndex, ALL, RIGHT,
  //                       Ver4.RenormVer4, Ver4.RenormVer4, true,
  //                       HasBeenBoxed);

  // // RU Vertex
  // LegK[6] = {G[4].K, G[3].K, G[5].K, OutR};
  // LegK[7] = {G[6].K, G[3].K, G[5].K, OutL};
  // LegK[8] = {G[4].K, G[3].K, InR, G[7].K};
  // LegK[9] = {G[6].K, G[3].K, InR, G[8].K};
  // for (int i = 6; i <= 9; i++)
  //   Env.Ver[i] = Vertex(LegK[i], RUInTL, 0, LoopIndex, ALL, RIGHT,
  //                       Ver4.RenormVer4, Ver4.RenormVer4, true,
  //                       HasBeenBoxed);

  // //////// T map (for all four envelope diagram) //////
  // // four diagrams have the same sub-vertex Tau configuration
  // // so here we just use the first diagram
  // Env.Map = CreateMapT4(Ver4, Env.Ver[0], Env.Ver[1], Env.Ver[3],
  // Env.Ver[6]);

  // Env.SymFactor = {-1.0, 1.0, -1.0, 1.0};
  // Ver4.Envelope.push_back(Env);
  return Ver4;
}

string verDiag::ToString(const ver4 &Ver4, string indent) {
  string SideStr;
  if (Ver4.Side == LEFT)
    SideStr = "LEFT";
  else
    SideStr = "RIGHT";
  // string Info =
  //     indent +
  //     "==============================================================\n";
  string Info =
      indent + fmt::format("├Level: {0}, {1}, LoopNum: {2}, Boxed: {3}\n",
                           Ver4.Level, SideStr, Ver4.LoopNum, Ver4.InBox);
  Info += indent + "├─T: ";
  for (auto &t : Ver4.T)
    Info +=
        fmt::format("({0}, {1}, {2}, {3}), ", t[INL], t[OUTL], t[INR], t[OUTR]);

  Info += "\n";
  Info += indent + fmt::format("└─\n");
  // Info += indent + ". │\n";
  for (int p = 0; p < Ver4.Pair.size(); p++) {
    pair pp = Ver4.Pair[p];
    Info += indent + fmt::format(". ├PAIR - Channel: {0}, LVerLoopNum: {1}\n",
                                 pp.Channel, pp.LVer.LoopNum);

    // Info += indent + " . │\n";
    Info += indent + fmt::format(". ├─Map:");
    for (auto &m : pp.Map)
      Info += fmt::format("({0},{1}):{2}, ", m.LVerTidx, m.RVerTidx, m.Tidx);

    Info += "\n" + indent + ". │\n";
    Info += ToString(pp.LVer, indent + ". ");
    // Info +=
    //     indent +
    //     ".....................................................\n";

    Info += indent + ". │\n";
    Info += ToString(pp.RVer, indent + ". ");
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
