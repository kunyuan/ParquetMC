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

int AddToTList(vector<array<int, 4>> &TList, const array<int, 4> &T) {
  // find the T array in the list, if failed, create a new array
  // cout << "AddTtoList" << endl;
  for (int i = 0; i < TList.size(); i++) {
    auto t = TList[i];
    // cout << "List: " << t[0] << ", " << t[1] << ", " << t[2] << ", " << t[3]
    //      << endl;

    ASSERT_ALLWAYS(t[INL] == T[INL],
                   "left Tin must be the same for all subvertex!"
                       << t[INL] << " vs " << T[INL]);
    if (t[OUTL] == T[OUTL] && t[INR] == T[INR] && t[OUTR] == T[OUTR])
      return i;
  }
  // cout << "Added: " << T[0] << ", " << T[1] << ", " << T[2] << ", " << T[3]
  //      << endl;
  TList.push_back(T);
  return TList.size() - 1;
}

int AddToGList(vector<green> &GList, const array<int, 2> &T) {
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
  Ver4.Channel = Channel;

  if (LoopNum == 0) {
    // the same for left and right vertex with loopnum=0
    Ver0(Ver4);
  } else {
    vector<channel> UST;
    vector<channel> II;
    for (auto &chan : Channel) {
      if (Para.Type == BARE && chan >= 4)
        // if one wants bare diagrams, filter all counter diagrams!
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
    // normal diagram
    ChanI(Ver4, II);
    ChanUST(Ver4, UST);
  }

  Ver4.Weight.resize(Ver4.T.size());
  return Ver4;
}

void verDiag::Ver0(ver4 &Ver4) {
  ////////////// bare interaction ///////////
  int InTL = Ver4.InTL;

  // cout << "Ver0: " << InTL << ", Order: " << Ver4.LoopNum
  //      << ", Level: " << Ver4.Level << ", Tsize: " << Ver4.T.size() << endl;

  auto &LegK = Ver4.LegK;

  // cout << "ver0LegK0: " << LegK[INL] << endl;
  // cout << "ver0LegK1: " << LegK[OUTL] << endl;
  // cout << "ver0LegK2: " << LegK[INR] << endl;
  // cout << "ver0LegK3: " << LegK[OUTR] << endl;

  AddToTList(Ver4.T, {InTL, InTL, InTL, InTL});
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

      switch (Chan) {
      case T:
        LegT = {LvT[INL], LvT[OUTL], RvT[INR], RvT[OUTR]};
        GTpair = {RvT[OUTL], LvT[INR]};
        break;
      case U:
        LegT = {LvT[INL], RvT[OUTR], RvT[INR], LvT[OUTL]};
        GTpair = {RvT[OUTL], LvT[INR]};
        break;
      case S:
        LegT = {LvT[INL], RvT[OUTL], LvT[INR], RvT[OUTR]};
        GTpair = {LvT[OUTL], RvT[INR]};
        break;
      case TC:
      case UC:
        // counterterms are equal-time
        LegT = {LvT[INL], LvT[INL], LvT[INL], LvT[INL]};
        GTpair = {RvT[OUTL], LvT[INR]};
        break;
      default:
        ABORT("The channel does not exist! " << Chan);
        break;
      }

      GT = AddToGList(Ver4.G[Chan], GTpair);

      // add T array into the T pool of the vertex
      // cout << Ver4.InTL << "; " << LegT[0] << ", " << LegT[1] << ", " <<
      // LegT[2]
      //      << ", " << LegT[3] << endl;
      // int tail = Ver4.T.size();
      // if (tail > 0)
      //   cout << Ver4.InTL << "; " << Ver4.T[tail - 1][0] << ", " << endl;
      Tidx = AddToTList(Ver4.T, LegT);
      Map.push_back(indexMap{lt, rt, Tidx, GT0, GT});
    }
  return Map;
}

void verDiag::ChanUST(ver4 &Ver4, const vector<channel> &Channel) {
  int InTL = Ver4.InTL;
  auto &LegK = Ver4.LegK;
  auto &K = Ver4.K;

  vector<channel> FULL = {I, T, U, S, TC, UC};
  vector<channel> F = {I, U, S, TC, UC};
  vector<channel> V = {I, T, U, TC, UC};
  vector<channel> FULL_CT = {I, T, TC};
  vector<channel> F_CT = {I, TC};

  array<momentum *, 4> LLegK, RLegK;

  for (auto &chan : Channel) {
    switch (chan) {
    case T:
    case TC:
      LLegK = {LegK[INL], LegK[OUTL], &K[chan], &K[0]};
      RLegK = {&K[0], &K[chan], LegK[INR], LegK[OUTR]};
      break;
    case U:
    case UC:
      LLegK = {LegK[INL], LegK[OUTR], &K[chan], &K[0]};
      RLegK = {&K[0], &K[chan], LegK[INR], LegK[OUTL]};
      break;
    case S:
      LLegK = {LegK[INL], &K[chan], LegK[INR], &K[0]};
      RLegK = {&K[0], LegK[OUTL], &K[chan], LegK[OUTR]};
      break;
    default:
      ABORT("Channel does not exist! " << chan);
      break;
    }
    // cout << "LegK0: " << LegK[INL] << endl;
    // cout << "LegK1: " << LegK[OUTL] << endl;
    // cout << "LegK2: " << LegK[INR] << endl;
    // cout << "LegK3: " << LegK[OUTR] << endl;

    for (int ol = 0; ol < Ver4.LoopNum; ol++) {
      bubble bub;
      ////////////////////   Right SubVer  ///////////////////
      int Level = Ver4.Level + 1;
      int oR = Ver4.LoopNum - 1 - ol;
      int RInTL = Ver4.InTL + (ol + 1);
      int Llopidx = Ver4.Loopidx + 1;
      int Rlopidx = Ver4.Loopidx + 1 + ol;

      switch (chan) {
      case T:
      case U:
        bub.LVer = Vertex(Level, ol, Llopidx, InTL, F, LLegK, LEFT, Ver4.InBox);
        bub.RVer =
            Vertex(Level, oR, Rlopidx, RInTL, FULL, RLegK, RIGHT, Ver4.InBox);
        break;
      case S:
        bub.LVer = Vertex(Level, ol, Llopidx, InTL, V, LLegK, LEFT, Ver4.InBox);
        bub.RVer =
            Vertex(Level, oR, Rlopidx, RInTL, FULL, RLegK, RIGHT, Ver4.InBox);
        break;
      case TC:
      case UC:
        bub.LVer = Vertex(Level, ol, Llopidx, InTL, F_CT, LLegK, LEFT, true);
        bub.RVer =
            Vertex(Level, oR, Rlopidx, RInTL, FULL_CT, RLegK, RIGHT, true);
        break;
      default:
        ABORT("The channel does not exist!");
        break;
      }

      bub.Channel = chan;
      bub.Map = CreateIndexMap(Ver4, bub.LVer, bub.RVer, chan);
      Ver4.Bubble.push_back(bub);
    }
  }
}

void verDiag::ResetMomMap(ver4 &Ver4, const array<momentum *, 4> &LegK) {
  Ver4.LegK = LegK;
  auto &K = Ver4.K;
  array<momentum *, 4> LLegK, RLegK;

  for (auto &bub : Ver4.Bubble) {
    auto chan = bub.Channel;
    switch (chan) {
    case T:
    case TC:
      LLegK = {LegK[INL], LegK[OUTL], &K[chan], &K[0]};
      RLegK = {&K[0], &K[chan], LegK[INR], LegK[OUTR]};
      break;
    case U:
    case UC:
      LLegK = {LegK[INL], LegK[OUTR], &K[chan], &K[0]};
      RLegK = {&K[0], &K[chan], LegK[INR], LegK[OUTL]};
      break;
    case S:
      LLegK = {LegK[INL], &K[chan], LegK[INR], &K[0]};
      RLegK = {&K[0], LegK[OUTL], &K[chan], LegK[OUTR]};
      break;
    default:
      ABORT("Channel does not exist! " << chan);
      break;
    }
    ResetMomMap(bub.LVer, LLegK);
    ResetMomMap(bub.RVer, RLegK);
  }
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

void verDiag::ChanI(ver4 &Ver4, const vector<channel> &Channel) {
  if (Channel.size() == 0)
    return;

  if (Ver4.LoopNum != 3)
    return;

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
      indent +
      fmt::format("├Level: {0}, {1}, LoopNum: {2}, Boxed: {3}, G0:{4}\n",
                  Ver4.Level, SideStr, Ver4.LoopNum, Ver4.InBox,
                  Ver4.G[0].size());
  Info += indent + fmt::format("├─Tau[{0}]: ", Ver4.T.size());
  for (auto &t : Ver4.T)
    Info +=
        fmt::format("({0}, {1}, {2}, {3}), ", t[INL], t[OUTL], t[INR], t[OUTR]);
  Info += "\n";
  Info += indent + fmt::format("├─K : ");
  auto &LegK = Ver4.LegK;
  Info += fmt::format("({0}, {1}, {2}, {3}),  K0: {4}", fmt::ptr(LegK[INL]),
                      fmt::ptr(LegK[OUTL]), fmt::ptr(LegK[INR]),
                      fmt::ptr(LegK[OUTR]), fmt::ptr(&Ver4.K[0]));
  Info += "\n";
  // Info += indent + fmt::format("└─\n");
  for (int p = 0; p < Ver4.Bubble.size(); p++) {
    Info += indent + ". │\n";
    bubble pp = Ver4.Bubble[p];
    Info += indent +
            fmt::format(". ├PAIR - Channel: {0}, LVerLoopNum: {1}, G: {2}\n",
                        ChanName[pp.Channel], pp.LVer.LoopNum,
                        Ver4.G[pp.Channel].size());

    // Info += indent + " . │\n";
    Info += indent + fmt::format(". ├─G[{0}, {1}]:", ChanName[pp.Channel],
                                 Ver4.G[pp.Channel].size());
    for (auto &g : Ver4.G[pp.Channel])
      Info += fmt::format("({0},{1}), ", g.T[IN], g.T[OUT]);

    Info += "\n";
    Info += indent + fmt::format(". ├─Map:");
    for (auto &m : pp.Map)
      Info += fmt::format("({0},{1}):{2}, ", m.LVerTidx, m.RVerTidx, m.Tidx);

    Info += "\n";
    // Info += "\n" + indent + ". │\n";
    Info += ToString(pp.LVer, indent + ". ");
    // Info +=
    //     indent +
    //     ".....................................................\n";

    // Info += indent + ". │\n";
    Info += ToString(pp.RVer, indent + ". ");

    ASSERT_ALLWAYS(pp.LVer.LegK[INL] == Ver4.LegK[INL],
                   "INL K address does not match! "
                       << pp.LVer.LegK[INL] << " vs " << Ver4.LegK[INL]);
    // Info += indent + ".  \n";
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
