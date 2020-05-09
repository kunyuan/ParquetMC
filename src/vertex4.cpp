#include "vertex4.h"
#include "utility/fmt/format.h"
#include <iostream>

using namespace diag;
using namespace std;

extern parameter Para;
extern variable Var;
extern propagator Prop;

int green::AddTidxPair(const array<int, 2> &T) {
  // find the T array in the list, if failed, create a new array
  for (int i = 0; i < _Tpair.size(); ++i) {
    auto &t = _Tpair[i];
    if (t[OUT] == T[OUT] && t[IN] == T[IN])
      return i;
  }
  _Tpair.push_back(T);
  //   cout << _Weight.size() << endl;
  _Weight.push_back(0.0);
  return _Tpair.size() - 1;
}

void green::Evaluate(bool IsAnomal) {
  int Size = _Tpair.size();
  for (int i = 0; i < Size; ++i) {
    auto &T = _Tpair[i];
    if (!IsAnomal)
      _Weight[i] = Prop.Green(Var.Tau[T[OUT]] - Var.Tau[T[IN]], K, UP);
    else
      _Weight[i] = Prop.F(Var.Tau[T[OUT]] - Var.Tau[T[IN]], K, UP);
    // if (Var.CurrOrder == 2)
    //   cout << T[IN] << "->" << T[OUT] << ": " << _Weight[i] << endl;
  }
}

void green::Evaluate(const momentum &_K, bool IsAnomal) {
  int Size = _Tpair.size();
  for (int i = 0; i < Size; ++i) {
    auto &T = _Tpair[i];
    if (!IsAnomal)
      _Weight[i] = Prop.Green(Var.Tau[T[OUT]] - Var.Tau[T[IN]], _K, UP);
    else
      _Weight[i] = Prop.F(Var.Tau[T[OUT]] - Var.Tau[T[IN]], _K, UP);
  }
}

void vertex4::Build(int level, int order, int loopIdx, int inTL,
                    const vector<channel> &chan, int side, bool inBox) {
  Level = level;
  Side = side;
  Tidx = inTL;
  LoopIdx = loopIdx;
  _InBox = inBox;
  Channel = chan;
  Order = order;

  if (LoopNum() == 0) {
    // the same for left and right vertex with loopnum=0
    _AddTidxPair({Tidx, Tidx, Tidx, Tidx});
  } else {
    vector<channel> UST, II;
    for (auto &c : Channel) {
      if (c == I)
        II.push_back(c);
      else
        UST.push_back(c);
    }
    // normal diagram
    // _BuildI(II);
    // UST channel
    for (auto c : UST)
      for (int ol = 0; ol < LoopNum(); ol++) {
        auto bubble = _BuildBubble(c, ol);
        if (bubble.Map.size() > 0)
          _UST.push_back(bubble);
      }
  }

  Weight.resize(Tpair.size());
  if (Level == 0 && Order >= 1)
    ChanWeight.resize(4);
}

int vertex4::_AddTidxPair(const array<int, 4> &T) {
  // find the T array in the list, if failed, create a new array
  for (int i = 0; i < Tpair.size(); i++) {
    auto t = Tpair[i];

    ASSERT_ALLWAYS(t[INL] == T[INL],
                   "left Tin must be the same for all subvertex!"
                       << t[INL] << " vs " << T[INL]);

    if (t[OUTL] == T[OUTL] && t[INR] == T[INR] && t[OUTR] == T[OUTR])
      return i;
  }
  Tpair.push_back(T);
  return Tpair.size() - 1;
}

bubble vertex4::_BuildBubble(channel chan, int ol) {
  vector<channel> FULL = {I, T, U, S, TC, UC};
  vector<channel> F = {I, U, S, TC, UC};
  vector<channel> V = {I, T, U, TC, UC};
  vector<channel> FULL_CT = {I, T, TC};
  vector<channel> F_CT = {I, TC};

  ASSERT_ALLWAYS(chan != I, "BuildUST can not process I channel!");
  ASSERT_ALLWAYS(ol < LoopNum(),
                 "LVer order must be smaller than the Ver order!");

  bubble bub;
  bub.Channel = chan;
  ////////////////////   Right SubVer  ///////////////////
  int lvl = Level + 1;
  int oR = Order - 1 - ol;
  int RInT = Tidx + (ol + 1);
  int Llopidx = LoopIdx + 1, Rlopidx = LoopIdx + 1 + ol;

  switch (chan) {
  case T:
    if (!_InBox) {
      bub.LVer.Build(lvl, ol, Llopidx, Tidx, F, LEFT, false);
      bub.RVer.Build(lvl, oR, Rlopidx, RInT, FULL, RIGHT, false);
    } else {
      bub.LVer.Build(lvl, ol, Llopidx, Tidx, F_CT, LEFT, true);
      bub.RVer.Build(lvl, oR, Rlopidx, RInT, FULL_CT, RIGHT, true);
    }
    break;
  case U:
    // continue;
    ASSERT_ALLWAYS(!_InBox,
                   "Ver4 in box can't have U diagram! Level: " << Level);
    bub.LVer.Build(lvl, ol, Llopidx, Tidx, F, LEFT, false);
    bub.RVer.Build(lvl, oR, Rlopidx, RInT, FULL, RIGHT, false);
    break;
  case S:
    // continue;
    ASSERT_ALLWAYS(!_InBox,
                   "Ver4 in box can't have S diagram! Level: " << Level);
    bub.LVer.Build(lvl, ol, Llopidx, Tidx, V, LEFT, false);
    bub.RVer.Build(lvl, oR, Rlopidx, RInT, FULL, RIGHT, false);
    break;
  case TC:
  case UC:
    // continue;
    bub.LVer.Build(lvl, ol, Llopidx, Tidx, F_CT, LEFT, true);
    bub.RVer.Build(lvl, oR, Rlopidx, RInT, FULL_CT, RIGHT, true);
    break;
  default:
    ABORT("The channel does not exist!");
    break;
  }

  ////// construct a map from LVer and RVer Tidx to the Ver Tidx //////
  for (int lt = 0; lt < bub.LVer.Tpair.size(); ++lt)
    for (int rt = 0; rt < bub.RVer.Tpair.size(); ++rt) {

      auto &LvT = bub.LVer.Tpair[lt];
      auto &RvT = bub.RVer.Tpair[rt];

      int GT0idx = G[0].AddTidxPair({LvT[OUTR], RvT[INL]});
      int GTidx, VerTidx;

      switch (chan) {
      case T:
        VerTidx = _AddTidxPair({LvT[INL], LvT[OUTL], RvT[INR], RvT[OUTR]});
        GTidx = G[chan].AddTidxPair({RvT[OUTL], LvT[INR]});
        break;
      case U:
        VerTidx = _AddTidxPair({LvT[INL], RvT[OUTR], RvT[INR], LvT[OUTL]});
        GTidx = G[chan].AddTidxPair({RvT[OUTL], LvT[INR]});
        break;
      case S:
        VerTidx = _AddTidxPair({LvT[INL], RvT[OUTL], LvT[INR], RvT[OUTR]});
        GTidx = G[chan].AddTidxPair({LvT[OUTL], RvT[INR]});
        break;
      case TC:
      case UC:
        // counterterms are equal-time
        VerTidx = _AddTidxPair({LvT[INL], LvT[INL], LvT[INL], LvT[INL]});
        GTidx = G[chan].AddTidxPair({RvT[OUTL], LvT[INR]});
        break;
      default:
        ABORT("The channel does not exist! " << chan);
        break;
      }

      bub.Map.push_back({lt, rt, GT0idx, GTidx, VerTidx});
    }

  return bub;
}

string vertex4::ToString(string indent) {
  string SideStr;
  if (Side == LEFT)
    SideStr = "LEFT";
  else
    SideStr = "RIGHT";
  // string Info =
  //     indent +
  //     "==============================================================\n";
  string Info =
      indent +
      fmt::format("├Level: {0}, {1}, LoopNum: {2}, Boxed: {3}, G0:{4}\n", Level,
                  SideStr, LoopNum(), _InBox, G[0].Size());
  Info += indent + fmt::format("├─Tau[{0}]: ", Tpair.size());
  for (auto &t : Tpair)
    Info +=
        fmt::format("({0}, {1}, {2}, {3}), ", t[INL], t[OUTL], t[INR], t[OUTR]);
  Info += "\n";

  ASSERT_ALLWAYS(Tpair.size() == Weight.size(),
                 "Tpair size must be equal to Weight size!");
  // Info += indent + fmt::format("├─LegK : ");
  // auto &LegK = Ver4.LegK;
  // Info +=
  //     fmt::format("({}, {}, {}, {})", fmt::ptr(LegK[INL]),
  //     fmt::ptr(LegK[OUTL]),
  //                 fmt::ptr(LegK[INR]), fmt::ptr(LegK[OUTR]));
  // Info += "\n";
  // if (Ver4.LoopNum > 0) {
  //   Info += indent + fmt::format("├─InterK : ");
  //   auto &K = Ver4.K;
  //   Info += fmt::format("c {}, T {}, U {}, S {}", fmt::ptr(&K[0]),
  //                       fmt::ptr(&K[1]), fmt::ptr(&K[2]), fmt::ptr(&K[3]));
  //   Info += "\n";
  // }
  Info += indent + fmt::format("└─\n");
  for (int p = 0; p < _UST.size(); p++) {
    Info += indent + ". │\n";
    bubble pp = _UST[p];
    Info += indent +
            fmt::format(". ├PAIR - Channel: {0}, LVerLoopNum: {1}, G: {2}\n",
                        ChanName[pp.Channel], pp.LVer.LoopNum(),
                        G[pp.Channel].Size());

    // Info += indent + " . │\n";
    Info += indent + fmt::format(". ├─G[{0}, {1}]:", ChanName[pp.Channel],
                                 G[pp.Channel].Size());
    for (auto &t : G[pp.Channel]._Tpair)
      Info += fmt::format("({0},{1}), ", t[IN], t[OUT]);

    Info += "\n";
    Info += indent + fmt::format(". ├─Map:");

    ASSERT_ALLWAYS(pp.Map.size() == pp.LVer.Tpair.size() * pp.RVer.Tpair.size(),
                   "Map size != LVer.Tpair size *RVer.Tpair size!");

    for (auto &m : pp.Map)
      Info += fmt::format("({0},{1}):{2}, ", m[LVERT], m[RVERT], m[VERT]);

    Info += "\n";
    // Info += "\n" + indent + ". │\n";
    Info += pp.LVer.ToString(indent + ". ");
    // Info +=
    //     indent +
    //     ".....................................................\n";

    // Info += indent + ". │\n";
    Info += pp.RVer.ToString(indent + ". ");

    // ASSERT_ALLWAYS(pp.LVer.LegK[INL] == Ver4.LegK[INL],
    //                "INL K address does not match! "
    //                    << pp.LVer.LegK[INL] << " vs " << Ver4.LegK[INL]);
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