#define FMT_HEADER_ONLY
#include "vertex4.h"
#include "utility/fmt/format.h"
#include <iostream>

using namespace tree;
using namespace std;

extern parameter Para;
extern variable Var;

kvector tree::KVector(int loopidx) {
  kvector k;
  k.Zero();
  k[loopidx] = 1;
  return k;
}

void vertex4::Build(int level, const set<channel> &chan,
                    const array<kvector, 4> &legK, int loopNum, int loopidx,
                    int inTL, int side) {
  Level = level;
  Side = side;
  Tidx = inTL;
  LoopIdx = loopidx;
  LoopNum = loopNum;
  LegK = legK;
  Channel = {};
  ChannelCT = {};
  if (LoopNum == 0) {
    // the same for left and right vertex with loopnum=0
    _AddTidxPair({Tidx, Tidx, Tidx, Tidx});
    return;
  }

  // push counter-term Tpair first
  for (auto &c : chan) {
    if (c == TC || c == UC) {
      ChannelCT.insert(c);
      _AddTidxPair({Tidx, Tidx, Tidx, Tidx});
    } else if (c == I) {
      Channel.insert(c);
    } else {
      // channel U, S, T
      Channel.insert(c);
      for (int ol = 0; ol < LoopNum; ol++)
        UST.push_back(_BuildBubble(c, ol));
    }
  }
}

int vertex4::_AddTidxPair(const t4pair &T) {
  // find the T array in the list, if failed, create a new array
  for (int i = 0; i < T4.size(); i++) {
    auto t = T4[i];
    ASSERT_ALLWAYS(t[INL] == T[INL],
                   "left Tin must be the same for all subvertex!"
                       << t[INL] << " vs " << T[INL]);
    if (t[OUTL] == T[OUTL] && t[INR] == T[INR] && t[OUTR] == T[OUTR])
      return i;
  }
  T4.push_back(T);
  return T4.size() - 1;
}

bubble vertex4::_BuildBubble(channel chan, int ol) {
  set<channel> FULL = {I, T, U, S, TC, UC};
  set<channel> F = {I, U, S, TC, UC};
  set<channel> V = {I, T, U, TC, UC};
  array<kvector, 4> LLegK, RLegK;

  ASSERT_ALLWAYS(chan != I, "BuildUST can not process I channel!");
  ASSERT_ALLWAYS(ol < LoopNum,
                 "LVer order must be smaller than the Ver order!");

  bubble bub;
  bub.Channel = chan;
  ////////////////////   Right SubVer  ///////////////////
  int lvl = Level + 1;
  int oR = LoopNum - 1 - ol;
  int RInT = Tidx + (ol + 1);
  int Llopidx = LoopIdx + 1, Rlopidx = LoopIdx + 1 + ol;

  auto K0 = KVector(LoopIdx);
  if (chan == T) {
    auto Kt = LegK[OUTL] + K0 - LegK[INL];

    LLegK = {LegK[INL], LegK[OUTL], Kt, K0};
    RLegK = {K0, Kt, LegK[INR], LegK[OUTR]};

    bub.LVer.Build(lvl, F, LLegK, ol, Llopidx, Tidx, LEFT);
    bub.RVer.Build(lvl, FULL, RLegK, oR, Rlopidx, RInT, RIGHT);
  } else if (chan == U) {
    auto Ku = LegK[OUTR] + K0 - LegK[INL];

    LLegK = {LegK[INL], LegK[OUTR], Ku, K0};
    RLegK = {K0, Ku, LegK[INR], LegK[OUTL]};

    bub.LVer.Build(lvl, F, LLegK, ol, Llopidx, Tidx, LEFT);
    bub.RVer.Build(lvl, FULL, RLegK, oR, Rlopidx, RInT, RIGHT);
  } else if (chan == S) {
    auto Ks = LegK[INL] + LegK[INR] - K0;

    LLegK = {LegK[INL], Ks, LegK[INR], K0};
    RLegK = {K0, LegK[OUTL], Ks, LegK[OUTR]};

    bub.LVer.Build(lvl, V, LLegK, ol, Llopidx, Tidx, LEFT);
    bub.RVer.Build(lvl, FULL, RLegK, oR, Rlopidx, RInT, RIGHT);
  } else
    ABORT("The channel does not exist!");

  auto &G0 = bub.G0;
  auto &Gx = bub.Gx;

  ////// construct a map from LVer and RVer Tidx to the Ver Tidx //////
  for (int lt = 0; lt < bub.LVer.T4.size(); ++lt)
    for (int rt = 0; rt < bub.RVer.T4.size(); ++rt) {

      auto &LvT = bub.LVer.T4[lt];
      auto &RvT = bub.RVer.T4[rt];
      int VerTidx;

      G0.push_back({LvT[OUTR], RvT[INL]});

      switch (chan) {
      case T:
        VerTidx = _AddTidxPair({LvT[INL], LvT[OUTL], RvT[INR], RvT[OUTR]});
        Gx.push_back({RvT[OUTL], LvT[INR]});
        break;
      case U:
        VerTidx = _AddTidxPair({LvT[INL], RvT[OUTR], RvT[INR], LvT[OUTL]});
        Gx.push_back({RvT[OUTL], LvT[INR]});
        break;
      case S:
        VerTidx = _AddTidxPair({LvT[INL], RvT[OUTL], LvT[INR], RvT[OUTR]});
        Gx.push_back({LvT[OUTL], RvT[INR]});
        break;
      default:
        ABORT("The channel does not exist! " << chan);
        break;
      }

      bub.Map.push_back({lt, rt, static_cast<int>(G0.size() - 1),
                         static_cast<int>(Gx.size() - 1), VerTidx});
    }

  return bub;
}

string compatK(const kvector &K) {
  stringstream ss;
  for (int i = 0; i < K.size(); ++i)
    ss << K[i];
  return ss.str();
}

string vertex4::_ToString(const set<channel> &chan) {
  stringstream ss;
  ss << "(";
  for (auto &c : chan)
    ss << c << ", ";
  ss << ")";
  return ss.str();
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
      indent + fmt::format("├Level: {0}, {1}, LoopNum: {2}, Counterterm: {3}\n",
                           Level, SideStr, LoopNum, _ToString(ChannelCT));
  Info += indent + fmt::format("├─Tau[{0}]: ", T4.size());
  for (auto &t : T4)
    Info +=
        fmt::format("({0}, {1}, {2}, {3}), ", t[INL], t[OUTL], t[INR], t[OUTR]);
  Info += "\n";

  Info += indent + fmt::format("├─LegK : ");
  Info +=
      fmt::format("({}, {}, {}, {})", compatK(LegK[INL]), compatK(LegK[OUTL]),
                  compatK(LegK[INR]), compatK(LegK[OUTR]));
  Info += "\n";
  // if (Ver4.LoopNum > 0) {
  //   Info += indent + fmt::format("├─InterK : ");
  //   auto &K = Ver4.K;
  //   Info += fmt::format("c {}, T {}, U {}, S {}", fmt::ptr(&K[0]),
  //                       fmt::ptr(&K[1]), fmt::ptr(&K[2]),
  //                       fmt::ptr(&K[3]));
  //   Info += "\n";
  // }
  Info += indent + fmt::format("└─\n");
  for (int p = 0; p < UST.size(); p++) {
    Info += indent + ". │\n";
    bubble pp = UST[p];
    Info += indent + fmt::format(". ├PAIR - Channel: {0}, LVerLoopNum: {1}\n",
                                 ChanName[pp.Channel], pp.LVer.LoopNum);

    // Info += indent + " . │\n";
    Info += indent + fmt::format(". ├─G0:");
    for (auto &t : pp.G0)
      Info += fmt::format("({0},{1}), ", t[IN], t[OUT]);
    Info += "\n";
    Info += indent + fmt::format(". ├─Gx:");
    for (auto &t : pp.Gx)
      Info += fmt::format("({0},{1}), ", t[IN], t[OUT]);
    Info += "\n";
    Info += indent + fmt::format(". ├─Map:");

    ASSERT_ALLWAYS(pp.Map.size() == pp.LVer.T4.size() * pp.RVer.T4.size(),
                   "Map size != LVer.Tpair size *RVer.Tpair size!");

    for (auto &m : pp.Map)
      Info += fmt::format("({0},{1}):{2}, ", m[LVER], m[RVER], m[VER]);

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