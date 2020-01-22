//
//  markov.cpp
//
//  Created by Kun Chen on 1/21/19.
//  Copyright (c) 2019 Kun Chen. All rights reserved.
//
#include "markov.h"
#include "utility/fmt/format.h"
#include "utility/fmt/printf.h"
#include <iostream>

extern parameter Para;
extern RandomFactory Random;

using namespace mc;
using namespace diag;
using namespace std;

int GetTauNum(int Order) { return Order + 1; }
int GetLoopNum(int Order) { return Order; }

void markov::ChangeOrder() {
  Updates Name;
  double Prop = 1.0;
  int NewOrder;

  if (Random.urn() < 0.5) {
    // increase order
    Name = INCREASE_ORDER;
    if (Var.CurrOrder == Para.Order)
      return;
    NewOrder = Var.CurrOrder + 1;
    static momentum NewMom;
    double NewTau;
    // Generate New Tau
    Prop = GetNewTau(NewTau);
    int NewTauIndex = GetTauNum(Var.CurrOrder);
    Var.Tau[NewTauIndex] = NewTau;
    // Generate New Mom
    Prop *= GetNewK(NewMom);
    Var.LoopMom[GetLoopNum(Var.CurrOrder)] = NewMom;

    // if the current order is zero, then set channel of order 1 at T
    if (Var.CurrOrder == 0)
      Var.CurrChannel = dse::T;

  } else {
    // decrese order
    if (Var.CurrOrder == 0)
      return;

    // if the current order is one, then decrease order is possible only for T
    if (Var.CurrOrder == 1 && Var.CurrChannel != dse::T)
      return;

    Name = DECREASE_ORDER;
    NewOrder = Var.CurrOrder - 1;
    // Remove OldTau
    int TauToRemove = GetTauNum(Var.CurrOrder) - 1;
    Prop = RemoveOldTau(Var.Tau[TauToRemove]);
    // Remove OldMom
    int LoopToRemove = GetLoopNum(Var.CurrOrder) - 1;
    Prop *= RemoveOldK(Var.LoopMom[LoopToRemove]);
  }
  Proposed[Name][Var.CurrOrder] += 1;

  // Weight.ChangeGroup(NewGroup);
  NewWeight = Weight.Evaluate(NewOrder, Var.CurrChannel);
  NewAbsWeight = fabs(NewWeight.Sum());
  double R = Prop * NewAbsWeight * Para.ReWeight[NewOrder] / Var.CurrAbsWeight /
             Para.ReWeight[Var.CurrOrder];

  // if (Name == INCREASE_ORDER) {
  //   cout << "time:" << Var.Tau[2] << ", " << Var.Tau[3] << endl;
  // cout << NewWeight << endl;
  // }
  // if (NewGroup.Order == 2)
  //   cout << NewWeight << endl;

  if (Random.urn() < R) {
    Accepted[Name][Var.CurrOrder]++;
    Var.CurrVersion++;
    Var.CurrOrder = NewOrder;
    Var.CurrWeight = NewWeight;
    Var.CurrAbsWeight = NewAbsWeight;
  }
  return;
};

void markov::ChangeTau() {
  int TauIndex = Random.irn(0, GetTauNum(Var.CurrOrder) - 1);

  Proposed[CHANGE_TAU][Var.CurrOrder]++;

  double CurrTau = Var.Tau[TauIndex];
  double NewTau;
  double Prop = ShiftTau(CurrTau, NewTau);

  Var.Tau[TauIndex] = NewTau;

  NewWeight = Weight.Evaluate(Var.CurrOrder, Var.CurrChannel);
  NewAbsWeight = fabs(NewWeight.Sum());
  double R = Prop * NewAbsWeight / Var.CurrAbsWeight;

  if (Random.urn() < R) {
    Accepted[CHANGE_TAU][Var.CurrOrder]++;
    Var.CurrVersion++;
    Var.CurrWeight = NewWeight;
    Var.CurrAbsWeight = NewAbsWeight;
  } else {
    // retore the old Tau if the update is rejected
    // if TauIndex is external, then its partner can be different
    Var.Tau[TauIndex] = CurrTau;
  }
};

void markov::ChangeMomentum() {
  int LoopIndex = Random.irn(0, GetLoopNum(Var.CurrOrder) - 1);
  // int LoopIndex = int(Random.urn() * (Var.CurrGroup->Order + 3));

  // InL momentum is locked
  if (LoopIndex == 1)
    return;

  Proposed[CHANGE_MOM][Var.CurrOrder]++;

  double Prop;
  int NewExtMomBin;
  static momentum CurrMom;

  CurrMom = Var.LoopMom[LoopIndex];

  if (LoopIndex == 0) {
    // transfer momentum
    Prop = ShiftExtTransferK(Var.CurrExtMomBin, NewExtMomBin);
    Var.LoopMom[LoopIndex] = Para.ExtMomTable[NewExtMomBin];
    if (Var.LoopMom[LoopIndex].norm() > Para.MaxExtMom) {
      Var.LoopMom[LoopIndex] = CurrMom;
      return;
    }
  } else if (LoopIndex == 2) {
    // InR momentum
    Prop = ShiftExtLegK(CurrMom, Var.LoopMom[LoopIndex]);
  } else {
    Prop = ShiftK(CurrMom, Var.LoopMom[LoopIndex]);
  }

  NewWeight = Weight.Evaluate(Var.CurrOrder, Var.CurrChannel);
  NewAbsWeight = fabs(NewWeight.Sum());
  double R = Prop * NewAbsWeight / Var.CurrAbsWeight;

  if (Random.urn() < R) {
    Accepted[CHANGE_MOM][Var.CurrOrder]++;
    Var.CurrVersion++;
    Var.CurrWeight = NewWeight;
    Var.CurrAbsWeight = NewAbsWeight;
    if (LoopIndex == 0)
      Var.CurrExtMomBin = NewExtMomBin;
  } else {
    Var.LoopMom[LoopIndex] = CurrMom;
  }
};

// void markov::ChangeScale() {
//   if (Para.Type != RG)
//     return;
//   double Prop = 1.0;
//   double OldScale = Var.CurrScale;
//   if (Random.urn() < 0.5)
//     Var.CurrScale += Random.urn() * 0.5;
//   else
//     Var.CurrScale -= Random.urn() * 0.5;

//   if (Var.CurrScale < 0.0) {
//     Var.CurrScale = OldScale;
//     return;
//   }

//   Proposed[CHANGE_SCALE][Var.CurrGroup->ID] += 1;

//   // force to change the group weight
//   // Weight.ChangeGroup(*(Var.CurrGroup), true);
//   double NewWeight = Weight.GetNewWeight(*Var.CurrGroup);

//   double R = Prop * fabs(NewWeight) / fabs(Var.CurrGroup->Weight);
//   if (Random.urn() < R) {
//     Accepted[CHANGE_SCALE][Var.CurrGroup->ID]++;
//     Weight.AcceptChange(*Var.CurrGroup);
//   } else {
//     Var.CurrScale = OldScale;
//     Weight.RejectChange(*Var.CurrGroup);
//   }
//   return;
// }

void markov::ChangeChannel() {
  if (Var.CurrOrder == 0)
    return;
  double Prop = 1.0;
  int NewChannel = int(Random.urn() * 4);
  // Var.CurrChannel = dse::T;
  // if (Var.CurrChannel == dse::U) {
  //   Var.CurrChannel = OldChannel;
  //   return;
  // }

  Proposed[CHANGE_CHANNEL][Var.CurrOrder] += 1;

  NewWeight = Weight.Evaluate(Var.CurrOrder, NewChannel);
  NewAbsWeight = fabs(NewWeight.Sum());
  double R = Prop * NewAbsWeight / Var.CurrAbsWeight;
  if (Random.urn() < R) {
    Accepted[CHANGE_CHANNEL][Var.CurrOrder]++;
    Var.CurrVersion++;
    Var.CurrWeight = NewWeight;
    Var.CurrAbsWeight = NewAbsWeight;
    Var.CurrChannel = NewChannel;
  }
  return;
}
