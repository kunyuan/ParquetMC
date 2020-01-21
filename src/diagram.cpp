#include "diagram.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/utility.h"
#include <iostream>
#include <sstream>

using namespace diag;
using namespace std;

extern parameter Para;

/*Find Case Insensitive Sub String in a given substring */
size_t _findCaseInsensitive(std::string data, std::string toSearch,
                            size_t pos = 0) {
  // Convert complete given String to lower case
  std::transform(data.begin(), data.end(), data.begin(), ::tolower);
  // Convert complete given Sub String to lower case
  std::transform(toSearch.begin(), toSearch.end(), toSearch.begin(), ::tolower);
  // Find sub string in given string
  return data.find(toSearch, pos);
}

string _CheckKeyWord(istream &file, string KeyWord) {
  // a line only contains space will be skipped
  string buff;
  do {
    getline(file, buff);
    if (file.bad())
      ABORT("Fail to read the file!");
  } while (buff.find_first_not_of(' ') == buff.npos);

  auto found = _findCaseInsensitive(buff, KeyWord);
  ASSERT_ALLWAYS(found != string::npos,
                 fmt::format("{0} is not in: {1}", KeyWord, buff));

  return buff;
}

template <typename T>
vector<T> _ExtractNumbers(istream &file, string KeyWord = "") {
  // This function extracts the first T type number in one line in the file
  // stream.
  string line;
  getline(file, line);

  if (KeyWord.size() > 0) {
    auto found = _findCaseInsensitive(line, KeyWord);
    ASSERT_ALLWAYS(found != string::npos,
                   fmt::format("{0} is not in: {1}", KeyWord, line));
  }

  stringstream ss;
  vector<T> IntList;
  /* Storing the whole string into string stream */
  ss << line;
  /* Running loop till the end of the stream */
  string temp;
  T found;
  while (!ss.eof()) {
    /* extracting word by word from stream */
    ss >> temp;
    /* Checking the given word is integer or not */
    if (stringstream(temp) >> found)
      IntList.push_back(found);
    /* To save from space at the end of string */
    temp = "";
  }
  return IntList;
}

vector<loop> _Transpose(vector<vector<double>> &Basis) {
  vector<loop> NewBasis;
  if (Basis.size() == 0)
    return NewBasis;

  for (int i = 0; i < Basis[0].size(); i++) {
    loop TempLoopBasis;
    TempLoopBasis.fill(0);

    for (int j = 0; j < Basis.size(); j++)
      TempLoopBasis[j] = Basis[j][i];

    NewBasis.push_back(TempLoopBasis);
  }
  return NewBasis;
}

group diag::ReadOneGroup(istream &DiagFile) {
  group Group;
  _CheckKeyWord(DiagFile, "Type"); // group type, simply skip

  Group.HugenNum = _ExtractNumbers<int>(DiagFile, "DiagNum")[0];

  Group.Order = _ExtractNumbers<int>(DiagFile, "Order")[0];

  Group.GNum = _ExtractNumbers<int>(DiagFile, "GNum")[0];

  Group.Ver4Num = _ExtractNumbers<int>(DiagFile, "Ver4Num")[0];

  Group.LoopNum = _ExtractNumbers<int>(DiagFile, "LoopNum")[0];

  vector<int> ExtLoop = _ExtractNumbers<int>(DiagFile, "ExtLoopIndex");
  Group.ExtLoopNum = ExtLoop.size();
  Group.IsExtLoop.fill(false);
  for (auto index : ExtLoop)
    Group.IsExtLoop[index] = true;

  vector<int> ExtTransferLoop =
      _ExtractNumbers<int>(DiagFile, "ExtTransferLoopIndex");
  Group.ExtTransferLoopNum = ExtTransferLoop.size();
  Group.IsExtTransferLoop.fill(false);
  for (auto index : ExtTransferLoop)
    Group.IsExtTransferLoop[index] = true;

  vector<int> ExtLegLoop = _ExtractNumbers<int>(DiagFile, "ExtLegLoopIndex");
  Group.ExtLegLoopNum = ExtLegLoop.size();
  Group.IsExtLegLoop.fill(false);
  for (auto index : ExtLegLoop)
    Group.IsExtLegLoop[index] = true;

  vector<int> LockedLoop = _ExtractNumbers<int>(DiagFile, "DummyLoopIndex");
  Group.IsLockedLoop.fill(false);
  for (auto index : LockedLoop)
    Group.IsLockedLoop[index] = true;

  Group.TauNum = _ExtractNumbers<int>(DiagFile, "TauNum")[0];
  ASSERT_ALLWAYS(Group.TauNum <= MaxTauNum,
                 "Tau Number must be smaller than " << MaxTauNum);

  cout << Group.Order << ", " << Group.TauNum << " " << Group.LoopNum << endl;

  vector<int> ExtTau = _ExtractNumbers<int>(DiagFile, "ExtTauIndex");
  Group.ExtTauNum = ExtTau.size();
  Group.IsExtTau.fill(false);
  for (auto index : ExtTau)
    Group.IsExtTau[index] = true;

  vector<int> LockedTau = _ExtractNumbers<int>(DiagFile, "DummyTauIndex");
  Group.IsLockedTau.fill(false);
  for (auto index : LockedTau)
    Group.IsLockedTau[index] = true;

  Group.InternalLoopNum = Group.LoopNum - Group.ExtLoopNum;
  Group.InternalTauNum = Group.TauNum - Group.ExtTauNum;

  // diag::Test(Group);
  //   cout << "Group" << endl;
  //   for (int i = 0; i < Group.GNum; i++)
  //   cout << "in Group" << endl;
  //   cout << ToString(Group) << endl;
  //   cout << ToString(*(Group.DiagList[0].GIndex[0])) << endl;

  return Group;
}

std::string ToString(const group &Group) {
  return fmt::format("GroupID: {0}\n Weight={1}\n NewWeight={2}\n "
                     "HugenNum={3}\n ReWeight={4}\n",
                     Group.ID, Group.Weight, Group.NewWeight, Group.HugenNum,
                     Group.ReWeight);
}