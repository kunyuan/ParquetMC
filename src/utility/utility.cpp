//---------------------------------------------------------------------------
// �Nikolai V. Shokhirev, 2004-2008  <nikolai@shokhirev.com>
// http://www.shokhirev.com/nikolai.html
// Reduced demo version
//---------------------------------------------------------------------------

#include "utility.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>

using namespace std;

//---------------------------------------------------------------------------

bool IsEqual(double x1, double x2, double eps) { return (fabs(x1 - x2) < eps); }
bool IsEqual(uint x1, uint x2, double eps) { return x1 == x2; }
bool IsEqual(int x1, int x2, double eps) { return x1 == x2; }
// FORTRAN abs
// double fabs(double x) { return ((x >= 0.0) ? x : -x); }

// FORTRAN iabs
// int iabs(int x) { return ((x >= 0.0) ? x : -x); }
// more functions ...

bool IsZero(double x, double eps) { return (fabs(x) < eps); }

bool CleanFile(const string &FileName) {
  ofstream ofs(FileName, std::ofstream::out | std::ofstream::trunc);
  if (ofs.is_open()) {
    ofs.close();
    return true;
  } else {
    ofs.close();
    return false;
  }
}

bool DoesFileExist(const string &FileName) {
  ofstream ofs(FileName, std::ofstream::in);
  if (ofs.is_open()) {
    ofs.close();
    return true;
  } else {
    ofs.close();
    return false;
  }
}

std::string ToString(const double x, const int width, const int decDigits) {
  stringstream ss;
  ss << fixed << right;
  ss.fill(' ');            // fill space around displayed #
  ss.width(width);         // set  width around displayed #
  ss.precision(decDigits); // set # places after decimal
  ss << x;
  return ss.str();
}
std::string Center(const string s, const int w) {
  stringstream ss, spaces;
  auto padding = w - s.size(); // count excess room to pad
  for (int i = 0; i < padding / 2; ++i)
    spaces << " ";
  ss << spaces.str() << s << spaces.str(); // format with padding
  if (padding > 0 && padding % 2 != 0)     // if odd #, add 1 space
    ss << " ";
  return ss.str();
}
string Right(const string s, const int w) {
  stringstream ss, spaces;
  auto padding = w - s.size(); // count excess room to pad
  for (int i = 0; i < padding; ++i)
    spaces << " ";
  ss << spaces.str() << s; // format with padding
  return ss.str();
}

string Left(const string s, const int w) {
  stringstream ss, spaces;
  auto padding = w - s.size(); // count excess room to pad
  for (int i = 0; i < padding; ++i)
    spaces << " ";
  ss << s << spaces.str(); // format with padding
  return ss.str();
}

std::string ProgressBar(double progress) {
  std::ostringstream oss;
  int barWidth = 70;

  oss << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos)
      oss << "I";
    else if (i == pos)
      oss << "I";
    else
      oss << " ";
  }
  oss << "] " << int(progress * 100.0) << " %\r";
  oss.flush();
  return oss.str();
}

/**
 * Remove surrounding whitespace from a std::string.
 * @param s The string to be modified.
 * @param t The set of characters to delete from each end
 * of the string.
 * @return The same string passed in as a parameter reference.
 */
std::string &trim(std::string &s, const char *t) {
  s.erase(0, s.find_first_not_of(t));
  s.erase(s.find_last_not_of(t) + 1);
  return s;
}

stringstream GetLine(ifstream &File) {
  string line;
  while (true) {
    getline(File, line);
    line = trim(line);
    if (line.size() > 0 && line[0] != '#') {
      replace(line.begin(), line.end(), ',', ' ');
      return stringstream(line);
    }
  }
}
