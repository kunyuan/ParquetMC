#ifndef grid_H
#define grid_H

#include <array>
#include <string>
#include <vector>

class grid {
public:
  std::array<double, 2> Bound;
  int Size;
  bool Dense2Sparse;
  double Lambda;

  void Initialize(std::array<double, 2> bound, int size, bool dense2Sparse,
                  double lambda);
  std::array<int, 2> Index(double x); // return the index of a given value
  // return the Grid with a given index
  double &Grid(int index);
  std::string ToString();

private:
  double _Gamma;
  double _Factor;
  std::vector<double> _Grid;
};

class tauGrid {
public:
  int Size;
  void Initialize(double beta, int size, double lambda);
  std::array<int, 2> Index(double x); // return the index of a given value
  double Grid(int Idx);               // return the Grid with a given index
  std::string ToString();

private:
  grid _Grid[2];
};

template <int FermiBose> class kGrid {
public:
  void Initialize(double kF, double MaxK);
  std::array<int, 2> Index(double x); // return the index of a given value
  double Grid(int Idx);               // return the Grid with a given index
  std::string ToString();

private:
  grid _Grid[3];
};

#endif