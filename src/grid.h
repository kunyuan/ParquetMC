#ifndef grid_H
#define grid_H

#include <Eigen/Dense>
#include <array>
#include <string>
#include <vector>

class logGrid {
public:
  std::array<double, 2> Bound; // the dense and sparse side bounds
  std::array<double, 2> Idx;   // the dense and sparse side index

  double Lambda; //>0 for dense to sparse; <0 for sparse to dense

  void Initialize(std::array<double, 2> bound, std::array<double, 2> idx,
                  double lambda, bool dense2sparse);
  int Floor(double x); // return the index of a given value
  // return the Grid with a given index
  double Grid(int index);

private:
  double _a;
  double _b;
};

class tauGrid {
public:
  int Size;
  void Initialize(double beta, int size, double scale);
  int Floor(double x); // return the index of a given value
  Eigen::VectorXd Grid;
  Eigen::VectorXd Weight;
  std::string ToString();

private:
  logGrid _Grid0;
  logGrid _Grid1;
};

class kFermiGrid {
public:
  int Size;
  double MaxK;
  double kF;
  int kFIdx;
  void Initialize(double kF, double MaxK, int size, double scale);
  int Floor(double x); // return the index of a given value
  Eigen::VectorXd Grid;
  std::string ToString();

private:
  logGrid _Grid0;
  logGrid _Grid1;
};

class kBoseGrid {
public:
  int Size;
  double MaxK;
  double kF;
  int kFIdx;
  int TwokFIdx;
  void Initialize(double kF, double MaxK, int size, double scale);
  int Floor(double x); // return the index of a given value
  Eigen::VectorXd Grid;
  std::string ToString();

private:
  logGrid _Grid0;
  logGrid _Grid1;
  logGrid _Grid2;
};

class uniformGrid {
public:
  int Size;
  double Delta;
  double LowerBound;
  void Initialize(std::array<double, 2> bounds, int size);
  int Floor(double x);
  Eigen::VectorXd Grid;
  std::string ToString();
};

#endif