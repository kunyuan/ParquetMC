#ifndef grid_H
#define grid_H

#include <array>
#include <string>
#include <vector>

namespace grid {
  using namespace std;

  class Coeff {
  public:
    void init(array<double, 2> _bound, array<double, 2> _idx, double _lambda,
              bool dense2sparse);

    // return the index of a given value
    int floor(double x) const;
    // return the Grid with a given index
    double grid(int i) const;

  private:
    double _a;
    double _b;
    array<double, 2> bound; // the dense and sparse side bounds
    array<double, 2> idx;   // the dense and sparse side index
    double lambda;          //>0 for dense to sparse; <0 for sparse to dense
  };

  class UniLog {
    // generating a unilog grid
    // N=n*(m+1) points
    // the interval is first divided into a log grid,
    //   0, alpha^m, alpha^(m-1), ... , alpha, 1
    // then each subinterval is divided into a uniform grid of n
    //
    // the actual interval is 2n/(2n+alpha-1) the given length
  public:
    void init(array<double,2> _bound, int _init, int _m, int _n, double _alpha);

    // return the index of a given value
    int floor(double x) const;
    // return the Grid with a given index
    double grid(int i) const;
    // for d2s, give bound[0] at i=idx[0];  bound[1]-(1-alpha)*length/n at i = idx[1] 

    //private:
    double length;
    int m;
    int n;
    array<double, 2> bound; 
    array<int, 2> idx;   
    double alpha;
    bool d2s=true;
  };


  class Tau {
  public:
    int size;
    vector<double> grid;
    vector<double> weight;

    vector<double> build(double beta, int _size, double scale);
    int floor(double x) const; // return the index of a given value
    string str();

  private:
    Coeff _coeff0;
    Coeff _coeff1;
  };

  class TauUL {
  public:
    int size;
    vector<double> grid;
    vector<double> weight;

    vector<double> build(double beta, int m,int n, double scale); // scale gives minimum interval
    int floor(double x) const; // return the index of a given value
    string str();

  private:
    UniLog _unilog0;
    UniLog _unilog1;
  };


  class FermiK {
  public:
    int size;
    int kFidx;
    vector<double> grid;
    void build(double kF, double maxK, int _size, double scale);
    int floor(double x) const; // return the index of a given value
    string str();

  private:
    Coeff _coeff0;
    Coeff _coeff1;
  };

  class FermiKUL {
  public:
    int size;
    int kFidx;
    vector<double> grid;
    void build(double kF, double maxK, int m, int n, double scale);
    int floor(double x) const; // return the index of a given value
    string str();

  private:
    UniLog _unilog0;
    UniLog _unilog1;
  };


  class BoseK {
  public:
    int size;
    int kFidx;
    int twokFidx;
    vector<double> grid;
    void build(double kF, double maxK, int _size, double scale);
    int floor(double x) const; // return the index of a given value
    string str();

  private:
    Coeff _coeff0;
    Coeff _coeff1;
    Coeff _coeff2;
  };

  class BoseKUL {
  public:
    int size;
    int kFidx;
    int twokFidx;
    vector<double> grid;
    void build(double kF, double maxK, int m,int n, double scale);
    int floor(double x) const; // return the index of a given value
    string str();

  private:
    UniLog _unilog0;
    UniLog _unilog1;
    UniLog _unilog2;
  };

  class Uniform {
  public:
    int size;
    double delta;
    double lowerBound;
    vector<double> weight;
    vector<double> grid;
    void build(std::array<double, 2> bounds, int _size);
    int floor(double x) const;
    std::string str();
  };

  vector<double> logGrid(const vector<Coeff> &Coeff,
                         const vector<array<int, 2>> &Range);

  vector<double> uniLogGrid(const vector<UniLog> &Unilog, const vector<array<int,2>> &Range);

} // namespace grid

#endif
