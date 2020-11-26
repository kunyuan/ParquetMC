#include "testcode.h"

extern parameter Para;
namespace testcode{


void TestCode(){

    std::vector<double> res;
    int num = 10000;
    for (int i = 0; i <= num; i++)
    {
        double T = 0.0 + 1.0*i/num*Para.Beta;
        int Tidx = Para.TauGrid.floor(T);
        double interpT = Para.TauGrid.grid[Tidx];
        res.push_back(interpT);
    }
    
    std::string fname = "tauinterp.data";
    std::ofstream interpFile;
    interpFile.open(fname, std::ios::out|std::ios::trunc);
    if (interpFile.is_open()){
        for (int i = 0; i < res.size(); i++)
            interpFile << res[i] << "   ";
    }
    interpFile.close();
    

    std::vector<double> resK;
    num = 10000;
    double kMax = Para.FermiKGrid.grid.back();
    for (int i = 0; i <= num; i++)
    {
        double k = -1.5*kMax + 1.0*i/num*(3*kMax);
        int Kidx = Para.FermiKGrid.floor(k);
        double interpK = Para.FermiKGrid.grid[Kidx];
        resK.push_back(interpK);
    }

    fname = "kinterp.data";
    interpFile;
    interpFile.open(fname, std::ios::out|std::ios::trunc);
    if (interpFile.is_open()){
        for (int i = 0; i < res.size(); i++)
            interpFile << resK[i] << "   ";
    }
    interpFile.close();


    std::exit(0);
}


template <typename T> int GetArraySize(T& array)
{
    return (sizeof(array)/sizeof(array[0]));
}


template <typename T> void Print1DArray(T& array)
{
    for (size_t i = 0; i < GetArraySize(array); i++)
        std::cout << i << "th:  " << array[i] << std::endl;
}


template <typename T> void Print2DArray(T& array)
{
    int iLen = GetArraySize(array);
    int jLen = GetArraySize(array[0]);
    for (size_t i = 0; i < iLen; i++)
        for (size_t j = 0; j < jLen; j++)
            std::cout << "["<<i<<", "<<j<<"]: "  << array[i][j] << std::endl;
}

};   // namespace testcode;
