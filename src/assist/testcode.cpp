#include "testcode.h"

extern parameter Para;
namespace testcode{


void TestCode(){
    std::cout << "Beta = " << Para.Beta << std::endl;
    std::cout << "TauGrid.size = " << Para.TauGrid.size << std::endl;

    for (size_t i = 0; i < Para.TauGrid.size; i++)
    {
        std::cout << Para.TauGrid.grid[i] << std::endl;
    }
    

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
