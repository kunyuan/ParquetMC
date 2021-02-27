#ifndef testcode_H
#define testcode_H

#include "../global.h"
#include <iostream>

namespace testcode{
    
    void TestCode();
    template <typename T> int GetArraySize(T& array);
    template <typename T> void Print1DArray(T& array);
    template <typename T> void Print2DArray(T& array);


};  //namespace testcode

#endif