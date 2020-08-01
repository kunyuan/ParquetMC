#!/bin/bash
#compiler=clang++
compiler=icpc
#compiler=g++
#type=Debug
type=Release
mkdir build
cd build
if [ -n "$1" ]
  then
  if [ $1 = "-n" ] || [ $1 = "--new" ]; then
    make clean
    rm -rf * #force to rerun cmake configuration
    cmake -DCMAKE_BUILD_TYPE=$type -DCMAKE_CXX_COMPILER=$compiler ../src
  fi
fi
make -j
make install
cd -

g++ -I./lib -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` ./src/lib/grid.cpp ./pybind/grid.cpp -o ./grid.so
g++ -I./lib -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` ./src/lib/green.cpp ./pybind/green.cpp -o ./green.so
