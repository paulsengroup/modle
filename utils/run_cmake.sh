#!/usr/bin/env bash

set -e

wd="$(dirname "$(dirname "$(readlink -f "$0")")")"

cd "$wd"
mkdir -p cmake-build-debug-gcc10 cmake-build-release-gcc10 cmake-build-relwithdebinfo-gcc10 \
         cmake-build-debug-llvm10 cmake-build-release-llvm10 cmake-build-relwithdebinfo-llvm10

cd cmake-build-debug-gcc10 && rm -rf ./*
/usr/bin/cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DMODLE_BUILD_DOCS=ON -G "Ninja" ..
cd ../cmake-build-release-gcc10 && rm -rf ./*
/usr/bin/cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DMODLE_BUILD_DOCS=ON -G "Ninja" ..
cd ../cmake-build-relwithdebinfo-gcc10 && rm -rf ./*
/usr/bin/cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DMODLE_BUILD_DOCS=ON -G "Ninja" ..


cd ../cmake-build-debug-llvm10 && rm -rf ./*
/usr/bin/cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DMODLE_BUILD_DOCS=ON -G "Ninja" ..
cd ../cmake-build-release-llvm10 && rm -rf ./*
/usr/bin/cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DMODLE_BUILD_DOCS=ON -G "Ninja" ..
cd ../cmake-build-relwithdebinfo-llvm10 && rm -rf ./*
/usr/bin/cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DMODLE_BUILD_DOCS=ON -G "Ninja" ..
