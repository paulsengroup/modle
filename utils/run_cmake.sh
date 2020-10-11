#!/usr/bin/env bash

set -e

wd="$(dirname "$(dirname "$(readlink -f "$0")")")"

cd "$wd"
mkdir -p cmake-build-debug-gcc10 cmake-build-release-gcc10 cmake-build-debug-llvm10 cmake-build-release-llvm10

cd cmake-build-debug-gcc10
/usr/bin/cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DMODLE_BUILD_DOCS=ON -G "CodeBlocks - Unix Makefiles" ..
cd ../cmake-build-release-gcc10
/usr/bin/cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DMODLE_BUILD_DOCS=ON -G "CodeBlocks - Unix Makefiles" ..


cd ../cmake-build-debug-llvm10
/usr/bin/cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DMODLE_BUILD_DOCS=ON -G "CodeBlocks - Unix Makefiles" ..
cd ../cmake-build-release-llvm10
/usr/bin/cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DMODLE_BUILD_DOCS=ON -G "CodeBlocks - Unix Makefiles" ..
