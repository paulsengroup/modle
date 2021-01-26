#!/usr/bin/env bash

set -e

CMAKE_VER='3.11.4'
#CMAKE_VER='3.19.2'

rsync -aPv hovig:scratch/clion-eap/modle/cmake-builds/cmake-build-release-gcc8-hovig/CMakeCache.txt ~/github/modle/cmake-builds/cmake-build-release-gcc8-hovig/
rsync -aPv hovig:scratch/clion-eap/modle/cmake-builds/cmake-build-release-gcc8-hovig/CMakeFiles/TargetDirectories.txt ~/github/modle/cmake-builds/cmake-build-release-gcc8-hovig/CMakeFiles/
rsync -aPv hovig:scratch/clion-eap/modle/cmake-builds/cmake-build-release-gcc8-hovig/CMakeFiles/Makefile.cmake ~/github/modle/cmake-builds/cmake-build-release-gcc8-hovig/CMakeFiles/
rsync -aPv hovig:scratch/clion-eap/modle/cmake-builds/cmake-build-release-gcc8-hovig/CMakeFiles/$CMAKE_VER/CMakeCXXCompiler.cmake ~/github/modle/cmake-builds/cmake-build-release-gcc8-hovig/CMakeFiles/$CMAKE_VER/

rsync -aPv hovig:scratch/clion-eap/modle/cmake-builds/cmake-build-debug-gcc8-hovig/CMakeCache.txt ~/github/modle/cmake-builds/cmake-build-debug-gcc8-hovig/
rsync -aPv hovig:scratch/clion-eap/modle/cmake-builds/cmake-build-debug-gcc8-hovig/CMakeFiles/TargetDirectories.txt ~/github/modle/cmake-builds/cmake-build-debug-gcc8-hovig/CMakeFiles/
rsync -aPv hovig:scratch/clion-eap/modle/cmake-builds/cmake-build-debug-gcc8-hovig/CMakeFiles/Makefile.cmake ~/github/modle/cmake-builds/cmake-build-debug-gcc8-hovig/CMakeFiles/
rsync -aPv hovig:scratch/clion-eap/modle/cmake-builds/cmake-build-debug-gcc8-hovig/CMakeFiles/$CMAKE_VER/CMakeCXXCompiler.cmake ~/github/modle/cmake-builds/cmake-build-debug-gcc8-hovig/CMakeFiles/$CMAKE_VER/

rsync -aPv hovig:scratch/clion-eap/modle/cmake-builds/cmake-build-debug-gcc8-waserr-hovig/CMakeCache.txt ~/github/modle/cmake-builds/cmake-build-debug-gcc8-waserr-hovig/
rsync -aPv hovig:scratch/clion-eap/modle/cmake-builds/cmake-build-debug-gcc8-waserr-hovig/CMakeFiles/TargetDirectories.txt ~/github/modle/cmake-builds/cmake-build-debug-gcc8-waserr-hovig/CMakeFiles/
rsync -aPv hovig:scratch/clion-eap/modle/cmake-builds/cmake-build-debug-gcc8-waserr-hovig/CMakeFiles/Makefile.cmake ~/github/modle/cmake-builds/cmake-build-debug-gcc8-waserr-hovig/CMakeFiles/
rsync -aPv hovig:scratch/clion-eap/modle/cmake-builds/cmake-build-debug-gcc8-waserr-hovig/CMakeFiles/$CMAKE_VER/CMakeCXXCompiler.cmake ~/github/modle/cmake-builds/cmake-build-debug-gcc8-waserr-hovig/CMakeFiles/$CMAKE_VER/
