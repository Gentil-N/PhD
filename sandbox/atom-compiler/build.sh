#!/bin/sh

echo "------------------------------ Building ------------------------------"

if [ "$1" = "debug" ]
then
    echo "------------------------------ Debug ------------------------------"
    mkdir -p build/debug
    cd build/debug
    cmake ../../. -DCMAKE_BUILD_TYPE="Debug"
elif [ "$1" = "release" ]
then
    echo "------------------------------ Release ------------------------------"
    mkdir -p build/release
    cd build/release
    cmake ../../. -DCMAKE_BUILD_TYPE="Release"
else
    echo "Build type is not valid. Only debug or release."
    exit 1
fi

make -j6

echo "------------------------------ Running ------------------------------"

./atomcomp
