#!/bin/bash

#Preparing directories
mkdir -p ~/include/
mkdir -p ~/lib/
mkdir -p ~/bin/
mkdir -p build
cd build
platform=$(uname)

#Parsing arguments
if [ $# -gt 1 ]
then
    echo "Usage : sh install.sh [-o]"
    exit 1
fi
if [ $# -eq 1 ]
then
    if [ $1 = "-o" ]
    then
		rm -rf ~/lib/libElas.*
		rm -rf ~/lib/libCommons.*
	else
		echo "Usage : sh install.sh [-o]"
		exit 1
    fi
fi

#Installing Commons
if [ \( ! -f ~/lib/libCommons.so \) -a \( ! -f ~/lib/libCommons.dylib \) ]
then
    echo "-- Installing Commons"
    git clone https://github.com/ICStoolbox/Commons.git Commons
    mkdir Commons/build
    cd Commons/build
    cmake ..
    make install
    cd -
    rm -rf Commons/
else
	echo "-- Commons already installed. Skipping..."
fi

#Installing Elasticity library and executable
if [ \( ! -f ~/lib/libElas.so \)  -a \( ! -f ~/lib/libElas.dylib \) ]
then
    echo "-- Installing Elasticity"
    cmake ..
    make
    make install
else
	echo "-- Elasticity already installed. Skipping..."
fi

cd ../
