#!/bin/bash

# increase verbosity for debugging
set -x

mkdir -p $PREFIX/bin

conda env config vars set AWKLIBPATH="$PREFIX/lib/gawk"
conda env config vars set AWKPATH="$PREFIX/bin/"
cp $SRC_DIR/scripts/* $PREFIX/bin/* 
chmod +x $PREFIX/bin/*

