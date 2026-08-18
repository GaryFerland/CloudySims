#!/bin/sh
g++ $1/$1.cpp -I ../../source  ../../source/sys_gcc/libcloudy.a -o $1/$1.exe
