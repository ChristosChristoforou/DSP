#!/bin/bash

echo -e "\n" >> debug.txt &&
date >> debug.txt &&
gcc -g my_fft.c -o my_fft -lm 2>> debug.txt &&

echo -e "\n" >> debug.txt &&
valgrind ./my_fft --leak-check=full --track-origins=yes 2>> debug.txt &&

echo -e "\n" >> FFT.txt &&
date >> FFT.txt &&
./my_fft >> FFT.txt &&
clear


