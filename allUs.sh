#!/bin/bash 

rm U_vals.csv

for dir in freq*/; do
    wait `python find_U.py $dir`    
done
