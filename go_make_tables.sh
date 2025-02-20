#!/bin/bash

# ./make_data_tables.py < edge_chromatic_number.dat

f=edge_chromatic_number_connected.dat
for (( n=1; n<=12; n+=1 )) do
  echo n=$n >> $f
  geng -c $n | ./vn_graph_chi_dash >> $f
done

#f=edge_chromatic_number.dat
#for (( n=1; n<=10; n+=1 )) do
#  echo n=$n >> $f
#  geng $n | ./vn_graph_chi_dash >> $f
#done
