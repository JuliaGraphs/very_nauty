#/bin/bash

for f in *c;
do
  #echo $f
  rm -f ${f%.c}
done
