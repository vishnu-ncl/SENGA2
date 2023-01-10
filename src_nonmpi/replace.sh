#!/bin/bash

array=(store7 trun transp prun itndex urhs vrhs wrhs drhs erhs yrhs)

arraylength=${#array[@]}

for (( i=0; i<${arraylength}; i++ ));
do
    x=$((i+1))
    echo "index: $x, value: ${array[$i]}"
    
    sed -i "s/${array[$i]}(ic,jc,kc)/${array[$i]}(OPS_ACC$x(0,0,0))/g" test
done
