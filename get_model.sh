# -*- mode: sh; -*-
#!/bin/sh

echo $1
echo "\n" > field_temp.cl
sed -n "/#pragma start_$1/,/#pragma end_$1/p" src/model.h >> field_temp.cl
echo "\n" >> field_temp.cl
sed -n "/#pragma start_$1/,/#pragma end_$1/p" src/model.c >> field_temp.cl

sed -i "/start_$1/d" field_temp.cl
sed -i "/end_$1/d" field_temp.cl

echo "\n" >> field_temp.cl
sed -n -i "/#pragma start_opencl/,/#pragma end_opencl/p" field_temp.cl
sed -i "/start_opencl/d" field_temp.cl
sed -i "/end_opencl/d" field_temp.cl

sed -i -e "s/$1//g" field_temp.cl

rm field_temp.cl
return 0
