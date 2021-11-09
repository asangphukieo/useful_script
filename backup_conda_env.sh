#!/bin/bash

time_stamp=`date "+%m-%d-%Y %H:%M:%S"`

conda env list | awk '{print $1}' | grep -v -E "^$|^#" | while read -r env_name
do
    yml_file=${env_name}.yml
    echo "[$time_stamp] Exporting conda environment $env_name to $yml_file ..."
    conda env export -n $env_name > $yml_file
done

echo "[$time_stamp] Done!"
