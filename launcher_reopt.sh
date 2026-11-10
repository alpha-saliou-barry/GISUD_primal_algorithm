#!/bin/bash

compiled_executable_path="/path/to/executable"
root_instances_folder="/path/to/folder"

script_path=$compiled_executable_path

paths=("$root_instances_folder/NW_320_0.25_0.25_0.25_0.25_1_day_16_day_23"
"$root_instances_folder/NW_320_0.25_0.25_0.25_0.25_1_day_1_day_9"
"$root_instances_folder/NW_320_0.34_0.33_0.33_1_day_1_day_9"
"$root_instances_folder/NW_320_0.34_0.33_0.33_1_day_16_day_23"
"$root_instances_folder/NW_757_0.25_0.25_0.25_0.25_1_day_1_day_11"
"$root_instances_folder/NW_757_0.25_0.25_0.25_0.25_1_day_12_day_22"
"$root_instances_folder/NW_757_0.34_0.33_0.33_1_day_12_day_22"
"$root_instances_folder/NW_757_0.34_0.33_0.33_1_day_21_day_31")

levels=("perturbation_0.2" "perturbation_0.3" "perturbation_0.4")
instances=("instance_0" "instance_1" "instance_2" "instance_3")

for i in ${!paths[@]}; do
    for j in ${!levels[@]}; do
        for k in ${!instances[@]}; do
            path="${paths[$i]}/${levels[$j]}/${instances[$k]}"
            cmd="${script_path} ${path} 0 1 0"
            echo $cmd
            $cmd

            cmd="${script_path} ${path} 0 0 0"
            echo $cmd
            $cmd

            cmd="${script_path} ${path} 0 1 1"
            echo $cmd
            $cmd
        done
    done
done