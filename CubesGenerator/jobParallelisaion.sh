#!/bin/bash

pathFiles="/trinity/home/pturk/ScriptCubeGenerator"

for file in $pathFiles/cube*.sh
do

sbatch $file

done