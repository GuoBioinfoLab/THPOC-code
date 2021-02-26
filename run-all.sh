#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-03-05 20:24:23
# @DESCRIPTION:

# Number of input parameters

# jrocker ex 8686
cd /home/liucj/github/THPOC-code
echo "Workcing dir $PWD in docker container ${HOSTNAME}"

for r in `ls modeling/*.R`;
do

  logname=`basename ${r}`
  echo "Start running ${r}."
  cmd="Rscript ${r} 1> logs/${logname}.log 2> logs/${logname}.err"
  eval ${cmd}
  echo "${r} running done."
done
