#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-03-06 19:31:50
# @DESCRIPTION:

# Number of input parameters
root_dir=/home/liucj/github/THPOC-code
data_dir=/workspace/liucj/project/09-THPOC

cd ${root_dir}

rm -rf logs/0*
rm data
ln -s ${data_dir} data

Rscript ${root_dir}/analysis/01-counts2se.R 1> logs/01-counts2se.R.log 2> logs/01-counts2se.R.err

function fn_run_modeling {
  cf=$1
  for r in `ls modeling/*.R`;
  do

    logname=`basename ${r}`
    echo "Start running ${r} ${cf}."
    cmd="Rscript ${r} ${cf} 1> logs/${logname}.log 2> logs/${logname}.err"
    eval ${cmd}
    echo ${cmd}
    echo "${r} running done."
  done
}

fn_run_modeling 0.30