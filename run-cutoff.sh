#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-03-05 20:24:23
# @DESCRIPTION:

# Number of input parameters

# jrocker ex 8686

root_dir=/home/liucj/github/THPOC-code
cutoff_dir=/home/liucj/tmp/THPOC-cutoff
data_dir=/workspace/liucj/project/09-THPOC
cd ${root_dir}


function fn_run_modeling {
  cf=$1
  for r in `ls modeling/*.R`;
  do

    logname=`basename ${r}`
    echo "Start running ${r} ${cf}."
    cmd="Rscript ${r} ${cf} 1> data/logs/${logname}.log 2> data/logs/${logname}.err"
    eval ${cmd}
    # echo ${cmd}
    echo "${r} running done."
  done
}


for cutoff in `seq 0.3 0.01 0.5`;
do
  rm ${root_dir}/data
  rm -rf ${cutoff_dir}/${cutoff}
  mkdir -p ${cutoff_dir}/${cutoff}/{logs,output,rda}
  ln -sf ${cutoff_dir}/${cutoff} ${root_dir}/data
  ln -sf ${data_dir}/rda/wuhan.se.rds.gz ${root_dir}/data/rda/wuhan.se.rds.gz
  ln -sf ${data_dir}/rda/tom.se.rds.gz ${root_dir}/data/rda/tom.se.rds.gz

  fn_run_modeling ${cutoff}
done

