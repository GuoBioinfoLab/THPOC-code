#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-06-25 10:11:58
# @DESCRIPTION:

# Number of input parameters
path_root=/home/liucj/github/THPOC-code
path_analysis=${path_root}/analysis
path_modeling=${path_root}/modeling
path_src=${path_root}/src
rmd=${path_root}/test/scripts-word.Rmd
license=${path_root}/LICENSE

echo '
---
title: "TEP-code"
output:
  html_document:
    df_print: paged
  word_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

' > ${rmd}

echo '# LICENSE' >> ${rmd}
cat ${license} >> ${rmd}


srcs=`find ${path_analysis} ${path_modeling} ${path_src} -name "*.R" -type f | sort`

for src in ${srcs[@]};
do
  filename=`basename ${src}`

  header="# @AUTHOR: Chun-Jie Liu\n# @CONTACT: chunjie.sam.liu.at.gmail.com\n# @DATE: 2021-03-1 10:00:58\n# @DESCRIPTION: ${filename}\n"

  sed -i "1s/^/${header}/" ${src}

done



