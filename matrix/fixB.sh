#!/bin/bash

# run as
#$./fixB B3.c

# input file from command line
IN=$1
BAK=bak/${IN}_bak
OUT=$IN

mv $IN $BAK
awk -f combine.awk $BAK > tmp.txt

# remove all whitespace (including tabs) from left to first word
# http://www.cyberciti.biz/tips/delete-leading-spaces-from-front-of-each-word.html
awk -f awk.awk tmp.txt | sed -e 's/^[ \t]*//' > ${OUT}
echo '' >> ${OUT}
cat $BAK >> ${OUT}

rm tmp.txt
