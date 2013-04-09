#!/bin/bash

OUT=CD.c

mv CD.c bak/CD_bak.c
awk -f combine.awk CD_bak.c > tmp.txt

# remove all whitespace (including tabs) from left to first word
# http://www.cyberciti.biz/tips/delete-leading-spaces-from-front-of-each-word.html
awk -f awk.awk tmp.txt | sed -e 's/^[ \t]*//' > ${OUT}
echo '' >> ${OUT}
cat bak/CD_bak.c >> ${OUT}

rm tmp.txt
