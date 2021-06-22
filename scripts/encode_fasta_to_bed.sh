#!/usr/bin/env bash

set -e

# TODO: fix code to support chrX/Y

zcat "$1" |
    grep '>' |
    awk -v OFS='\t' '{ print $3,$1 }' |
    sed 's/>//g' |
    awk -v OFS='\t' -F ":" '{ print $3,$4,$5,$6}' |
    awk -v OFS='\t' '{ print $1,$2,$3,$5,"0.0",$4 }' |
    awk -v OFS='\t' '{ gsub("-1", "-", $6); gsub("1", "+", $6); print }' |
    sed -E 's/^([[:digit:]]+.*)/chr\1/g'
