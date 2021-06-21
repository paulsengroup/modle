#!/usr/bin/env bash

set -e

awk -v OFS='\t' '{ if ($6 == "+")  print $1,$2,$2,$4,$5,$6; else print $1,$3,$3,$4,$5,$6 }' "$1"
