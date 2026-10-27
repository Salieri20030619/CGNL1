#!/bin/bash
cd '/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/epithelial_independence_20261012_v1'
/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006/env/methods/bin/python scripts/01_extract.py > logs/extract.log 2>&1
rc=$?
printf '%s\n' "$rc" > logs/extract.exit
exit "$rc"
