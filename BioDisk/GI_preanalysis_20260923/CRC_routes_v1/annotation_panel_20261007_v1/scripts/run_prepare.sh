#!/bin/bash
/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006/env/methods/bin/python -u /mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/annotation_panel_20261007_v1/scripts/01_prepare.py > /mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/annotation_panel_20261007_v1/logs/prepare.log 2>&1
status=$?
printf "%s\n" "$status" > /mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/annotation_panel_20261007_v1/logs/prepare.exit
exit "$status"
