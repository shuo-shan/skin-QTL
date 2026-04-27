#!/bin/bash
cmd=$(sed -n "${LSB_JOBINDEX}p" "$1")
echo "RUN: $cmd"
eval $cmd
