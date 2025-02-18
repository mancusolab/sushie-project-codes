#!/bin/bash

for idx in `seq 1 5`
do
    source score_warmup_group.sh 17 $idx
done
