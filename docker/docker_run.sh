#!/bin/bash

docker run \
--platform linux/amd64 \
--rm \
--name bioc_3_15 \
-it \
-d \
-e USER=nshen7 \
-e PASSWORD=Sn199789. \
-p 8787:8787 \
-v /Users/sn/Volumes/sockeye_scratch:/src \
nshen7/bioc_3_15
