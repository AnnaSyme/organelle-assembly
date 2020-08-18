#!/usr/bin/env bash

#convert gfa files into png images
#using Bandage image
#formatting set so that image is not cropped off page

#in: ~/Bandage 
 
```
for file in ~/gfas-directory/*;\
do Bandage.app/Contents/MacOS/Bandage image $file ${file}.png \
--lengths --depth --fontsize 6 --edgelen 50 --toutline 2 --minnodlen 50 --colour blastsolid; \
done
```


