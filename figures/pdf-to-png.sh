#!/bin/sh

#convert PDFs to PNGs for preview in Google docs

for i in *.pdf; do

    echo $i
    
    name=$i;
    name=${name%.*};
    convert -density 512 -quality 100 -colorspace RGB $i ${name}.png;

done