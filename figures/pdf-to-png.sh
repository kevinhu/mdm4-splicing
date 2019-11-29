#!/bin/sh

#convert PDFs to PNGs for preview in Google docs

for i in *.pdf; do
    
  name=$i;
  name=${name%.*};
  sips -s format png $i --out ${name}.png;

done