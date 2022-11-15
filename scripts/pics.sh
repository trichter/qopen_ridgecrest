#!/usr/bin/env bash
cp ../qopen/02_sites/plots/sites.pdf ../figs/
cp ../qopen/04_source_nconst/plots/sds.pdf ../figs/

for f in Q l sites sds8 eqparams mags2; do
  convert -density 300x300 ../figs/$f.pdf ../figs/png/$f.png
done

