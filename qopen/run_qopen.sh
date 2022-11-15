
qopen go -e 38474119 --prefix "00_test/"
qopen go --prefix "01_go/" --overwrite-conf '{"skip": {"coda_window": 5, "num_pairs": 5}}'  # first run of Qopen with 55 events
qopen fixed --input "01_go/results.json" --prefix "02_sites/" --no-plot-results --align-sites --align-sites-station "CI.MPM,CI.WBM,CI.WMF"  # rerun and fix attenuation
qopen source --events "../data/cat55+1.csz" --input "01_go/results.json" --input-sites "02_sites/results.json" --prefix "03_source/" --no-plot-results  # rerun and fix attenuation and site amplification, use 55+1 catalog
qopen recalc_source --events "../data/cat55+1.csz" --input "03_source/results.json" --prefix "04_source_nconst/" --seismic-moment-options '{"fc": null, "n": 2.58, "gamma": 2}'  # refit source spectra using constant n
