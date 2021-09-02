for fn in `ls data/MR/MR_sumstats/meta/*.bgz`; do mv ${fn} ${fn//.bgz/.gz}; done
for fn in `ls data/MR/MR_sumstats/*.bgz`; do mv ${fn} ${fn//.bgz/.gz}; done
for fn in `ls data/MR/MR_sumstats/*mtag_meta.txt`; do gzip ${fn}; done
