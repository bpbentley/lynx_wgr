source('~/R_functions.r')
Ne<- get_Ne("./SSLR")
write.table(Ne$Ne_est, "./SSLR.Ne", sep = '\t')

#bsub -q long -R rusage[mem=60000] -n 4 -W 12:00 -R span\[hosts=1\] "Rscript ./NFLD_estimate_Ne.r"
