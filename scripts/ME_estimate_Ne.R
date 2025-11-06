source('~/R_functions.r')
Ne<- get_Ne("./ME")
write.table(Ne$Ne_est, "./ME.Ne", sep = '\t')
