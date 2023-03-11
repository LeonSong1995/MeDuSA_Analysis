#! /usr/lib/R/bin/Rscript
library(peer)
library(argparser, quietly=TRUE)
library(data.table)

p = arg_parser("Prepare Gene expression")
p = add_argument(p, "exp.file", help="")
p = add_argument(p, "--output_dir", short="-o", help="Output directory", default=".")
argv = parse_args(p)


RINT = as.data.frame(fread(argv$exp.file))
RINT = RINT[,-seq(1,4)]


if(ncol(RINT)<150){NK = 10}
if(ncol(RINT)<250 & ncol(RINT)>=150){NK = 30}
if(ncol(RINT)<350 & ncol(RINT)>=250){NK = 45}
if(ncol(RINT)>=350){NK = 60}

expr = t(RINT)
model = PEER()
PEER_setPhenoMean(model,as.matrix(expr))
PEER_setNmax_iterations(model,1e+4)
PEER_setNk(model,NK)
PEER_update(model)
PEER_factor = PEER_getX(model)
colnames(PEER_factor) = paste0('PEER',seq(1,ncol(PEER_factor)))
rownames(PEER_factor) = colnames(RINT)

write.table(PEER_factor,argv$output_dir,row.names=T)