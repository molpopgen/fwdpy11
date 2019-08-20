#!/usr/bin/env Rscript

library(dplyr)
library(DBI)

args <- commandArgs(TRUE)
dbfile <- args[1]
xmin <- as.numeric(args[2])
xmax <- as.numeric(args[3])

db = dbConnect(RSQLite::SQLite(), dbfile)
dbt = tbl(db, 'data')

x = collect(dbt)

pdf("stats.pdf")
par(mfrow=c(3,1))
plot(x$vg,type='l',xlab="", ylab="Genetic variance",xlim=c(xmin,xmax))
plot(x$zbar,type='l',xlab="", ylab="Mean trait value",xlim=c(xmin,xmax))
plot(x$wbar,type='l',xlab="Time (generations)", ylab="Mean fitness",xlim=c(xmin,xmax))
dev.off()


