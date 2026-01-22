#setwd("/public4/group_crf/home/cuirf/Rhinogobius/dp_scan");
library(data.table)
dat <- fread("normdp.wilcox.gz", header=T, sep="\t", showProgress=T, nThread=10, tmpdir="./tmp")
quantile(dat$WilcoxP, c(1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.025, 0.975, 0.9999), na.rm=T)
# 0.0000001%    0.000001%     0.00001%      0.0001%       0.001%        0.01%         0.1%           1%         2.5% 
# 4.960180e-08 6.989645e-08 4.810676e-07 3.273075e-06 3.288505e-05 2.311158e-04 1.585376e-03 1.165259e-02 2.674284e-02 
# 97.5%       99.99% 
#   9.782992e-01 1.000000e+00 

which(dat$WilcoxP == 4.960180e-08);
datChr <- dat[dat$CHROM=="Chromosome1",];
datChr$LogP <- -log10(datChr$WilcoxP)

plot(datChr$START , datChr$LogP, type = "l")
