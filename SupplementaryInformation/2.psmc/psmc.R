setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/plot/psmc");
library(ape)
library(geiger)
library(picante)
library(phytools)
library(brms)
library(BH)
library(ggplot2)
library(reshape2)
#缺z9
sMainFolder <- "./psmc_result/";
datSamples <- read.table("samples.txt", sep=' ', header=F);
arrSamples <- datSamples$V1;
samples_order <- read.table("samples.txt", header = FALSE, stringsAsFactors = FALSE)$V1
samples_order <- gsub(".bam", "", samples_order)
time_points <- seq(1e+05, 1e+06, length.out = 10)
# 初始化存储结果的数据框
all_species_data <- data.frame()
# 处理每个物种的数据
nMu <- 2.6440E-9 #4fold枝长0.1203，mcmctree分化时间0.455*100Mya, per year is 2.6440E-9,per generation per year
for (sPop in arrSamples) {
  sInFile1 <- paste(sMainFolder, "/", sPop, ".final.txt", sep = "")
  cat("Processing:", sInFile1, "\n")
  # 读取种群大小数据
  datMSMC1 <- read.table(sInFile1, header = TRUE, sep = "\t")
  datMSMC1$gen <- datMSMC1$left_time_boundary / nMu  # 计算世代数
  datMSMC1$popsize <- (1 / datMSMC1$lambda) / (2 * nMu)  # 计算种群大小
  # 过滤时间范围
  filtered_data <- datMSMC1[datMSMC1$gen > 1e+05 & datMSMC1$gen < 1e+06, ]
  filtered_data$Species <- gsub(".bam", "", sPop)  # 添加物种名称
  # 合并数据
  all_species_data <- rbind(all_species_data, filtered_data)
}

# 确保数据格式正确
all_species_data$gen <- as.numeric(all_species_data$gen)
all_species_data$popsize <- as.numeric(all_species_data$popsize)


# 对每个物种，找到最近的时间点的种群大小
simplified_data <- all_species_data %>%
  group_by(Species) %>%
  summarise(
    closest_gen = list(sapply(time_points, function(t) gen[which.min(abs(gen - t))])),
    closest_popsize = list(sapply(time_points, function(t) popsize[which.min(abs(gen - t))]))
  ) %>%
  unnest(cols = c(closest_gen, closest_popsize))

#重新整理数据，以便用于热图
heatmap_data <- data.frame(
  Species = rep(unique(simplified_data$Species), each = length(time_points)),
  gen = rep(time_points, times = 50),
  popsize = unlist(simplified_data$closest_popsize)
)
heatmap_data$Species <- factor(heatmap_data$Species, levels = samples_order)
#绘制热图
ggplot(heatmap_data, aes(x = factor(gen), y = Species, fill = popsize)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Species Population Size Over Time",
       x = "Time (Generations)",
       y = "Species",
       fill = "Population Size") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))









#计算种群大小变化斜率
fnCalPop <- function(sPop, i) {
  sInFile1 <- paste(sMainFolder,"/",sPop, ".final.txt", sep="");
  cat("Open ", sInFile1,"\n");
  datMSMC1 <- read.table(sInFile1, header=T, sep="\t" );
  datMSMC1$gen <- datMSMC1$left_time_boundary/nMu;
  datMSMC1$gen[1] <- 0.01;
  datMSMC1$popsize <- (1/datMSMC1$lambda)/(2*nMu);
  filtered_data <- datMSMC1[datMSMC1$gen > 1e+05 & datMSMC1$gen < 1e+06, ]
  model <- lm(filtered_data$popsize ~ filtered_data$gen)
  summary(model)
  slope <- coef(model)[2]
  cat("斜率为:", slope, "\n")
  
  time_points <- seq(1e+05, 1e+06, length.out = 10)
  pop_sizes <- sapply(time_points, function(t) {
    approx(filtered_data$gen, filtered_data$popsize, xout = t)$y
  })
  
  return(c(gsub(".bam", "", sPop), -slope, pop_sizes))
}

file_list <- list.files(sMainFolder)
results_df <- data.frame(Species = character(length(file_list)), slope = numeric(length(file_list)), stringsAsFactors = FALSE)
colnames(results_df) <- c("Species", "slope")

for(i in 1:length(arrSamples) ) {
  results_df[i, ] <- fnCalPop(arrSamples[i], i);
}
print(results_df)
write.table(results_df, file = "slope_traits.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

 
#phygenic analysis
oTree <- read.tree("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/cafe5/cafe/tree.tre")
plot(oTree) #Plot the phylogeny
datTraits <- read.table("slope_traits.tsv", header = TRUE, sep = "\t")
datTraits$lifestyle <- as.numeric(datTraits$lifestyle)
datTraits$slope <- as.numeric(datTraits$slope)
arrSlope <- datTraits$slope
names(arrSlope) <- datTraits$Species
arrLifestyle <- datTraits$lifestyle
names(arrLifestyle) <- datTraits$Species

fnFitAllEvolModels <- function(oTree, arrMeans) {
  oTreeLambda0 <- rescale(oTree, model="lambda", 0); #set lambda to 0 so that the tree becomes a star tree (no phylogenetic signal)
  lambdaModel <- fitContinuous(oTree, arrMeans, model = "lambda", SE=NA)
  brownianModel <- fitContinuous(oTree, arrMeans, model = "BM", SE=NA)
  nosigModel <- fitContinuous(oTreeLambda0, arrMeans, SE=NA)
  ouModel <- fitContinuous(oTree, arrMeans, model = "OU", SE=NA)
  ebModel <- fitContinuous(oTree, arrMeans, model = "EB", SE=NA)
  kappaModel <- fitContinuous(oTree, arrMeans, model = "kappa", SE=NA)
  deltaModel <- fitContinuous(oTree, arrMeans, model = "delta", SE=NA)
  ratetrendModel <- fitContinuous(oTree, arrMeans, model = "rate_trend", SE=NA)
  arrAICc <- c(lambdaModel$opt$aicc, brownianModel$opt$aicc, nosigModel$opt$aicc, ouModel$opt$aicc, 
               ebModel$opt$aicc, kappaModel$opt$aicc, deltaModel$opt$aicc, ratetrendModel$opt$aicc);
  names(arrAICc) <- c("lambda", "brownian", "nosig", "OU", "EB", "kappa", "delta", "rate_trend");
  return(arrAICc)
}
fnModelWeights <- function(arrAICc) {
  aiccD <- arrAICc - min(arrAICc);
  aw <- exp(-0.5 * aiccD)
  return(aw/sum(aw))
}

#Evolutionary modelling of single traits
#Bloomberg's K>1 suggests that closely related species are more similar than expected by Brownnian motion (drift), 
#K<1 suggests that related species are more dissimilar compared to BM
phylosignal(arrSlope, oTree)
phylosig(oTree, arrSlope, method = "K", test = T); #Bloomberg's K
phylosignal(arrLifestyle, oTree)
phylosig(oTree, arrLifestyle, method = "K", test = T); #Bloomberg's K

#Pagel's Lambda, lambda = 0 suggest no phylogenetic signal, lambda = 1 suggest trait evolves completely according to the phylogeny following the BM model
phylosig(oTree, arrSlope, method = "lambda", test = T) #Pagel's lambda
phylosig(oTree, arrLifestyle, method = "lambda", test = T) #Pagel's lambda
oTreeLambda0 <- rescale(oTree, model="lambda", 0); #set lambda to 0 so that the tree becomes a star tree (no phylogenetic signal)
plot(oTreeLambda0)

#Compare the Akaike information criterion (corrected) AICc
#Which is the best model? (AICc smallest?)
arrAICc <- fnFitAllEvolModels(oTree, arrSlope)
which(arrAICc==min(arrAICc)); #what is the best model?
fnModelWeights(arrAICc); #relative model weight

##### bug
#Analyze the correlation between two traits on a phylogeny
#Now, perform Felsenstein's phylogenetically independent contrast on the traits
pic.Slope <- pic(arrSlope, oTree);
pic.Lifestyle <- pic(arrLifestyle, oTree);

plot(pic.Lifestyle, pic.Slope)
summary(lm(pic.Slope ~ pic.Lifestyle))
#How did the p value change, compared to not correcting for the phylogeny?


##Use bayesian method to estimate correlated trait evolution:
oPhyloCovar <- ape::vcv.phylo(oTree)
View(as.data.frame(oPhyloCovar))
model_brms <- brm(
  slope ~ lifestyle + (1|gr(Species, cov = oPhyloCovar)),
  data = datTraits,
  family = gaussian(),
  data2 = list(oPhyloCovar = oPhyloCovar),
  iter = 2000000,
  thin = 100
)
summary(model_brms)
plot(model_brms, nvariables=2, ask = FALSE)
plot(conditional_effects(model_brms), points = TRUE); #Plot the phylogenetic corrected correlation
bayes_lambda <- "sd_Species__Intercept^2 / (sd_Species__Intercept^2 + sigma^2) = 0"
bayes_lambda <- hypothesis(model_brms, bayes_lambda, class = NULL)
bayes_lambda #this is the lambda estimated using the bayesian method




