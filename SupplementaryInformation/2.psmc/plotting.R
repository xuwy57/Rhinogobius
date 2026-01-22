#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/PLPv1.2LG/msmc2");
#pdf("allpops.pdf", width=7, height=6);
#pdf("praslin.pdf", width=7, height=6);
#pdf("mahesouth.pdf", width=7, height=6);
pdf("allwild.pdf", width=7, height=6);

sMainFolder <- "msmc2ret/";
sBSFolder <- "bootstrapped/"
nBSReps <- 30;
arrYLim <- c(1e4, 5e7);
arrXLim <- c(5e4, 5e6);

nMu <- 7.75273E-9; # estimated from divergence between BTP and Miichthys, this is per generation. assuming 10 years of generation time, per year is 7.75273E-10


fnPlotPop <- function(sPop, arrCol, bBSPlot = F, bAppendPlot=T, sThisBSFile1 = "") {
 
  sInFile1 <- paste(sMainFolder,"/",sPop, ".final.txt", sep="");
  arrCol <- add.alpha(arrCol, 0.8);
  
  if (bBSPlot) {
    cat("bs: ", sThisBSFile1,"\n");
    sInFile1 <- sThisBSFile1;
    arrCol <- add.alpha(arrCol, 0.2);
  }
  
  cat("Open ", sInFile1,"\n");
  datMSMC1 <- read.table(sInFile1, header=T, sep="\t" );
  datMSMC1$gen <- datMSMC1$left_time_boundary/nMu;
  datMSMC1$gen[1] <- 0.01;
  datMSMC1$popsize <- (1/datMSMC1$lambda)/(2*nMu);
  
  
  if (bBSPlot == F && (!bAppendPlot) ) {
    plot( datMSMC1$gen ,datMSMC1$popsize ,  type = 's' , log='xy', xlim=arrXLim, ylim = arrYLim, xlab="Generations" , ylab = "Pop size", lwd=3, col=arrCol[1])
  } else {
    lines( datMSMC1$gen ,datMSMC1$popsize ,  type = 's' , lwd=3, col=arrCol[1])
    
  }
  

  if (bBSPlot==TRUE) {
    #plot variation
    return();
  } else {
    cat("Plot bs...\n");
    for (nRep in 1:nBSReps) {
      sBSF1 <- paste(sBSFolder,'/',sPop, '/_', nRep, '/out.final.txt' , sep="");
      #cat(sBSF1);
      if (file.exists(sBSF1)) {
        fnPlotPop(sPop, arrCol, bBSPlot = T, sThisBSFile1=sBSF1 )
      }
    }
    
    # datReal <- read.table( paste(sRealFolder, '/',sRealFile ,sep=""), header=T, sep="\t");
    # points(datReal$Pop1_Gen, datReal$Pop1_Popsize, col=arrCol[2])
    # points(datReal$Pop2_Gen, datReal$Pop2_Popsize, col=arrCol[3])
    # 
  }
  
  #construct output table:
  return();
  
}

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#Wild Praslin:

datSamples <- read.table("samples.txt", sep=' ', header=F);
arrSamples <- datSamples$V1;
arrCol <- rainbow(length(arrSamples));
bAppend <-F;
for(i in 1:length(arrSamples) ) {
	fnPlotPop(arrSamples[i], arrCol[i], bAppendPlot = bAppend);
	bAppend <-T
}

dev.off();

