##############

#setwd("/mnt/Timina/dortega/hlopezh/data/Alessandro/")
library(deSolve)
library(rootSolve)

### Let's read all the stuff
Medians <- read.csv("median.csv")
args = commandArgs(trailingOnly=TRUE)

print (args)

SNPFrequencyFileNumber = args[1]
SNPFrequencyFileTxt = paste("SNP/SNPData",SNPFrequencyFileNumber, ".txt", sep = "")
SNPFrequencyFile = read.table(SNPFrequencyFileTxt)

head(SNPFrequencyFile)
head(Medians)
nrow(SNPFrequencyFile)
nrow(Medians)

## Now let's reorder the data a little bit
SNPFrequencyFile <- SNPFrequencyFile[SNPFrequencyFile$V3 != "no",]

colnames(SNPFrequencyFile) <- c('Chr','Pos','Cluster','A1', 'A2','AlleleFrequency','DerivedAlleles','AlleleCount')

nrow(SNPFrequencyFile)

MergedDataset <- merge(Medians, SNPFrequencyFile)

MergedDataset<-MergedDataset[MergedDataset$Lat>-25 & MergedDataset$Lat<63,]
SortedDatasetByAge <- MergedDataset[order(MergedDataset$Date_mean),]

### Get the tentative origin of the allele

LastOcurrenceData <- 0

for (i in 1:nrow(SortedDatasetByAge)){
    if (SortedDatasetByAge$DerivedAlleles[i] > 0){
        LastOcurrenceData <- i
    }
}

AlleleOriginLat <- SortedDatasetByAge$Lat[LastOcurrenceData]
AlleleOriginLong <- SortedDatasetByAge$Long[LastOcurrenceData]
AlleleOriginAge <- SortedDatasetByAge$Date_mean[LastOcurrenceData]

######### Change coordinates

AlleleXAxis <- round(((SortedDatasetByAge$Lat - min(SortedDatasetByAge$Lat))/(max(SortedDatasetByAge$Lat)- min(SortedDatasetByAge$Lat))) * 99 + 1)
AlleleYAxis <- round(((SortedDatasetByAge$Long - min(SortedDatasetByAge$Long))/(max(SortedDatasetByAge$Long)- min(SortedDatasetByAge$Long))) * 99 + 1)

SortedDatasetByAge <- cbind(SortedDatasetByAge,AlleleXAxis)
SortedDatasetByAge <- cbind(SortedDatasetByAge,AlleleYAxis)

AlleleOriginLatX <- round(((AlleleOriginLat - min(SortedDatasetByAge$Lat))/(max(SortedDatasetByAge$Lat)- min(SortedDatasetByAge$Lat))) * 99 + 1)
AlleleOriginLongY <- round(((AlleleOriginLong - min(SortedDatasetByAge$Long))/(max(SortedDatasetByAge$Long)- min(SortedDatasetByAge$Long))) * 99 + 1)


######### Check this reference to see how the equation is solved https://cran.r-project.org/web/packages/rootSolve/vignettes/rootSolve.pdf

diffusion2D <- function(t, conc, par) {
Conc <- matrix(nrow = n, ncol = n, data = conc) # vector to 2-D matrix
dConc <- Conc*(1-Conc)*(Conc*d+s*(1-2*Conc)) # consumption
#dConc <- -r*Conc*Conc
BND <- rep(1, n) # boundary concentration
# constant production in certain cells
# dConc[ii]<- dConc[ii]+ p
#diffusion in X-direction; boundaries=imposed concentration

Flux <- -Dx * rbind(rep(0, n), (Conc[2:n,]-Conc[1:(n-1),]),
rep(0, n) )/dx
dConc <- dConc - (Flux[2:(n+1),] - Flux[1:n,])/dx
#diffusion in Y-direction

Flux <- -Dy * cbind(rep(0, n), (Conc[,2:n]-Conc[,1:(n-1)]),
rep(0, n))/dy
dConc <- dConc - (Flux[,2:(n+1)]-Flux[,1:n])/dy
return(list(as.vector(dConc)))
}

# parameters
dy <- dx <- 1 # grid size
Dy <- Dx <- 1.5 # diffusion coeff, X- and Y-direction
n <- 100
d <- 0.0
s <- 0.0
pars<-c(n,d,s)
# 10 random cells where substance is produced at rate p
# ii <- trunc(cbind(runif(10)*n+1, runif(10)*n+1))
### LL calculations

DifussionValuesToCheck <- c(10)
Counter <- 0

TimesToTest <- c()
Conc0 <- matrix(nrow = n, ncol = n, 0.)
for (i in 1:nrow(SortedDatasetByAge)){
    Time = round((AlleleOriginAge - SortedDatasetByAge$Date_mean[i])/25)
    if (Time >= 0){
    TimesToTest <- c(TimesToTest, Time)
    if (Time == 0){
        Conc0[SortedDatasetByAge$Lat[i], SortedDatasetByAge$Long[i]] = round(SortedDatasetByAge$AlleleFrequency[i],3)
    }
    }
}

SelectionValuesToCheck <- c(-0.5,-0.4,-0.3,-0.2,-0.1,0)
if (TimesToTest[1] > 400){
LLFile <- paste("SNP/LLSel", SNPFrequencyFileNumber, ".txt", sep = "")
LLToPrint <- c(0,0,0,0,0)
LLToPrint <- t(LLToPrint)
write.table(LLToPrint,file=LLFile,append=TRUE,row.names=FALSE,col.names=FALSE)

}else{
LL<-c()
for (j in 1:6){
Dy <- Dx <- 1
d = SelectionValuesToCheck[j]
s = d / 2
LL<-c(LL,0)
N = 1000 ### Population size per deme
    print (i)

    print(system.time(
    ST3 <- ode.2D(y = Conc0, times = c(0,TimesToTest), func = diffusion2D, parms = pars,
    dimens = c(n, n), method = rkMethod("rk45ck"))
    ))

    for (i in 1:length(TimesToTest)){
        if ( (ST3[i+1,(SortedDatasetByAge$AlleleXAxis[i]-1)*100 + SortedDatasetByAge$AlleleYAxis[i]] > 0) && (SortedDatasetByAge$AlleleCount[i] > 0) && (ST3[i+1,(SortedDatasetByAge$AlleleXAxis[i]-1)*100 + SortedDatasetByAge$AlleleYAxis[$        # print (i, log (dbinom(SortedDatasetByAge$DerivedAlleles[i],SortedDatasetByAge$AlleleCount[i], ST3[i+1,(SortedDatasetByAge$AlleleXAxis[i]-1)*100 + SortedDatasetByAge$AlleleYAxis[i]])))
    LL[j] <- LL[j] + log (dbinom(SortedDatasetByAge$DerivedAlleles[i],SortedDatasetByAge$AlleleCount[i], ST3[i+1,(SortedDatasetByAge$AlleleXAxis[i]-1)*100 + SortedDatasetByAge$AlleleYAxis[i]]))
    }}
}
LLToPrint <- t(LL)

LLFile <- paste("SNP/LLSel", SNPFrequencyFileNumber, ".txt", sep = "")
write.table(LLToPrint,file=LLFile,append=TRUE,row.names=FALSE,col.names=FALSE)
}