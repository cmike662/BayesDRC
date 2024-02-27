
library(tidyverse)
library(ecotox)
library(rstan)
library(ggplot2)
library(gt)
library(rethinking)


#set.seed(123)
nPops =1
useLogit = FALSE
plottrace = F  #traceplots can slow down pdf viewing
currentResultsRow = 0
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

Results1 <- data.frame(Pop = character(),
                       Model = character(),
                       adjMort = character(),
                       OD = character(),
                       Pooled = character(),
                       group = numeric(),
                       LC50 = numeric(),
                       LCmin = numeric(),
                       LCmax = numeric(),
                       Slope = numeric(),
                       Slopemin = numeric(),
                       Slopemax = numeric(),
                       nSurv = numeric(),
                       nSurvmin = numeric(),
                       nSurvmax = numeric(),
                       kappa = numeric(),
                       kappamin = numeric(),
                       kappamax = numeric()
)


library(readxl)
fileName = file.choose()
namedPops = sort(excel_sheets(fileName))
nPops = length(namedPops)
Data <- data.frame()
plottrace = F  #traceplots can slow down pdf viewing
currentResultsRow = 0
ResultsOD <- setNames(data.frame(matrix(ncol = length(namedPops), nrow = 0)), c("Pop",namedPops[-1]))

for (k1 in 1:(length(namedPops)-1)){
  ResultsOD[k1,1] = namedPops[k1]
}
ResultsnoOD <- ResultsOD

for (j1 in 1:(length(namedPops))){
  namePop1 = namedPops[j1]
  Data1 <- read_excel(fileName, sheet=namePop1)
  nBlocks = 1
  Data1$Population <- j1
  Data1$block <- 1
  Data1$surv = Data1$tested - Data1$dead
  Data1$SRate <- Data1$surv/Data1$tested
  Data <- rbind(Data, Data1)
}


## ----Run Stan model} #, warning=FALSE, message=FALSE---------------------------------------
DataUncorrected <- Data
DataUncorrectedNoCntrl <- Data[Data$Dose!=0,] #Just to see what happens when
#we remove control doses
DataUncorrected$scaledDose = DataUncorrected$Dose/max(DataUncorrected$Dose)
sdata <- list(
  N = nrow(DataUncorrected),
  nPops = nPops,
  nBlocks = nBlocks,
  scaledDose = DataUncorrected$scaledDose,
  block = DataUncorrected$block,
  tested = DataUncorrected$tested,
  surv = DataUncorrected$surv,
  Population = DataUncorrected$Population,
  correction = max(DataUncorrected$Dose),
  uselogit = useLogit
)


Model43OD <- stan(
  file = "Model43v2-OD-reparamV2-probit.stan",
  data = sdata,    # named list of data
  chains = 6,             # number of Markov chains
  iter = 6000,            # total number of iterations per chain
  warmup = 1000,
  cores = 6,              # number of cores (could use one per chain)
  control = list(adapt_delta = 0.99, stepsize = 0.1, max_treedepth = 15),
  refresh = 500             # no progress shown
)

precis(Model43OD, depth=2, pars=c("lc50", "B", "nSurv", "theta"))
precis(Model43OD)
postOD <- extract.samples(Model43OD, n=20000, pars=c("lc50", "B", "nSurv", "theta"))
precis(postOD, depth=2, pars=c("lc50", "B", "nSurv", "theta"))



Model43NoOD <- stan(
  file = "Model43v2-noOD-Complete-probit.stan",
  data = sdata,    # named list of data
  chains = 6,             # number of Markov chains
  iter = 6000,            # total number of iterations per chain
  warmup=1000,
  cores = 6,              # number of cores (could use one per chain)
  refresh = 500 ,            # no progress shown
  control = list(adapt_delta = 0.99, stepsize = 0.1, max_treedepth = 15)
)

precis(Model43NoOD, depth=2, pars=c("lc50", "B", "nSurv"))#, "theta"))
precis(Model43NoOD)
postnoOD <- extract.samples(Model43NoOD, n=20000, pars=c("lc50", "B", "nSurv"))#, "theta"))
precis(postnoOD, depth=2, pars=c("lc50", "B", "nSurv"))#, "theta"))

Model43ODNoPool <- stan(
  file = "Model43v2-Complete-NothingPooled-probit.stan",
  data = sdata,    # named list of data
  chains = 6,             # number of Markov chains
  iter = 6000,            # total number of iterations per chain
  warmup = 1000,
  cores = 6,              # number of cores (could use one per chain)
  refresh = 500 ,            # no progress shown
  control = list(adapt_delta = 0.99, stepsize = 0.1, max_treedepth = 15)
)

precis(Model43ODNoPool, depth=2, pars=c("lc50", "B", "nSurv", "theta"))
postODNoPool <- extract.samples(Model43ODNoPool, n=20000, pars=c("lc50", "B", "nSurv","theta"))
precis(postODNoPool, depth=2, pars=c("lc50", "B", "nSurv", "theta"))

print(rethinking::compare(Model43OD, Model43NoOD, Model43ODNoPool))


## ----Analysis of stan model, echo=F--------------------------------------------------------
cat("\n\n===================================================\n\n")
cat("Analysis of model with overdispersion\n")
print(mR <- precis(Model43OD,depth=2, prob=0.95,pars=c("lc50","B", "nSurv", "theta")))
#postOD <- extract.samples(Model43OD, n=20000, pars=c("lc50", "B", "nSurv", "theta"))

currentResultsRow = 0

for (j1 in 1:nPops){
  currentResultsRow = currentResultsRow + 1
  cat("HPDI lc50,",j1,":", HPDI(postOD$lc50[,j1], p=0.95),"\n")
  cat("Peak mode pop",j1, ":", chainmode(postOD$lc50[,j1]), "\n")
  
  {
    Results1[currentResultsRow,] = list(namedPops[j1],
                                        "Bayes1", "T", "T","T", namedPops[j1],
                                        chainmode(postOD$lc50[,j1]),
                                        HPDI(postOD$lc50[,j1], 0.95)[1],
                                        HPDI(postOD$lc50[,j1], 0.95)[2],
                                        chainmode(postOD$B[,j1]),
                                        HPDI(postOD$B[,j1], 0.95)[1],
                                        HPDI(postOD$B[,j1], 0.95)[2],
                                        chainmode(postOD$nSurv[,j1]),
                                        HPDI(postOD$nSurv[,j1], 0.95)[1],
                                        HPDI(postOD$nSurv[,j1], 0.95)[2],
                                        #chainmode(postOD$theta[,j1],   bw=0.5, n=2^16),
                                        chainmode(postOD$theta[,j1], n=2^16),
                                        HPDI(postOD$theta[,j1], 0.95)[1],
                                        HPDI(postOD$theta[,j1], 0.95)[2] 
    )  
  }
  
  ### Generate graphs
  xmax = max(HPDI(postOD$lc50[,j1], 0.95)[2]) * 1.2
  ymax = max(density(postOD$lc50[,j1])$y)*1.2  
  
  dens(postOD$lc50[,j1], ylab="Posterior plausible probability density",
       main = "Model with overdispersion and pooling",
       xlim = c(min(postOD$lc50), xmax), 
       ylim = c(0, ymax),
       xlab = parse(text=paste0('"LC"[50]')), lty=1, col="red", lwd=2)
  legend("topright", inset=0.05,  legend=c(namedPops[j1]),
         col=c("red", "blue"), lty=1:2, lwd=2, cex=0.8)  
  
  
  
  
  xmax = max(HPDI(postOD$B[,j1], 0.95)[2]) * 1.2
  ymax = max(density(postOD$B[,j1])$y)*1.2  
  
  dens(postOD$B[,j1], ylab="Posterior plausible probability density",
       main = "Model with overdispersion and pooling",
       xlim = c(min(postOD$B),xmax),
       ylim = c(0, ymax),
       xlab=bquote("slope"),lty=1, col="red", lwd=2)
  
  legend("topright", inset=0.05,  legend=c(namedPops[j1]),
         col=c("red", "blue"), lty=1:2, lwd=2, cex=0.8)  
  
  
  
  
  xmax = max(HPDI(postOD$nSurv[,j1], 0.95)[2]) * 1.2
  ymax = max(density(postOD$nSurv[,j1])$y)*1.2  
  
  dens(postOD$nSurv[,j1], ylab="Posterior plausible probability density",
       main = "Model with overdispersion and pooling",
       xlim = c(0,1),  
       ylim = c(0, ymax),   
       xlab=bquote("Natural survivorship"),lty=1, col="red", lwd=2)
  legend("topleft", inset=0.05,  legend=c(namedPops[j1]),
         col=c("red", "blue"), lty=1:2, lwd=2, cex=0.8) 
  
  
  
  
  xmax = (HPDI(postOD$theta[,j1], 0.95)[2]) * 1.2
  ymax = max(density(postOD$theta[,j1])$y)*1.2  
  
  plot(density(postOD$theta[,j1], n=2^16),xlim=c(2,xmax),  
       main = "Model with overdispersion and pooling",
       ylab="Posterior plausible probability density",
       xlab=bquote("kappa"),lty=1, col="red", lwd=2)
  
  legend("topright", inset=0.05,  legend=c(namedPops[j1]),
         col=c("red", "blue"), lty=1:2, lwd=2, cex=0.8) 
  
  
  
  #if (plottrace)
  #  print(traceplot(Model43, pars=c("a","B", "nSurv", "theta")))
  #print(trankplot(Model43, pars=c("a","B", "nSurv", "theta")))
  par(mfrow=c(1,1))
}


## ----Analysis of stan model no OD, echo=F--------------------------------------------------
cat("\n\n===================================================\n\n")
cat("Analysis of model without overdispersion\n")
print(mR <- precis(Model43NoOD,depth=2, prob=0.95,pars=c("lc50","B", "nSurv") ))

#postnoOD <- extract.samples(Model43NoOD, n=20000)



for (j1 in 1:nPops){
  currentResultsRow = currentResultsRow + 1
  cat("HPDI lc50,",j1,":", HPDI(postnoOD$lc50[,j1], p=0.95),"\n")
  cat("Peak mode pop",j1, ":", chainmode(postnoOD$lc50[,j1]), "\n")
  
  Results1[currentResultsRow,] = list(namedPops[j1],
                                      "Bayes2", "T", "F","T",  namedPops[j1],
                                      chainmode(postnoOD$lc50[,j1]),
                                      HPDI(postnoOD$lc50[,j1], 0.95)[1],
                                      HPDI(postnoOD$lc50[,j1], 0.95)[2],
                                      chainmode(postnoOD$B[,j1]),
                                      HPDI(postnoOD$B[,j1], 0.95)[1],
                                      HPDI(postnoOD$B[,j1], 0.95)[2],
                                      chainmode(postnoOD$nSurv[,j1]),
                                      HPDI(postnoOD$nSurv[,j1], 0.95)[1],
                                      HPDI(postnoOD$nSurv[,j1], 0.95)[2],
                                      NA,NA,NA
  )  
  
  xmax = max(HPDI(postnoOD$lc50[,j1], 0.95)[2]) * 1.2
  ymax = max(density(postnoOD$lc50[,j1])$y)*1.2  
  
  dens(postnoOD$lc50[,j1], ylab="Posterior plausible probability density",
       main = "Model without overdispersion",
       xlim = c(min(postnoOD$lc50), xmax), 
       ylim = c(0, ymax),
       xlab = parse(text=paste0('"LC"[50]')), lty=1, col="red", lwd=2)
  legend("topright", inset=0.05,  legend=c(namedPops[j1]),
         col=c("red", "blue"), lty=1:2, lwd=2, cex=0.8)  
  
  xmax = max(HPDI(postnoOD$B[,j1], 0.95)[2]) * 1.2
  ymax = max(density(postnoOD$B[,j1])$y)*1.2  
  
  dens(postnoOD$B[,j1], ylab="Posterior plausible probability density",
       main = "Model without overdispersion",
       xlim = c(min(postnoOD$B),xmax),
       ylim = c(0, ymax),
       xlab=bquote("slope"),lty=1, col="red", lwd=2)
  
  legend("topright", inset=0.05,  legend=c(namedPops[j1]),
         col=c("red", "blue"), lty=1:2, lwd=2, cex=0.8)  
  
  xmax = max(HPDI(postnoOD$nSurv[,j1], 0.95)[2]) * 1.2
  ymax = max(density(postnoOD$nSurv[,j1])$y)*1.2  
  
  dens(postnoOD$nSurv[,j1], ylab="Posterior plausible probability density",
       main = "Model without overdispersion",
       xlim = c(0,1), #c(min(postnoOD$nSurv),xmax),   
       ylim = c(0, ymax),   
       xlab=bquote("Natural survivorship"),lty=1, col="red", lwd=2)
  legend("topleft", inset=0.05,  legend=c(namedPops[j1]),
         col=c("red", "blue"), lty=1:2, lwd=2, cex=0.8) 
  
  #if (plottrace)
  #  print(traceplot(Model43, pars=c("a","B", "nSurv")))
  #print(trankplot(Model43, pars=c("a","B", "nSurv")))
  par(mfrow=c(1,1))
}

## ----Analysis of stan model with OD but no pooling, echo=F--------------------------------------------------
cat("\n\n===================================================\n\n")
cat("Analysis of model with overdispersion but no pooling\n")
print(mR <- precis(Model43ODNoPool,depth=2, prob=0.95,pars=c("lc50","B", "nSurv", "theta") ))
#postODNoPool <- extract.samples(Model43ODNoPool, n=20000)

for (j1 in 1:nPops){
  currentResultsRow = currentResultsRow + 1
  cat("HPDI lc50,",j1,":", HPDI(postODNoPool$lc50[,j1], p=0.95),"\n")
  cat("Peak mode pop",j1, ":", chainmode(postODNoPool$lc50[,j1]), "\n")
  
  Results1[currentResultsRow,] = list(namedPops[j1],
                                      "Bayes3", "T", "T","F",  namedPops[j1],
                                      chainmode(postODNoPool$lc50[,j1]),
                                      HPDI(postODNoPool$lc50[,j1], 0.95)[1],
                                      HPDI(postODNoPool$lc50[,j1], 0.95)[2],
                                      chainmode(postODNoPool$B[,j1]),
                                      HPDI(postODNoPool$B[,j1], 0.95)[1],
                                      HPDI(postODNoPool$B[,j1], 0.95)[2],
                                      chainmode(postODNoPool$nSurv[,j1]),
                                      HPDI(postODNoPool$nSurv[,j1], 0.95)[1],
                                      HPDI(postODNoPool$nSurv[,j1], 0.95)[2],
                                      #chainmode(postODNoPool$theta[,j1]),
                                      chainmode(postODNoPool$theta[,j1],  n=2^16),
                                      HPDI(postODNoPool$theta[,j1], 0.95)[1],
                                      HPDI(postODNoPool$theta[,j1], 0.95)[2] 
  )  
  
  xmax = max(HPDI(postODNoPool$lc50[,j1], 0.95)[2]) * 1.2
  ymax = max(density(postODNoPool$lc50[,j1])$y)*1.2  
  
  dens(postODNoPool$lc50[,j1], ylab="Posterior plausible probability density",
       main = "Model with overdispersion but no pooling",
       xlim = c(0, xmax), 
       xlab = parse(text=paste0('"LC"[50]')), lty=1, col="red", lwd=2)
  legend("topright", inset=0.05,  legend=c(namedPops[j1]),
         col=c("red", "blue"), lty=1:2, lwd=2, cex=0.8)  
  
  xmax = max(HPDI(postODNoPool$B[,j1], 0.95)[2]) * 1.2
  ymax = max(density(postODNoPool$B[,j1])$y)*1.2  
  
  dens(postODNoPool$B[,j1], ylab="Posterior plausible probability density",
       main = "Model with overdispersion but no pooling",
       xlim = c(min(postODNoPool$B),xmax),
       ylim = c(0, ymax),
       xlab=bquote("slope"),lty=1, col="red", lwd=2)
  
  legend("topright", inset=0.05,  legend=c(namedPops[j1]),
         col=c("red", "blue"), lty=1:2, lwd=2, cex=0.8)  
  
  xmax = max(HPDI(postODNoPool$nSurv[,j1], 0.95)[2]) * 1.2
  ymax = max(density(postODNoPool$nSurv[,j1])$y)*1.2  
  
  dens(postODNoPool$nSurv[,j1], ylab="Posterior plausible probability density",
       main = "Model with overdispersion but no pooling",
       xlim = c(0,1), #c(min(postODNoPool$nSurv),xmax),   
       ylim = c(0, ymax),   
       xlab=bquote("Natural survivorship"),lty=1, col="red", lwd=2)
  legend("topleft", inset=0.05,  legend=c(namedPops[j1]),
         col=c("red", "blue"), lty=1:2, lwd=2, cex=0.8) 
  
  xmax = max(HPDI(postODNoPool$theta[,j1], 0.95)[2]) * 1.2
  ymax = max(density(postODNoPool$theta[,j1])$y)*1.2  
  
  plot(density(postODNoPool$theta[,j1],   n=2^16),xlim=c(2,xmax),  
       main = "Model with overdispersion but no pooling", 
       ylab="Posterior plausible probability density",
       xlab=bquote("kappa"),lty=1, col="red", lwd=2)
  
  legend("topright", inset=0.05,  legend=c(namedPops[j1]),
         col=c("red", "blue"), lty=1:2, lwd=2, cex=0.8) 
  
  #if (plottrace)
  #  print(traceplot(Model43, pars=c("a","B", "nSurv")))
  #print(trankplot(Model43, pars=c("a","B", "nSurv")))
  par(mfrow=c(1,1))
}

## ---- echo=F-------------------------------------------------------------------------------
library(BioRssay)
cat("\n\n===================================================\n\n")
cat("BioRssay analysis \n")
DataB <- Data
colnames(DataB) <- c("dose", "dead", "total","strain", "block", "SRate")
DataB$surv <- DataB$total - DataB$dead
DataT<-probit.trans(DataB, conf=0.05)
print(cm <-DataT$convrg)
head(DataT$tr.data)

dataBio<-DataT$tr.data
#dataBio contains a NAN, I think where it expects < 0% mort
dataBio <- na.omit(dataBio)
dataBio <- dataBio[dataBio$mort != 0,]
print(RR<-resist.ratio(dataBio ))
if (nPops > 1)
  print(model.signif(DataT$tr.data))
mort.plot(dataBio)
for (j1 in 1:nPops){
  currentResultsRow = currentResultsRow + 1
  Results1[currentResultsRow,] = list(namedPops[j1],
                                      "GLM-Bio", "T", "T","F", namedPops[j1],
                                      RR[j1, 13],
                                      RR[j1,14],
                                      RR[j1,15],
                                      RR[j1,1],
                                      RR[j1,1] - 1.96 * RR[j1,2],
                                      RR[j1,1] + 1.96 * RR[j1,2],
                                      #NA,NA,
                                      if (!is.null(DataT$convrg)){
                                      1-as.numeric(cm$ControlMortality[j1])}
                                      else { NA},
                                      NA,
                                      NA,
                                      NA,NA,NA
  )  
}



########Try ecotox####################### 
cSurv <- vector()
ABCorrected <- Data[-(1:nrow(Data)),]

#WHO gives correctedMort = (Mort-controlMort)/(1-controlMort)
#Another option is correctedMort = (controlSurv-Surv)/controlSurv
#Both of these give a corrected morality of 1.0 when there are no survivors
#BioRssay gives a corrected mortality of 0.9935 when there are no survivors
#The ABCorrected data uses a corrected mortality of 1 when there are no survivors

for (k1 in 1:nPops){
  #cSurv[k1] <- sum(as.numeric((Data[Data$Dose == 0 & Data$Population==k1,])['surv']))/
  #  sum(as.numeric((Data[Data$Dose == 0 & Data$Population==k1,])['tested']))
  cSurv[k1] <-sum((Data[Data$Dose==0 & Data$Population==k1,])['surv'])/sum((Data[Data$Dose==0 & Data$Population==k1,])['tested'])
}
cRow = 0

for (k1 in 1:nrow(Data)){
  if (Data$Dose[k1] > 0){
    cRow = cRow+1
    z = Data$surv[k1]/ (cSurv[Data$Population[k1]])
    z = ifelse(z < 0, 0,z)
    z = ifelse(z>Data$tested[k1], Data$tested[k1],z)
    ABCorrected[cRow, ] = Data[k1,]
    ABCorrected[cRow,'surv']=z
  }
}

for (j1 in 1:nPops){
  if (useLogit){
    conv <- LC_logit(cbind(tested-surv,surv) ~ log10(Dose), weights = tested, 
                    data = ABCorrected[ABCorrected$Population==j1,], p=50)
  } else {
    conv <- LC_probit(cbind(tested-surv,surv) ~ log10(Dose), weights = tested, 
                     data = ABCorrected[ABCorrected$Population==j1,], p=50)    
  }
  cat("Ecotox LC50:",conv$dose,"\n")
  cat("Ecotox Slope:", conv$slope, "\n")
  currentResultsRow = currentResultsRow + 1
  Results1[currentResultsRow,] = list(namedPops[j1],
                                      "GLM-Eco", "F", "F","F", namedPops[j1],
                                      conv$dose,
                                      conv$LCL,
                                      conv$UCL,
                                      conv$slope ,
                                      conv$slope - 1.96*conv$slope_se,
                                      conv$slope + 1.96*conv$slope_se,
                                      #NA,NA,
                                      cSurv[j1],
                                      NA,
                                      NA,
                                      NA,NA,NA
  ) 
}

colnames(Results1) <- c("Pop", "Model", "adjMort", "OD","Pooled", "group", "LC50", "LC50\nmin","LC50\nmax",
                        "slope", "slope\nmin","slope\nmax", "nSurv", "nsurv\nmin","nsurv\nmax",
                        "kappa", "kappa\nmin", "kappa\nmax")


SortedResults <- Results1[order(Results1$Pop),]
tab1 <- 
  SortedResults |>
  gt(
    groupname_col = "group") |>
  cols_hide(columns = c("Pop")) |>
  fmt_number(decimals = 3) |>
  fmt_number(columns=c("kappa\nmax"), decimals = 1) |>
  fmt_number(columns=c("kappa"), decimals = 1) |>
  fmt_number(columns=c("kappa\nmin"), decimals = 1) |>
  #page.orientation only works with .rtf output, fails here.  Perhaps will change someday?
  tab_options(
    page.orientation = "landscape"
    )  # |>

gtsave(tab1, "CompleteParamEstimates.docx" )


########Multiple comparison
if(length(namedPops) > 1){
  Rframe <- data.frame()
  RframenoOD <- data.frame()
  for (j1 in 1:(length(namedPops)-1)){
    namePop1 = namedPops[j1]
    for (j2 in (j1+1):length(namedPops)){
      namePop2 = namedPops[j2]
      postOD$diff <- (postOD$lc50[,j1] - postOD$lc50[,j2])
      ResultsOD[j1,j2] =  sum(postOD$diff > 0)/length(postOD$diff)
      postnoOD$diff <- (postnoOD$lc50[,j1] - postnoOD$lc50[,j2])
      ResultsnoOD[j1,j2] =  sum(postnoOD$diff > 0)/length(postnoOD$diff)
      Rframe <- rbind(Rframe, data.frame(max(ResultsOD[j1,j2], 1-ResultsOD[j1,j2]), j1, j2))
      RframenoOD <- rbind(RframenoOD, data.frame(max(ResultsnoOD[j1,j2], 1-ResultsnoOD[j1,j2]), j1, j2))
    }
  }
  colnames(Rframe) <-  c("p", "col", "row")
  colnames(RframenoOD) <-  c("p", "col", "row")
  
  tmp <- gt(ResultsOD)|>fmt_number(decimals = 3)
  tmp2 <- gt(ResultsnoOD)|>fmt_number(decimals = 3)
  
  sortedP <- order(Rframe$p, decreasing = T)
  sortednoODP <- order(RframenoOD$p, decreasing = T) 
  
  productP=1
  cat("Significance for OD data \n")
  for (j1 in 1:length(sortedP)){
    productP = productP * Rframe[sortedP[j1],1]
    if (productP < 0.95){break}
    cat(j1,": ", namedPops[Rframe[sortedP[j1],2]], " and ", 
        namedPops[Rframe[sortedP[j1],3]], "were significantly different, p=", productP," \n")
    tmp <- tmp |> tab_style(style=list(
      cell_fill(color="grey")),
      locations = cells_body(
        columns=Rframe[sortedP[j1],3], 
        rows=Rframe[sortedP[j1],2]))
  }
  
  gtsave(tmp, "CompletePops_OD_pvalue.docx", )#|>fmt_number(decimals = 3)
  
  productP=1
  cat("Significance for non-OD data \n")
  for (j1 in 1:length(sortednoODP)){
    productP = productP * RframenoOD[sortednoODP[j1],1]
    if (productP < 0.95){break}
    cat(j1,": ", namedPops[RframenoOD[sortednoODP[j1],2]], " and ", 
        namedPops[RframenoOD[sortednoODP[j1],3]], "were significantly different, p=", productP," \n")
    tmp2 <- tmp2 |> tab_style(style=list(
      cell_fill(color="grey")),
      locations = cells_body(
        columns=RframenoOD[sortednoODP[j1],3], 
        rows=RframenoOD[sortednoODP[j1],2]))
  }
  
  gtsave(tmp2, "CompletePops_noOD_pvalue.docx", )#|>fmt_number(decimals = 3)
  
}

############################################################
posterior = postOD
j1 = round(runif(1, 1, nPops))
if (useLogit){
  conv <- LC_logit(cbind(tested-surv,surv) ~ log10(Dose), weights = tested, 
                   data = ABCorrected[ABCorrected$Population==j1,], p=50)
} else {
  conv <- LC_probit(cbind(tested-surv,surv) ~ log10(Dose), weights = tested, 
                    data = ABCorrected[ABCorrected$Population==j1,], p=50)    
}

DataOrg <- Data
DataFiltered <- filter(Data, Population==j1)

controlSurv = cSurv[j1]

DataFiltered$probPred <- qnorm(1-(pnorm((log10(chainmode(posterior$lc50[,j1]))-log10(DataFiltered$Dose))*
                                  chainmode(posterior$B[,j1])))*chainmode(posterior$nSurv[,j1]))
DataFiltered$probEco = qnorm(1-(pnorm((log10(conv$dose)-log10(DataFiltered$Dose))*conv$slope))*controlSurv)

minX = min(filter(DataFiltered, Dose > 0)$Dose)
maxX = max(filter(DataFiltered, Dose > 0)$Dose)

maxY = max(c(DataFiltered$probEco,DataFiltered$probPred))
minY = min(filter(DataFiltered, Dose > 0)$probEco)
minY = min(minY, min(filter(DataFiltered, Dose > 0)$probPred))
secondYLabel = maxY - (maxY-minY)/6
if (useLogit) {
  subtit = "Using a logit model"
} else {
  subtit ="Using a probit model"
}
p1 <- ggplot()+
  geom_jitter(data = DataFiltered, aes(x = (Dose), y = qnorm(1-surv/tested)),color = "red", width=.02)+
  scale_x_continuous(trans='log10')+ 
  geom_smooth(method="lm", data = DataFiltered, aes(x = (Dose), y = probPred), color = "blue")+
  geom_smooth(method="lm", data = DataFiltered, aes(x = (Dose), y = probEco), color = "green")+
  ylab("Probits")+
  ggtitle(paste0("Bayes versus Ecotox: ", namedPops[j1]), subtitle = subtit)+
  annotate("text", x=minX*2.2, y=maxY, label= "Bayesian",
           col="blue", size=7, hjust=0)+
  annotate("text", x=minX*2.2, y=secondYLabel, label= "Ecotox",
         col="green", size=7, hjust=0)+
  annotate("segment", x = minX, xend = minX*2, y = maxY, yend = maxY,
             colour = "blue", linewidth=1.6)+
  annotate("segment", x = minX, xend = minX*2, y = secondYLabel, yend = secondYLabel,
         colour = "green", linewidth=1.6)
plot(p1)
 
