
# Chainsaw's Stat. eXtension Miscellaneous
# Data from:
# Leung SSF, et al. Annal Human Biol 1998; 25(2):169-174

# written by CH Chan aka Chainsaw Riot
# released under GPL 3

# Version history
#  1.00 - raw version
#  1.01 - Gender required to be factor
#         i.e. gender <- factor(gender,levels=c(0,1),labels=c("M","F"))
#  1.02 - add another function pair.t() for better pair.t test
#  1.03 - Function aor() and mega.confint() for easy display of glm results
#  1.04 - Function sens.cal() for quick and easy sensitivity, specificity, +LR and -LR calculation
#  1.05 - expired bmisds(), replaced with sds.cal()
#  1.06 - Added Wuhl et al's Ambulatory BP norm. Some clean up.
#  1.07 - Added layout suite. Enhanced meansd() to display also the range. A wrapper meansdr() was added.
#  1.08 - Added csx.mcnemar.test. Give an NA when 2x2 table cannot be formed.
#  1.09 - Added summ.stat()
#  1.10 - Added mm(), matlab style matrix creation.
#  1.11 - Fixed a incorrect entry in leung1998


pair.t <- function(var1,var2,rd=2) {
               return(data.frame(
                          mean=
                          c(
                            round(mean(var1),rd),
                            round(mean(var2),rd)),
                          sd=
                          c(
                            round(sd(var1),rd),
                            round(sd(var2),rd)),
                          p=c(t.test(var1,var2,paired=TRUE)$p.value,wilcox.test(var1,var2,paired=TRUE)$p.value)))
}

aor <- function (number, rd=2) {
        return(round(exp(number),rd))
}

mega.confint <- function(model,level=0.95,rd=2, aor=FALSE) {
    coef.model <- coef(model)
    confint.model <- confint(model,level=level)
    varnum <- dim(confint.model)[1]
    varnum1 <- varnum+1
    totalnum <- varnum*2
    if (aor==FALSE) {
    data.frame(
               name = names(coef.model),coef = round(as.numeric(coef.model),rd),upper = round(as.numeric(confint.model)[1:varnum],rd),lower = round(as.numeric(confint.model)[varnum1:totalnum],rd)
               ) }
    else if (aor==TRUE) {
    data.frame(
               name = names(coef.model),coef = sapply(as.numeric(coef.model),aor,rd=rd),upper = sapply(as.numeric(confint.model)[1:varnum], aor ,rd=rd),lower = sapply(as.numeric(confint.model)[varnum1:totalnum],aor,rd=rd)
               )
    }

}

bmi.cal <- function(wt, ht, htcm = TRUE) {
  if (htcm==TRUE) {
    ht.m <- ht/100
  }
  else ht.m <- ht
  return(wt/(ht.m)**2)
}

sens.cal <- function(diag,ques) {
    sens <- table(diag,ques)[4]/ (table(diag,ques)[4]+ table(diag,ques)[2])
    spec <-  table(diag,ques)[1]/ (table(diag,ques)[1]+ table(diag,ques)[3])
    ppv <-  table(diag,ques)[4]/( table(diag,ques)[4]+ table(diag,ques)[3])
    npv <-  table(diag,ques)[1]/( table(diag,ques)[1]+ table(diag,ques)[2])
    plr <- sens / (1-spec)
    nlr <- (1-sens)/spec
    return(c(sens,spec,ppv,npv,plr,nlr))
}

sds.cal <- function (value, age, male, lmsdata=leung1998) {
  output <- c()
  for(i in 1:length(value)) {
    query <- age[i] >= lmsdata$age.min & age[i] < lmsdata$age.max & male[i] ==lmsdata$male
    lms.list  <- c(lmsdata$l[query], lmsdata$m[query], lmsdata$s[query])
    output[i]<- ((((value[i]/lms.list[2])**lms.list[1] -1)/(lms.list[1]*lms.list[3])))
  }
  return(output)
}

#rev.sds.cal - to get the y value from z using binary search. For plotting the percentile chart
rsds.cal <- function (z, age, male, lmsdata=leung1998, epsilon=0.01, max.iters=50) {
  low.bound <- 0
  hi.bound <- 100
  theta <- (low.bound+hi.bound) / 2
  iters <- 0
  while (abs(sds.cal(theta, age, male, lmsdata) - z) > epsilon & iters <= max.iters) {
    if (sds.cal(theta, age, male, lmsdata) < z) {
      low.bound <- theta
    } else {
      hi.bound <- theta
    }
    theta <- (low.bound+hi.bound) / 2
    iters <- iters + 1
  }
  if (iters == max.iters) {
    return(NA)
  } else {
    return(theta)
  }
}
vec.rsds.cal <- function (z, agev, male, lmsdata=leung1998, epsilon=0.01, max.iters=50) {
  output <- c()
  for (i in 1:length(agev)) {
    output[i] <- rsds.cal(z, agev[i], male, lmsdata, epsilon, max.iters=50)
  }
  return(output)
}


deter.pct <- function(value, age, male, pctdata=sungWC) {
  output <- c()
  for(i in 1:length(value)) {
    query <- age[i] >= pctdata$criteria.min & age[i] < pctdata$criteria.max & male[i] ==pctdata$male
    cmp.pct  <- c(pctdata$p90[query])
    output[i]<- value[i] >= cmp.pct
  }
  return(output)
}

mm <- function(mj) {
  mj.sp <- unlist(strsplit(mj, ";"))
  for (i in 1:length(mj.sp)) {
    mj.row <- unlist(strsplit(mj.sp[i], " "))
    mj.row <- as.numeric(mj.row[mj.row!=""])
    mj.rowvec <- matrix(mj.row, byrow=TRUE, nrow=1)
    if (i == 1){
      output.mat <- mj.rowvec
    }
    else {
    output.mat <- rbind(output.mat, mj.rowvec)
  }
  }
  return(output.mat)
}


#layout suite

meansd <- function (x, sig.fig=3, na.rm=TRUE, rangeDis=FALSE) {
  values <- round(c(mean(x, na.rm=na.rm), sd(x, na.rm=na.rm), range(x, na.rm=na.rm)), sig.fig)
  outputText <- paste(values[1]," (",values[2],")", sep="")
  if (rangeDis == TRUE){
    outputText <- paste(outputText, " Range: ", values[3], " to ", values[4], sep="")
    }
  return(outputText)
}
meansdr <- function(x, sig.fig=3, na.rm=TRUE, rangeDis=TRUE){
  return(meansd(x, sig.fig=3, na.rm=TRUE, rangeDis=TRUE))
}

n.percent <- function(variable, na.count=TRUE, rd=1) {
  numTrue <- sum(variable, na.rm=TRUE)
  if (na.count==TRUE) {
    base <- length(variable)
  }
  else {
    base <- length(variable) - length(variable[is.na(variable)])
  }
  percentage <- round((numTrue/base)*100,rd)
  return(paste(numTrue, " (", percentage, "%)", sep=""))
}


summ.stat <- function(var1, sig.fig.meansd=3, na.rm.meansd=TRUE, rangeDis.meansd=FALSE, rd.n.percent=1, na.count.n.percent=TRUE){
  if (class(var1)=="numeric"){
    return(meansd(var1, sig.fig=sig.fig.meansd, na.rm=na.rm.meansd, rangeDis=rangeDis.meansd))
  }
  else if (class(var1)=="logical"){
    return(n.percent(var1, na.count=na.count.n.percent, rd=rd.n.percent))
  }
  else {
    return(NA)
  }
}

format.CI <- function(vector, rd=1) {
  v <- round(vector, rd)
  return(paste(v[1], " (", v[2], " to ", v[3], ")", sep=""))
}

beauti.OR <- function (x,y, sig.fig=2) {
  fish.prod <- fisher.test(x,y)
  OR <- c(fish.prod$estimate, fish.prod$conf.int[1],fish.prod$conf.int[2])
  return(format.CI(OR, rd=sig.fig))
}


beauti.meandiff <- function(x,y,sig.fig=2) {
  # y is a binary factor
  meanList <- by(x, y, FUN=mean, na.rm=TRUE)
  mean.diff <- meanList[1] - meanList[2]
  t.test.xy <- t.test(x~y)
  mean.diff.l <- t.test.xy$conf.int[1]
  mean.diff.u <- t.test.xy$conf.int[2]
  mean.diff.ci <- c(mean.diff, mean.diff.l, mean.diff.u)
  return(format.CI(mean.diff.ci, rd=sig.fig))
}


csx.mcnemar.test <- function(logicalV1, logicalV2) {
  dim.2by2 <- dim(table(logicalV1, logicalV2))
  if (dim.2by2[1] == 2 & dim.2by2[2] ==2) {
    return(mcnemar.test(logicalV1, logicalV2)[3])
  }
  else {
    return(NA)
  }
}

# validation
## generate test cases
## tc.bmi <- runif(10, 1.0, 40.0)
## tc.age <- runif(10, 0.0, 20.0)
## tc.male <- sample(c(TRUE,FALSE),10, replace=TRUE)
## sds.cal(tc.bmi, tc.age, tc.male, lmsdata)

## tc.bmi <- c(28.7946671976242, 6.35195504431613, 38.0401123196352, 35.2717898080591, 
## 29.7220048084855, 19.5666483980604, 11.5795691602398, 35.8160910827573, 
## 27.7417488575447, 39.4510587523691)

## tc.age <- c(14.1764045413584, 7.74132217280567, 19.3441010313109, 7.5099214585498, 
## 13.9197203237563, 7.28455384261906, 15.64328443259, 11.7556247673929, 
## 16.1663270508870, 5.41895815636963)

## tc. male <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE
## )

## test.results <- sds.cal(tc.bmi, tc.age, tc.male, lmsdata=lmsdata)
## epsilon <- 0.01
## gold.standard <- c(2.03, -24.64, 2.45, 3.20, 2.42, 1.74, -6.06, 2.64, 1.96, 4.33)
## diff.results <- test.results - gold.standard
## acceptable <- abs(diff.results) > epsilon
## if (sum(acceptable) == 0) {cat("all passed!\n")} else {cat("something wrong!")}

#### LMS data . norms


leung1998 <- structure(list(age.min = c(0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 
4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 
11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 
18, 0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 
6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 
13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18), age.max = c(0.25, 
0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 
8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 
15, 15.5, 16, 16.5, 17, 17.5, 18, 999, 0.25, 0.5, 1, 1.5, 2, 
2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 
10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 
16.5, 17, 17.5, 18, 999), male = c(TRUE, TRUE, TRUE, TRUE, TRUE, 
TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE), l = c(-1.868, -1.663, -1.66, -1.79, -1.931, -2.045, 
-2.133, -2.2, -2.251, -2.287, -2.313, -2.328, -2.336, -2.337, 
-2.332, -2.323, -2.311, -2.296, -2.279, -2.26, -2.241, -2.221, 
-2.2, -2.18, -2.159, -2.139, -2.12, -2.101, -2.082, -2.064, -2.047, 
-2.03, -2.013, -1.997, -1.982, -1.967, -1.952, -1.938, 0.346, 
-0.314, -0.641, -1.013, -1.243, -1.403, -1.522, -1.612, -1.681, 
-1.735, -1.776, -1.808, -1.833, -1.851, -1.868, -1.874, -1.879, 
-1.882, -1.882, -1.88, -1.877, -1.872, -1.866, -1.86, -1.854, 
-1.848, -1.841, -1.835, -1.828, -1.822, -1.816, -1.809, -1.803, 
-1.797, -1.791, -1.786, -1.781, -1.775), m = c(13.393, 16.538, 
17.055, 16.421, 16.037, 15.861, 15.741, 15.59, 15.426, 15.276, 
15.156, 15.072, 15.029, 15.026, 15.066, 15.149, 15.271, 15.423, 
15.601, 15.798, 16.009, 16.232, 16.463, 16.702, 16.945, 17.193, 
17.444, 17.697, 17.952, 18.205, 18.456, 18.704, 18.947, 19.186, 
19.418, 19.644, 19.865, 20.079, 13.125, 15.951, 16.562, 16.133, 
15.76, 15.58, 15.485, 15.39, 15.258, 15.097, 14.949, 14.84, 14.772, 
14.747, 14.77, 14.839, 14.951, 15.101, 15.284, 15.497, 15.735, 
15.993, 16.268, 16.557, 16.854, 17.155, 17.458, 17.757, 18.052, 
18.338, 18.646, 18.884, 19.143, 19.393, 19.635, 19.869, 20.096, 
20.317), s = c(0.10203, 0.07855, 0.0735, 0.07145, 0.07127, 0.07116, 
0.07153, 0.0729, 0.07546, 0.07908, 0.08355, 0.08867, 0.09419, 
0.09985, 0.10543, 0.11077, 0.11576, 0.12036, 0.12452, 0.12825, 
0.13155, 0.13443, 0.13691, 0.13902, 0.14081, 0.1423, 0.14355, 
0.14458, 0.14544, 0.14617, 0.14679, 0.14733, 0.14782, 0.14826, 
0.14868, 0.14907, 0.14945, 0.14982, 0.08773, 0.08044, 0.07914, 
0.07656, 0.07467, 0.07342, 0.07306, 0.0741, 0.07643, 0.07983, 
0.08409, 0.08897, 0.09419, 0.09948, 0.10465, 0.10956, 0.11409, 
0.11817, 0.12176, 0.12486, 0.12748, 0.12964, 0.13138, 0.13274, 
0.13376, 0.1346, 0.13499, 0.13528, 0.13542, 0.13544, 0.13538, 
0.13527, 0.13513, 0.13497, 0.13479, 0.13462, 0.13444, 0.13427
)), .Names = c("age.min", "age.max", "male", "l", "m", "s"), class = "data.frame", row.names = c(NA, 
-76L))


#< -

#dput(read.csv("~/Stat/csx/wuhl_sdbp.csv"))

wuhl.wsbp <- structure(list(age.min = c(0L, 125L, 130L, 135L, 140L, 145L, 
150L, 155L, 160L, 165L, 170L, 175L, 180L, 185L, 0L, 125L, 130L, 
135L, 140L, 145L, 150L, 155L, 160L, 165L, 170L, 175L), age.max = c(125L, 
130L, 135L, 140L, 145L, 150L, 155L, 160L, 165L, 170L, 175L, 180L, 
185L, 999L, 125L, 130L, 135L, 140L, 145L, 150L, 155L, 160L, 165L, 
170L, 175L, 999L), male = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE), l = c(-1.291, -1.1007, -0.71, -0.38, -0.075, 0.117, 0.125, 
-0.031, -0.251, -0.431, -0.463, -0.373, -0.244, -0.098, 2.107, 
1.947, 1.804, 1.686, 1.583, 1.48, 1.367, 1.259, 1.22, 1.261, 
1.332, 1.41), m = c(110.8, 111.1, 111.5, 112, 112.7, 113.7, 115.1, 
116.8, 118.6, 120.6, 122.6, 124.4, 126.2, 128, 110, 110.5, 111, 
111.6, 112.2, 113.1, 114.3, 115.6, 117, 118.3, 119.8, 121.2), 
    s = c(0.069, 0.069, 0.069, 0.068, 0.067, 0.067, 0.067, 0.066, 
    0.067, 0.068, 0.069, 0.069, 0.069, 0.069, 0.061, 0.062, 0.063, 
    0.064, 0.065, 0.066, 0.066, 0.066, 0.064, 0.06, 0.055, 0.05
    )), .Names = c("age.min", "age.max", "male", "l", "m", "s"
), class = "data.frame", row.names = c(NA, -26L))

wuhl.ssbp <- structure(list(age.min = c(0L, 125L, 130L, 135L, 140L, 145L, 
150L, 155L, 160L, 165L, 170L, 175L, 180L, 185L, 0L, 125L, 130L, 
135L, 140L, 145L, 150L, 155L, 160L, 165L, 170L, 175L), age.max = c(125L, 
130L, 135L, 140L, 145L, 150L, 155L, 160L, 165L, 170L, 175L, 180L, 
185L, 999L, 125L, 130L, 135L, 140L, 145L, 150L, 155L, 160L, 165L, 
170L, 175L, 999L), male = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE), l = c(-0.053, -0.314, -0.57, -0.807, -0.997, -1.106, 
-1.126, -1.068, -0.947, -0.795, -0.626, -0.451, -0.277, -0.1, 
1.565, 1.184, 0.823, 0.518, 0.292, 0.167, 0.186, 0.378, 0.745, 
1.272, 1.903, 2.579), m = c(93.6, 94.6, 95.6, 96.7, 97.9, 99, 
100.1, 101.3, 102.6, 104.1, 105.6, 107.2, 108.7, 110.2, 95, 95.7, 
96.4, 63.69, 645, 98.1, 98.9, 100, 101.1, 102.2, 103.4, 104.6
), s = c(0.077, 0.079, 0.08, 0.081, 0.082, 0.082, 0.081, 0.081, 
0.08, 0.079, 0.079, 0.078, 0.078, 0.078, 0.07, 0.071, 0.073, 
0.075, 0.077, 0.079, 0.08, 0.079, 0.077, 0.073, 0.07, 0.068)), .Names = c("age.min", 
"age.max", "male", "l", "m", "s"), class = "data.frame", row.names = c(NA, 
-26L))


wuhl.wdbp <- structure(list(age.min = c(0L, 125L, 130L, 135L, 140L, 145L, 
150L, 155L, 160L, 165L, 170L, 175L, 180L, 185L, 0L, 125L, 130L, 
135L, 140L, 145L, 150L, 155L, 160L, 165L, 170L, 175L), age.max = c(125L, 
130L, 135L, 140L, 145L, 150L, 155L, 160L, 165L, 170L, 175L, 180L, 
185L, 999L, 125L, 130L, 135L, 140L, 145L, 150L, 155L, 160L, 165L, 
170L, 175L, 999L), male = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE), l = c(1.345, 1.436, 1.531, 1.629, 1.711, 1.763, 1.777, 
1.74, 1.65, 1.509, 1.329, 1.136, 0.939, 0.741, 1.952, 1.915, 
1.881, 1.852, 1.832, 1.828, 1.836, 1.846, 1.835, 1.799, 1.739, 
1.671), m = c(72.3, 72.3, 72.2, 72.1, 72.1, 72.1, 72.1, 72.1, 
72.2, 72.3, 72.6, 72.8, 73.1, 73.4, 73.2, 72.8, 72.4, 72.1, 71.8, 
71.7, 71.8, 72, 72.4, 73.1, 73.9, 74.8), s = c(0.087, 0.086, 
0.086, 0.085, 0.083, 0.082, 0.081, 0.08, 0.081, 0.081, 0.082, 
0.082, 0.083, 0.083, 0.077, 0.07, 0.075, 0.08, 0.084, 0.087, 
0.089, 0.087, 0.083, 0.078, 0.071, 0.064)), .Names = c("age.min", 
"age.max", "male", "l", "m", "s"), class = "data.frame", row.names = c(NA, 
-26L))

wuhl.sdbp <- structure(list(age.min = c(0L, 125L, 130L, 135L, 140L, 145L, 
150L, 155L, 160L, 165L, 170L, 175L, 180L, 185L, 0L, 125L, 130L, 
135L, 140L, 145L, 150L, 155L, 160L, 165L, 170L, 175L), age.max = c(125L, 
130L, 135L, 140L, 145L, 150L, 155L, 160L, 165L, 170L, 175L, 180L, 
185L, 999L, 125L, 130L, 135L, 140L, 145L, 150L, 155L, 160L, 165L, 
170L, 175L, 999L), male = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE), l = c(0.44, 0.43, 0.421, 0.41, 0.398, 0.391, 0.395, 0.413, 
0.442, 0.487, 0.556, 0.647, 0.755, 0.871, 1.491, 1.276, 1.075, 
0.891, 0.705, 0.497, 0.279, 0.074, -0.091, -0.2, -0.261, -0.296
), m = c(54.3, 54.8, 55.1, 55.5, 55.8, 56, 56.2, 56.2, 56.3, 
56.5, 56.7, 56.9, 57.1, 57.3, 55.4, 55.3, 55.1, 54.8, 54.6, 54.4, 
54.3, 54.3, 54.6, 54.9, 55.1, 55.4), s = c(0.089, 0.092, 0.095, 
0.096, 0.1, 0.101, 0.101, 0.101, 0.1, 0.099, 0.097, 0.096, 0.094, 
0.093, 0.112, 0.115, 0.118, 0.12, 0.121, 0.121, 0.118, 0.114, 
0.11, 0.106, 0.147, 0.099)), .Names = c("age.min", "age.max", 
"male", "l", "m", "s"), class = "data.frame", row.names = c(NA, 
-26L))
