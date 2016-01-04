############# PART I Data Processing ###########################################

txt = readLines("/Users/chingman/Downloads/offline.final.trace.txt")

###### Function 1 ######
processLine = function(x){
  SPLIT = strsplit(x, "[;=,]")[[1]]
  if (length(SPLIT) > 10){
    MAT1 = matrix(SPLIT[-(1:10)], ncol = 4, byrow = TRUE)
    MAT2 = cbind(matrix(rep(SPLIT[c(2,4,6:8, 10)], nrow(MAT1)), ncol = 6,byrow = TRUE), MAT1)
    return(MAT2)
  }
  return(NULL)
}
tmp = lapply(txt[substr(txt, 1, 1) !="#"], processLine)
offline = as.data.frame(do.call("rbind", tmp), stringsAsFactors = F)
names(offline) = c("time", "scanMac", "posX", "posY", "posZ",
                   "orientation", "mac", "signal", "channel", "type")


###### Function 2 ######
cleanData = function(data, keepMacs = c("00:0f:a3:39:e1:c0", "00:0f:a3:39:dd:cd", "00:14:bf:b1:97:8a", "00:14:bf:3b:c7:c6", "00:14:bf:b1:97:90", "00:14:bf:b1:97:8d", "00:14:bf:b1:97:81")) {
  data[c("time", "posX", "posY", "posZ", "orientation", "signal")] = lapply(data[c("time", "posX", "posY", "posZ", "orientation", "signal")], as.numeric)
  data = data[,-c(2,5)]
  ANG = seq(0, by = 45, length = 9)
  data$RAngle = c(ANG[1:8], 0)[sapply(data$orientation, function(x) which.min(abs(x - ANG)))]
  data = data[data$type =="3",]
  data = data[,-c(8)]
  data$OTime = data$time
  data$time = data$time/1000
  class(data$time) = c("POSIXt", "POSIXct")
  data = data[data$mac %in% keepMacs, ]
  MacChannel = with(data, table(mac, channel))
  apply(MacChannel, 1, function(x) sum(x > 0))
  data = data[, "channel" != names(data)]
  return(data)
} 


offline2 = cleanData(offline, c("00:0f:a3:39:e1:c0", "00:14:bf:b1:97:8a", "00:14:bf:3b:c7:c6", "00:14:bf:b1:97:90", "00:14:bf:b1:97:8d", "00:14:bf:b1:97:81"))




############# PART II Visualization ###########################################
offline2$posXY = paste(offline2$posX, offline2$posY, sep = "-")  
InfoForSamePRM = with(offline2,
                      by(offline2, list(posXY, RAngle, mac),
                         function(x) x))


signalSummary =
  lapply(InfoForSamePRM,
         function(x) {
           info = x[1, ]
           info$medSignal = median(x$signal)
           info$avgSignal = mean(x$signal)
           info
         })


offlineSummary = do.call("rbind", signalSummary)

###### boxplot compare 2 Mac###############
library(lattice)
bwplot(signal ~ factor(angle) | mac, data = offline,
        subset = posX == 2 & posY == 12
        & mac == c("00:0f:a3:39:dd:cd", "00:0f:a3:39:e1:c0"),
        layout = c(2,1))
        
###### signal density plot function########
densityplot( ~ signal | mac + factor(angle), data = offline,
        subset = posX == 2 & posY == 2 &
               mac != "00:0f:a3:39:dd:cd" & mac !=  "00:14:bf:b1:97:81" & 
               mac != "00:14:bf:3b:c7:c6" & mac !=  "00:14:bf:b1:97:90",
        bw = 0.5, plot.points = FALSE)

###### signal distance plot function########
xyplot(signal ~ dist | factor(mac) + factor(angle),
       data = offlineSummary, pch = 19, cex = 0.3,
       xlab ="distance")[c(1,2,3), ]

###### heat plot function########
heatPlot = function(data , ma, an){
  
  data = data[data$mac == ma, ]
  data = data[data$angle == an, ]
  
  smoothSS = Tps(data[, c("posX","posY")], data$meanSignal)
  
  vizSmooth = predictSurface(smoothSS)
  
  plot.surface(vizSmooth, type = "C")
  
  points(data$posX, data$posY, pch=19, cex = 0.5)
}

parCur = par(mfrow = c(2,2), mar = rep(1, 4))
par(parCur)
a = heatPlot(offlineSet, subMacs[2], 0)
b = heatPlot(offlineSet, subMacs[2],90)
c = heatPlot(offlineSet, subMacs[1], 0)
d = heatPlot(offlineSet, subMacs[1],90)






############# PART III Best K and Cross-validatoin ###########################################

########### read online data ############
online = readLines("/Users/chingman/Downloads/online.final.trace.txt")
tmp = lapply(online[substr(online, 1, 1) !="#"], processLine)
online1 = as.data.frame(do.call("rbind", tmp), stringsAsFactors = F)
names(online1) = c("time", "scanMac", "posX", "posY", "posZ",
                   "orientation", "mac", "signal", "channel", "type")
online2 = cleanData(online1, keepMacs = unique(offlineSummary$mac))
online2$posXY = paste(online2$posX, online2$posY, sep = "-")


length(unique(online2$posXY))

tonlineXYA = table(online2$posXY, online2$RAngle)


########### create online data summary set ############
ONLINEInfoForSameRM = with(online2,
                           by(online2, list(posXY),
                              function(x) {
                                info = x[1, c("posXY", "posX","posY", "orientation", "RAngle")]
                                sigs = tapply(x$signal, x$mac, mean)
                                y = matrix(sigs, nrow = 1, ncol = 6,
                                           dimnames = list(info$posXY, names(sigs)))
                                cbind(info, y)
                              }))

onlineSummary = do.call("rbind", ONLINEInfoForSameRM)




########### create training set #######################
findTrain5 = function(fixedAng, signals){
  
  nearestAngle = function(fixedAng) {
    refs = seq(0, by = 45, length = 8)
    q = sapply(fixedAng, function(o) which.min(abs(o - refs)))
    c(refs[1:8], 0)[q]
  }
  NearestAngle =  nearestAngle(fixedAng)
  NearestAngle[NearestAngle < 0] = NearestAngle[ NearestAngle < 0 ] + 360
  NearestAngle[NearestAngle > 360] = NearestAngle[ NearestAngle > 360 ] - 360
  
  offlineSubset = signals[ signals$RAngle %in% NearestAngle, ]
  
  reshapeSS = function(data, varSignal = "signal",
                       keepVars = c("posXY", "posX","posY")) {
    SixS =
      with(data, by(data, list(posXY),
                    function(x) {
                      info = x[1, keepVars]
                      sigs = tapply(x[ , varSignal ], x$mac, mean)
                      y = matrix(sigs, nrow = 1, ncol = 6,
                                 dimnames = list(info$posXY,
                                                 names(sigs)))
                      cbind(info, y)}))
    newDataSS = do.call("rbind", SixS)
    return(newDataSS)
  }
  train6S = reshapeSS(offlineSubset, varSignal = "avgSignal")    
  return(train6S)
}



################# predict location function ####################
predXY = function(newSig, newAngles, trainData, k = 3){
  findClose = function(newSig, trainSubset) {
    diffs = apply(trainSubset[ , 4:9], 1,
                  function(x) x - newSig)
    dists = apply(diffs, 2, function(x) sqrt(sum(x^2)) )
    closest = order(dists) 
    CLOSE = trainSubset[closest, 1:3 ]
    return(CLOSE)
  }
  
  closeXY = list(length = nrow(newSig))
  for (i in 1:nrow(newSig)) {
    train6S = findTrain5(newAngles[i], trainData)
    closeXY[[i]] =
      findClose(newSig = as.numeric(newSig[i, ]), train6S)
  }
  estPO = lapply(closeXY, function(x) sapply(x[ , 2:3],
                                             function(x) mean(x[1:k])))
  estPO = do.call("rbind", estPO)
  return(estPO)
}


##############calculate error rate function #####################
calError =
  function(estPO, actualPO)
    sum((estPO - actualPO)^2)

############## data frame with 6 signals #######################
Reshape = function(data, varSignal = "signal", sampAngle = FALSE, 
                   keepVars = c("posXY", "posX","posY", "orientation", "RAngle")) {
  SixSig = with(data, 
                by(data, list(posXY),
                   function(x) {
                     
                     if(sampAngle == TRUE){
                       x = x[x$RAngle == sample(seq(0, by = 45, length = 8), size = 1), ]
                     }
                     info = x[1, keepVars]
                     SS = tapply(x[ , varSignal ], x$mac, mean)
                     y = matrix(SS, nrow = 1, ncol = 6,
                                dimnames = list(info$posXY,
                                                names(SS)))
                     cbind(info, y)
                   }))
  newDataSet = do.call("rbind", SixSig)
  return(newDataSet)
}

############ 2 fold online and offline #########################
v = 5
offline6SigSummary = Reshape(offline2, sampAngle = TRUE, keepVars = c("posXY", "posX","posY", "orientation", "RAngle")) 
offlineSummary
RanOff = sample(unique(offlineSummary$posXY))
RanOff = matrix(RanOff, ncol = v,
                nrow = floor(length(RanOff)/v))

########### calculate 20 k's error ###########################
K = 20
err = rep(0, K)
for (j in 1:v) {
  onlineFold = subset(offline6SigSummary,
                      posXY %in% RanOff[ , j])
  offlineFold = subset(offlineSummary,
                       posXY %in% RanOff[ , -j])
  actualFold = onlineFold[ , c("posX", "posY")]
  for (k in 1:K) {
    estFold = predXY(newSig = onlineFold[ , 6:11],
                     newAngles = onlineFold[ , 4],
                     offlineFold , k = k)
    err[k] = err[k] + calError(estFold, actualFold)
  }
}


estOn7 = predXY(newSig = onlineSummary[ , 6:11],
                newAngles = onlineSummary[ , 4],
                offlineSummary, k = 6)

actualOn = onlineSummary[ , c("posX", "posY")]
ErrorOn = sapply(list(estOn6), calError, actualOn)




