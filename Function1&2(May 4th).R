## Function 1
processLine2 = function(x){
  SPLIT = strsplit(x, "[;=,]")[[1]]
  
  if (length(SPLIT) > 10){
    MAT1 = matrix(SPLIT[-(1:10)], ncol = 4, byrow = TRUE)
    MAT2 = cbind(matrix(rep(SPLIT[c(2,4,6:8, 10)], nrow(MAT1)), ncol = 6,byrow = TRUE), MAT1)
    return(MAT2)
  }
  return(NULL)
}

## Function 2

cleanData = function(data, keepMacs = c("00:0f:a3:39:e1:c0", 
                                        "00:0f:a3:39:dd:cd", 
                                        "00:14:bf:b1:97:8a", 
                                        "00:14:bf:3b:c7:c6", 
                                        "00:14:bf:b1:97:90", 
                                        "00:14:bf:b1:97:8d", 
                                        "00:14:bf:b1:97:81")) {
  
  numericVariable = c("time","posX", "posY", "posZ",
                      "orientation",  "signal")
  data[numericVariable] = lapply(data[numericVariable], as.numeric)
  data = data[data$mac %in% keepMacs, ]
  data = data[data$type == "3", ] 
  data = data[, "type" != names(data)]  
  data = data[, !(names(data) %in% c("scanMac", "posZ"))]
  OldTime = data$time
  NewTime = data$time/1000
  class(NewTime) = c("POSIXt", "POSIXct")
  rf = seq(0, by = 45, length = 9)
  q = sapply(data$orientation, function(o) which.min(abs(o - rf)))
  angles= c(rf[1:8], 0)[q]
  
  DF = data.frame(NewTime, data$posX, data$posY, data$orientation, data$mac, data$signal, OldTime, angles, stringsAsFactors = FALSE)
  names(DF) = c("time", "posX", "posY", "orientation", "mac", "signal", "oringalTime", "angle")
  
  return(DF)
}