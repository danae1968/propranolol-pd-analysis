pupInterp <- function(data,samplingRateHz,interp_method,msToSmooth)
{
  
library('readr')
library('dplyr')
library('tidyr')
library('zoo')
library('car')
library('ggplot2')
library('splines')
library('data.table')
library('reshape2')
library('e1071')
library('pracma')
library('forecast')

#samplingRateHz <- 500
samplingInterval <- 1000/samplingRateHz 

     # data = gaze_data[gaze_data$subject == subNo & gaze_data$trial == trialNo, ]# & gaze_data$phase == phaseNo, ]
    
    if (nrow(data)>0) {
      
# if (pup_type == 'resi'){
#   data$pupil <- data$pup_resi
# }
      
      findBlinks <- data$pup_resi==0
      #findBlinks <- is.na(data$pup_resi)
      pupilSize <- data$pup_resi
      
      findBlinks[findBlinks==TRUE] <- 1
      findBlinks[findBlinks==FALSE] <- 0
      findBlinkStarts <- diff(findBlinks)
      
      blinkStarts <- -1*(which(findBlinkStarts %in% 1)) #negative index of blink start (start is negative, end of blink positive)
      blinkEnds <- (which(findBlinkStarts %in% -1))
      blinkCandidates <- c(blinkStarts, blinkEnds)
      
      blinkCandidates <- blinkCandidates[order(abs(blinkCandidates))]
      
      if (isempty(blinkCandidates)){pup_interp <- data$pup_resi}
      else {
      # edge conditions - if recording starts or ends with a blink:
      # starts with a blink?
      if(blinkCandidates[1]>0 && data$pup_resi[1]==0) {
        blinkCandidates <- c(-1, blinkCandidates)
      }
      
      # ends with a blink?
      endCan <- length(blinkCandidates)
      endSamples <- dim(data)[1]
      if(blinkCandidates[endCan]<0 && pupilSize[endSamples]==0) {
        blinkCandidates <- c(blinkCandidates, endSamples)
      }
      
      
      ## smooth data for algorithm (not final smoothing)
      #msToSmooth <- 1500
      samplesToSmooth <- ceiling(msToSmooth/samplingRateHz)
      
      gapInterval <- 250 #min number of ms between consecutive blinks (otherwise conctenate blinks) was 100 switched 01/12
      
      # smoothing with moving window
      smooth_data <- moving_average_pupil(data$pup_resi, n = samplesToSmooth, centered = TRUE)
      ## smoothing as hershman
      #smooth_data    <- ma(data$pup_resi, samplesToSmooth)   
      #smooth_data[1] <- data$pup_resi[1]

      smooth_data[smooth_data==0] <- NA;  # replace zeros with NaN values
      
      dpupil <- diff(smooth_data)
      
      
      ## finding blinks onset and offset
      blink <- 1
      blinkPositions <- data.frame(blinkCandidates)
      prevOffset <- -1
      
      while (blink < nrow(blinkPositions)) {
        # position of onset candidate
        onsetCandidate <- blinkPositions$blinkCandidates[blink]
        blink <- blink + 1
        
        # position of offset candidate
        offsetCandidate <- blinkPositions$blinkCandidates[blink]
        blink=blink + 1
        
        # find blink onset
        dataBeforeBlink <- dpupil[2:abs(onsetCandidate)]
        blinkOnset <- max(which(dataBeforeBlink>0)[(which(dataBeforeBlink>0))!= max(which(dataBeforeBlink>0))])
        blinkOnset <- blinkOnset
        
        # for the case that the sample begins with a blink:
        if (!(is.finite(blinkOnset))) {
          blinkOnset <- abs(onsetCandidate)
        } 
        
        # find blink offset
        dataAfterBlink <- dpupil[abs(offsetCandidate):endSamples]
        blinkOffset <- offsetCandidate + min(which(dataAfterBlink<0)[(which(dataAfterBlink<0)) != min(which(dataAfterBlink<0))])
        blinkOffset <- blinkOffset + 25
        
        # for the case that the sample ends with a blink
        if (!(is.finite(blinkOffset))) {
          blinkOffset <- endSamples
        } 
        
        # set onset to be equal to previous offset in cases where blinks are too close together:
        if (blink>3 & samplingInterval*blinkOnset - samplingInterval*prevOffset <= gapInterval) {
          blinkOnset <- prevOffset
        }
        
        prevOffset <- blinkOffset - 1
        # insert onset into results array
        blinkPositions$onOff[blink-2] <- blinkOnset
        # insert offset into array
        blinkPositions$onOff[blink-1] <- blinkOffset - 1
        
        
      }
      
      # find cases in which there are an odd number of repetitions in order to keep the first one
      findOddCases <- as.data.frame(table(blinkPositions$onOff))
      oddCases <- filter(findOddCases, Freq %in% seq(from=3, to=99, by=2))
      numToRm <- oddCases$Freq -1
      keep <- rep(1, nrow(oddCases))
      repVec <- c(rbind(keep, numToRm))
      basVec <- rep(c(0,1), times=length(numToRm))
      keepVec <- rep(basVec, repVec)
      rmVec <- which(keepVec==1)
      allOddCases <- blinkPositions$blinkCandidates[blinkPositions$onOff %in% oddCases$Var1]
      rmVals <- allOddCases[rmVec]
      blinkPositions <- filter(blinkPositions, !(blinkCandidates %in% rmVals))
      
      # removing duplications (i.e. concatenated blinks)
      blinkPositionsClean <- blinkPositions[!(duplicated(blinkPositions$onOff) | duplicated(blinkPositions$onOff, fromLast=TRUE)),]
      
      if (isempty(blinkPositionsClean)){ pup_interp <- data$pup_resi}
      
      else {
      # classifying as onset or offset
      blinkPositionsClean$type[blinkPositionsClean$blinkCandidates<0] <- 'start'
      blinkPositionsClean$type[blinkPositionsClean$blinkCandidates>0] <- 'end'
      blinkPositionsClean <- select(blinkPositionsClean, onOff, type, blinkCandidates)
      nrows <- nrow(blinkPositionsClean)/2
      blinkPositionsClean$index <- rep(seq(1, nrows), each=2)
      blinkPos<-blinkPositionsClean
      #maxInd <- nrow(pupilSize)
      maxInd=endSamples

      pup_blink <- data$pup_resi

      # values between blink onset and offset are turned to NA
      
      for (i in seq(1,length(blinkPos$index),2)){
        pup_blink[blinkPos$onOff[i]:blinkPos$onOff[i+1]] <- NA
        
      }
      #use na.approx for interpolation or na.spline for cubic interpolation
      
      if(interp_method== 'linear'){
        pup_interp <- na.approx(pup_blink,maxgap = Inf)}
        else if (interp_method== 'cubic'){
          pup_interp <- na.spline(pup_blink, maxgap = Inf)
      }
      }
      }
    } else { warning('no data for subject and trial')
    pup_interp <- NA}
      
      
return(pup_interp)

}
 
