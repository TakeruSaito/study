#��Ǝ菇
#�P�D�K�͂̔��萔�������߂�
#�������߂邽�߁A�ŏ����@��p����
#�ŏ����@�͊֐�lsfit(x,y)��p����B�@x�F�����ϐ��̍s��@�@y�G������ϐ��̃x�N�g��
#f <- lsfit(�\���H��, �K��)
#f$coefficients �ŕ\�������@Intercept  V1��y = ax + b��b,a�ɊY������B

StepOne <- function(pbl){ #csv�̌��f�[�^���̂܂�

  TriAlfa <- lsfit(pbl$Estimated, pbl$FPTrial)$coeff
  AppAlfa <- lsfit(pbl$Estimated, pbl$FPApproxi)$coeff 

  return (invisible(TriAlfa))
}

#�Q�D�e�ϓ��v���Ԃ̑��֐��𒲂ׁA���d���`����r������
#source("http://aoki2.si.gunma-u.ac.jp/R/src/tolerance.R", encoding="euc-jp")
#��R�R���\�[���Ƀy�[�X�g�����̂��Atolerance(x)�@x�F�����ϐ������̃f�[�^�s��@���g��
#tolerance(x)�̏o�͌��ʂ̂����Atolerance�̒l���傫�����̂��폜����

StepTwo <- function(AllFact){

  data <- AllFact  #����AllFact�͍H���E�K�͂��������A�S�����v��

  data[is.na(data)] <- 0
  buffer <- data

  for(i in 1:length(data)){
   if(sum(data[,i]) == 0){
     buffer <- buffer[, -i]
   }
  }

  data <- buffer

  vect <- 0
  vname <- 0
  for(i in 1:(length(data)-1)){
    for(j in (i + 1):length(data)){
      vect <- c(vect, cor(data[,i], data[,j]))
      vname <- c(vname, paste(as.character(i), "-", as.character(j)))
    }
  }

  vect <- vect[-1]
  vname <- vname[-1]
  names(vect) <- vname
  sortList <- order(-abs(vect))

  return(invisible( vect[sortList]))

}

#�v�Z�p�ɁA�\�[�g���ꂽ�f�[�^�̂�������Ȃ������폜����悤�ɏC������StepTwo
Multico <- function(AllFact){

  data <- AllFact  #����AllFact�͍H���E�K�͂��������A�S�����v��

  data[is.na(data)] <- 0
  buffer <- data
  ts<- 0
  for(i in 5:length(data)){
   if(sum(data[,i]) == 0){
    ts <- c(ts, i)
   }
  }

  ts <- ts[-1]
  if(length(ts) != 0){
    for(i in length(ts):1){
      buffer <- buffer[,-ts[i]]
    }
  }
  data <- buffer

  vect <- 0
  for(i in 5:(length(data)-1) ){
    for(j in (i + 1):length(data) ){
      x <- abs( cor(data[,i] , data[,j]) );
      if(is.na(x) || x >= 0.8){
        vect <- rbind(vect, c(i,j));
      }
    }
  }
  
  if(length(vect) != 0){
    vname <- 0;  
    for(i in 2:length(vect[1,])){
      data[,vect[i,2]] <- 0 
    }
  }
  for(i in length(data):5){
   if(sum(data[,i]) == 0){
    data <- data[,-i]
   }
  }
#  return (vect)
  return(data)
}



#�R�D�H���Ɗe�ϓ��v���̕ϐ��̊Ԃ̉�A�W�������߂�
#��A�W����lm(formula, data)��p����
#�g������Flm(�H��?�ϓ��v��, data=(�H���ƕϓ��v�����i�[�����s��))
#��A�W����lm()���������ϐ�����Coefficients�̒��AEstimate�s�ϓ��v������̒l

StepThree <- function(SelectedData){
  kaiki <- 0
  data2 <- SelectedData�@�@#���d���`����r�������S�f�[�^
  data2[is.na(data2)] <- 0
  
  for(i in 5:length(data2)){
    kaiki <- rbind(kaiki, lm((Actual-Estimated)?data2[,i], data2)$coefficients)
  }
 
  kaiki <- kaiki[-1,]
  kaiki[is.na(kaiki)] <- 0

  kname <- names(data2)[5:length(data2)]
  dimnames(kaiki) <- list(kname, names(kaiki[1,]))
  return(invisible(kaiki))
}

#3.�ŉ�A�W���ł͂Ȃ��C���֌W����Ԃ��p�^�[��
StepThreeCor <- function(SelectedData){
  soukan <- 0
  data2 <- SelectedData
  data2[is.na(data2)] <- 0

  for(i in 5:length(data2)){
    soukan <- rbind(soukan , cor(data2$Actual - data2$Estimated, data2[,i]))
  }

  soukan <- soukan[-1,]
  soukan[is.na(soukan)] <- 0

  sname <- names(data2)[5:length(data2)]
  #dimnames(soukan) <- list(sname, )
  names(soukan) <- sname
  return (invisible(soukan))

}

#3�D�ŒP��A���͂ł͂Ȃ��C�d��A���͂������ꍇ
StepThreeMulti <- function(SelectedData){
# kaiki <- 0
  library("QuantPsyc", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
  data2 <- SelectedData
  data2[is.na(data2)] <- 0
 
  kaiki <- lm.beta(lm((data2$Actual - data2$Estimated)?., data2[-1:-4]))

  kaiki[is.na(kaiki)] <- 0
  
  return(invisible(kaiki))
}


#�S�D�ŏ����@��p���āA���ς��胂�f�������쐬����
#���ݑz�肵�Ă��郂�f������(�H��)�@= ���@�~�@FP(�K��)�@�~�@(�P�{�@���@�~�@��{(�ϓ��v��)�@�~�@(�e�ϓ��v���̉�A�W��)})
#�Ȃ̂ŁA(�H��)- ���~FP�@=�@���@�~�@FP�~�@���@�~�@(�ϓ��v���Ƃ��̉�A�W���̐ρ@�̑��a)�Ǝ��ό`���A
#��lsfit(��FP(�ϓ��v���Ƃ��̉�A�W���̐ρ@�̑��a), (�H��-��FP))�Ƃ��A�P.�Ɠ��l�ɂ��āA�W�������Z�o����

#(�H��)- ���~FP = data$Actual - (TriAlfa[2] * data$FPTrial)

StepFour <- function(data, TriAlfa, kaiki){

  FactMulCoeff <- 0
  for(i in 5:length(data)){ #����data�͑��d���`����r�������H�������̃v���W�F�N�g�f�[�^
    FactMulCoeff <- c(FactMulCoeff, c(data[,i] * abs(kaiki[i - 4, 2])))
  }

  FactMulCoeff <- FactMulCoeff[-1]

  sig <- matrix(FactMulCoeff, length(data[,1]), (length(data[1,]) - 4))

  sigma <- 0
  for(i in 1:length(sig[1,])){
    sigma <- sigma + sig[,i]
  }

  leftFormula <- data$Actual - (TriAlfa[2] * data$FPTrial)
  rightFormula <- TriAlfa[2] * data$FPTrial * sigma

  x <- data.frame(LEFT=leftFormula, RIGHT=rightFormula);
  return( invisible(lm(leftFormula?rightFormula, data=x)$coeff ))

}


StepFourCor <- function(data, TriAlfa, soukan){

  FactMulCoeff <- 0
  for(i in 5:length(data)){ #����data�͑��d���`����r�������H�������̃v���W�F�N�g�f�[�^
    FactMulCoeff <- c(FactMulCoeff, c(data[,i] * soukan[i - 4]))
  }

  FactMulCoeff <- FactMulCoeff[-1]

  sig <- matrix(FactMulCoeff, length(data[,1]), (length(data[1,]) - 4))

  sigma <- 0
  for(i in 1:length(sig[1,])){
    sigma <- sigma + sig[,i]
  }

  leftFormula <- data$Actual - (TriAlfa[2] * data$FPTrial)
  rightFormula <- TriAlfa[2] * data$FPTrial * sigma

  x <- data.frame(LEFT=leftFormula, RIGHT=rightFormula);
  return( invisible(lm(leftFormula?rightFormula, data=x)$coeff ))

}


StepFourMulti <- function(data, TriAlfa, kaiki){

  FactMulCoeff <- 0
  for(i in 5:length(data)){ #����data�͑��d���`����r�������H�������̃v���W�F�N�g�f�[�^
    FactMulCoeff <- c(FactMulCoeff, c(data[,i] * kaiki[i - 4]))
  }

  FactMulCoeff <- FactMulCoeff[-1]

  sig <- matrix(FactMulCoeff, length(data[,1]), (length(data[1,]) - 4))

  sigma <- 0
  for(i in 1:length(sig[1,])){
    sigma <- sigma + sig[,i]
  }

  leftFormula <- data$Actual - (TriAlfa[2] * data$FPTrial)
  rightFormula <- TriAlfa[2] * data$FPTrial * sigma

  x <- data.frame(LEFT=leftFormula, RIGHT=rightFormula);
  return( invisible(lm(leftFormula?rightFormula, data=x)$coeff ))

}

#�T�D�N���X�o���f�[�V�����Ő��萸�x���Z�o����B
#R�ł̃N���X�o���f�[�V�����̃��\�b�h�͂���悤�����A�g�������킩���̂ŁA�ň���ł���

StepFive <- function(pbl, checkData){ # pbl:�����Ώۃv���W�F�N�g�폜�ς݂����d���`���r���ς݂̑S�f�[�^
  TriAlfa <- StepOne(pbl)
  pbl <- Choice(pbl)
  kaiki <- abs(StepThree(pbl))
  predict <- StepFour(pbl, TriAlfa, kaiki)
  FactMulCoeff <- 0

  for(i in 5:length(pbl)){ 
    FactMulCoeff <- c(FactMulCoeff, c(checkData[i] * kaiki[i - 4, 2]))
  }

  FactMulCoeff <- FactMulCoeff[-1]
  sigma <- 0
  for(i in 1:length(FactMulCoeff)){
    sigma <- sigma + as.numeric(FactMulCoeff[i])
  }

#  sig <- matrix(FactMulCoeff, length(pbl[,1]), (length(pbl[1,]) - 4))

#  sigma <- 0
#  for(i in 1:length(sig[1,])){
#    sigma <- sigma + sig[,i]
#  }

  ret <- TriAlfa[2] * checkData$FPTrial * (1 + predict[2] * sigma)
  ret <- abs(ret)
  ret <- c(ret, checkData$Actual, ret - checkData$Actual, TriAlfa[2], checkData$FPTrial, predict[2], sigma)
  return(invisible(ret)) 
}


StepFiveCor <- function(pbl, checkData){ # pbl:�����Ώۃv���W�F�N�g�폜�ς݂����d���`���r���ς݂̑S�f�[�^
  TriAlfa <- StepOne(pbl)
  pbl <- Choice(pbl)
  soukan <- StepThreeCor(pbl)
  predict <- StepFourCor(pbl, TriAlfa, soukan)
  FactMulCoeff <- 0

  for(i in 5:length(pbl)){ 
    FactMulCoeff <- c(FactMulCoeff, c(checkData[i] * soukan[i - 4]))
  }

  FactMulCoeff <- FactMulCoeff[-1]
  sigma <- 0
  for(i in 1:length(FactMulCoeff)){
    sigma <- sigma + as.numeric(FactMulCoeff[i])
  }

  ret <- TriAlfa[2] * checkData$FPTrial * (1 + predict[2] * sigma)
  ret <- abs(ret)
  ret <- c(ret, checkData$Actual, ret - checkData$Actual, TriAlfa[2], checkData$FPTrial, predict[2], sigma)
  return(invisible(ret)) 
}


StepFiveMulti <- function(pbl, checkData){ # pbl:�����Ώۃv���W�F�N�g�폜�ς݂����d���`���r���ς݂̑S�f�[�^
  TriAlfa <- StepOne(pbl)
  pbl <- Choice(pbl)
  kaiki <- StepThreeMulti(pbl)
  predict <- StepFourMulti(pbl, TriAlfa, kaiki)
  FactMulCoeff <- 0

  for(i in 5:length(pbl)){ 
    FactMulCoeff <- c(FactMulCoeff, c(checkData[i] * kaiki[i - 4]))
  }

  FactMulCoeff <- FactMulCoeff[-1]
  sigma <- 0
  for(i in 1:length(FactMulCoeff)){
    sigma <- sigma + as.numeric(FactMulCoeff[i])
  }

  ret <- TriAlfa[2] * checkData$FPTrial * (1 + predict[2] * sigma)
  ret <- abs(ret)
  ret <- c(ret, checkData$Actual, ret - checkData$Actual, TriAlfa[2], checkData$FPTrial, predict[2], sigma)
  return(invisible(ret)) 
}


Cross <- function(pblData){ #pblData�͑��d���`���r���ς݂̑S�v���W�F�N�g�f�[�^
�@�@ret <- 0; 

  for(i in 1:length(pblData[,1])){
    Study <- pblData[-(i),];
    Est <- pblData[i, ];
    ret <-  rbind(ret, StepFive(Study, Est));
  }
  ret <- ret[-1,]
  return(ret)
}

CrossCor <- function(pblData){ #pblData�͑��d���`���r���ς݂̑S�v���W�F�N�g�f�[�^
�@�@ret <- 0; 

  for(i in 1:length(pblData[,1])){
    Study <- pblData[-(i),];
    Est <- pblData[i, ];
    ret <-  rbind(ret, StepFiveCor(Study, Est));
  }
  ret <- ret[-1,]
  return(ret)
}


CrossMulti <- function(pblData){ #pblData�͑��d���`���r���ς݂̑S�v���W�F�N�g�f�[�^
�@�@ret <- 0; 

  for(i in 1:length(pblData[,1])){
    Study <- pblData[-(i),];
    Est <- pblData[i, ];
    ret <-  rbind(ret, StepFiveMulti(Study, Est));
  }
  ret <- ret[-1,]
  return(ret)
}


Choice <- function(pbl){ #���d���`���r���ς݂̑S�f�[�^
�@�@Data <- pbl; Data[is.na(Data)] <- 0 #�����l��r��
  ret <- 0 #�Ԃ�l�p�̕ϐ�

  for(i in 5:length(Data)){
    ret <- c(ret, cor((Data$Actual -  Data$Estimated), Data[, i])) #ret�ɍH���Ɨv���Ԃ̑��֌W������
  }
  
  ret <- ret[-1]�@#�������̎���0���폜
  names(ret) <- names(pbl[1,5:length(pbl)]) #���x���t

  sortList <- order(-abs(ret)) #�傫�����Ƀ\�[�g
  ret <- Data
  s <- c(1:4); sortList <- c(s, sortList+4)
  ret <- ret[sortList]
  return(invisible(ret))
}

#plot(SecondResult[SecondSortList,3], xaxt="n", xlab="" , ylab="")
#axis(1,at=1:11, labels=label[SecondSortList])
#abline(SecondMean, 0, col="blue")

#plot(SecondResult[SecondSortList,3], xaxt="n", xlab="" , ylab="")
#axis(1,at=1:11, labels=label[SecondSortList])
#abline(Secondmedian, 0, col="red")

ChoiceMKtwo <- function(pbl){ #���d���`���r���ς݂̑S�f�[�^
�@�@Data <- pbl; Data[is.na(Data)] <- 0 #�����l��r��
  ret <- 0 #�Ԃ�l�p�̕ϐ�

  for(i in 5:length(Data)){
    ret <- c(ret, cor((Data$Actual -  Data$Estimated), Data[, i])) #ret�ɍH���Ɨv���Ԃ̑��֌W������
  }
  
  ret <- ret[-1]�@#�������̎���0���폜
  names(ret) <- names(pbl[1,5:length(pbl)]) #���x���t

  sortList <- order(-abs(ret)) #�傫�����Ƀ\�[�g

  return(invisible(ret))
}


Bunsan <- function(Data){
  ret <- 0
  for(i in 6:length(Data)){
    ret <- c(ret , sqrt(variance(CrossMulti(Data[,1:i])[,3])))
  }
  #ret <- ret[-1]
  return (ret)
}

variance <- function(x) var(x)*(length(x)-1)/length(x)

MakeModel <- function(PBL){
  SortedData <- Choice(PBL) #�H���ϓ��v�����H���덷�Ƃ̑��֌W���̍������Ƀ\�[�g
  ChoiceCoeff <- ChoiceMKtwo(PBL) #�e�ϓ��v���̑��֌W������
  
  multico <- Multico(SortedData) #���d���`����r������
  
#  write.table(ChoiceCoeff, file = "output.txt", append = TRUE, quote = FALSE);
 write.csv(ChoiceCoeff, file = "SubData/MetrixInfluence.csv", quote = FALSE, col.names = FALSE);
  
  MeanRes <- 0
  MedianRes <- 0
  VarRes <- 0
  for(i in 6:length(Data)){ #���ς���H���̌덷��result�ɑ��
    result <- CrossMulti(Data[, 1:i])[, 3]
#    MeanRes <- c(MeanRes, mean(result))
    MedianRes <- c(MedianRes, median(result))
    VarRes <- c(VarRes, sqrt(variance(result)))
  }
  
  MeanRes <- MeanRes[-1]
  MedianRes <- MedianRes[-1]
  VarRes <- VarRes[-1]

  Num <- c(2:(length(VarRes)+1))
  names(VarRes) <- Num
  names(MedianRes) <- Num

  VarSort <- order(abs(VarRes))

  VarRes <- VarRes[VarSort]
  lim <- 0
  for(i in 1:length(VarRes)){
    if(VarRes[i] <= 500){ lim <- i }
  } 

  MedianSort <- VarSort[1:lim]
  
  MedianRes <- MedianRes[MedianSort]
  MedianSort <- order(abs(MedianRes))
  MedianRes <- MedianRes[MedianSort]
  return (MedianRes[1])
}

CalcManHour <- function(PBL,i){
#  error <- MakeModel(PBL[-length(PBL[,1]),])

  error <- MakeModel(PBL)
  SortedData <- Choice(PBL)
  multico <- Multico(SortedData) #���d���`����r������
#  ret <- StepFiveMulti(multico[-length(PBL[,1]),],multico[length(PBL[,1]),])
  num <- as.integer(names(error))+4

  Study <- multico[-(i),1:num];
  Est <- multico[i, 1:num];
  ret <-  StepFiveMulti(Study, Est);
  
#  names(ret) <- c("Estimate", "Actual","Error" , "Alfa", "FPTrial", "predict", "Sigma");
#  write.table(ret, file = "output.txt", append = TRUE, quote = FALSE);
#  write.table(ret, file = "SubData/HyperParameters.csv", append = TRUE, quote = FALSE)
  names(ret) <- c(NULL, NULL,NULL , NULL, NULL, NULL, NULL);
  AandB <- 0
  AandB <- rbind(AandB, ret)
  write.table(AandB[-1,], file = "SubData/HyperParameters.csv", append = TRUE, quote = FALSE)
  
options(scipen=5);return (ret[3]  - error)
#  return(ret[,3] - error)
#  return (length(multico))
#  return (ret[3]- error)
}

#MakeModel(PBLData)
#CalcManHour(PBLData,3)
#system.time(CalcManHour(PBLData))
m <- 0
for(i in 1:length(PBLData[,1])){
#  write.table(label[i], file = "output.txt", append = TRUE, quote = FALSE);
  m <- rbind(m , CalcManHour(PBLData,i))
}
m <- m[-1]

plot(m, xaxt="n", pch = 20, xlab="projects", ylab="Error of scheduled man-hours and actual man-hours");
axis(side = 1, at = 1:length(m), labels = label)
text(1:length(m), m - 30, ceiling(m))