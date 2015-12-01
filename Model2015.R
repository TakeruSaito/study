#作業手順
#１．規模の比例定数αを求める
#αを求めるため、最小二乗法を用いる
#最小二乗法は関数lsfit(x,y)を用いる。　x：説明変数の行列　　y；被説明変数のベクトル
#f <- lsfit(予測工数, 規模)
#f$coefficients で表示される　Intercept  V1がy = ax + bのb,aに該当する。

StepOne <- function(pbl){ #csvの元データそのまま

  TriAlfa <- lsfit(pbl$Estimated, pbl$FPTrial)$coeff
  AppAlfa <- lsfit(pbl$Estimated, pbl$FPApproxi)$coeff 

  return (invisible(TriAlfa))
}

#２．各変動要因間の相関性を調べ、多重線形性を排除する
#source("http://aoki2.si.gunma-u.ac.jp/R/src/tolerance.R", encoding="euc-jp")
#をRコンソールにペーストしたのち、tolerance(x)　x：説明変数だけのデータ行列　を使う
#tolerance(x)の出力結果のうち、toleranceの値が大きいものを削除する

StepTwo <- function(AllFact){

  data <- AllFact  #このAllFactは工数・規模を除いた、全調整要因

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

#計算用に、ソートされたデータのうちいらない分を削除するように修正したStepTwo
Multico <- function(AllFact){
  library("car", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
  data <- AllFact  #このAllFactは全調整要因

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
  
  data <- PrepForVIF(data)
  
  vect <- 0
  for(i in 5:(length(data) - 3) ){
#    for(j in (i + 1):(length(data) - 1) ){
    x <- 0
    loop <- floor( (length(data) - i)/3 )
    for(j in 1:(loop) ){
#      x <- abs( cor(data[,i] , data[,j]) );
      if(j >= loop){
        x <- c(x, vif( lm( data[,i] ? . , data[, (i + (3*j - 2)):length(data)] ) ) )
      }else{
        x <- c(x, vif( lm( data[,i] ? . , data[, (i + (3*j - 2)):(i + (3*j))] ) ) )
      }
    }
    x <- x[-1]
    
    for(k in 1:length(x)){
      if(length(x) != 0){
        if(is.na(x[k]) || x[k] >= 10){
          vect <- rbind(vect, c(i,(j + k -1)))
        }
      }
    }
  }
  
  
  if(length(vect) >= 2 ){
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

PrepForVIF <- function(PBLData){
  buffer <- PBLData
  DeleteIndex<- 0
  RegrCoeff <- lm(buffer[,5] ? . , buffer[,5:length(buffer)])$coeff
  for(i in 1:length(RegrCoeff)){
    if(is.na(RegrCoeff[i])){
      DeleteIndex <- c(DeleteIndex, i + 3)
    }
  }
  
  DeleteIndex <- DeleteIndex[-1]
  if(length(DeleteIndex) != 0){
    for(i in length(DeleteIndex):1){
      buffer <- buffer[,-DeleteIndex[i]]
    }
  }
  
  return(buffer)
}

#３．工数と各変動要因の変数の間の回帰係数を求める
#回帰係数はlm(formula, data)を用いる
#使い方例：lm(工数?変動要因, data=(工数と変動要因を格納した行列))
#回帰係数はlm()を代入した変数中のCoefficientsの中、Estimate行変動要因名列の値

StepThree <- function(SelectedData){
  kaiki <- 0
  data2 <- SelectedData　　#多重線形性を排除した全データ
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

#3.で回帰係数ではなく，相関係数を返すパターン
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

#3．で単回帰分析ではなく，重回帰分析をした場合
StepThreeMulti <- function(SelectedData){
# kaiki <- 0
  library("QuantPsyc", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
  data2 <- SelectedData
  data2[is.na(data2)] <- 0
 
  kaiki <- lm.beta(lm((data2$Actual - data2$Estimated)?., data2[-1:-4]))

  kaiki[is.na(kaiki)] <- 0
  
  return(invisible(kaiki))
}


#４．最小二乗法を用いて、見積もりモデル式を作成する
#現在想定しているモデル式は(工数)　= α　×　FP(規模)　×　(１＋　β　×　Σ{(変動要因)　×　(各変動要因の回帰係数)})
#なので、(工数)- α×FP　=　α　×　FP×　β　×　(変動要因とその回帰係数の積　の総和)と式変形し、
#をlsfit(αFP(変動要因とその回帰係数の積　の総和), (工数-αFP))とし、１.と同様にして、係数βを算出する

#(工数)- α×FP = data$Actual - (TriAlfa[2] * data$FPTrial)

StepFour <- function(data, TriAlfa, kaiki){

  FactMulCoeff <- 0
  for(i in 5:length(data)){ #このdataは多重線形性を排除した工数等込のプロジェクトデータ
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
  for(i in 5:length(data)){ #このdataは多重線形性を排除した工数等込のプロジェクトデータ
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
  for(i in 5:length(data)){ #このdataは多重線形性を排除した工数等込のプロジェクトデータ
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

#５．クロスバリデーションで推定精度を算出する。
#Rでのクロスバリデーションのメソッドはあるようだが、使い方がわからんので、最悪手打ちで

StepFive <- function(pbl, checkData){ # pbl:検査対象プロジェクト削除済みかつ多重線形性排除済みの全データ
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


StepFiveCor <- function(pbl, checkData){ # pbl:検査対象プロジェクト削除済みかつ多重線形性排除済みの全データ
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


StepFiveMulti <- function(pbl, checkData){ # pbl:検査対象プロジェクト削除済みかつ多重線形性排除済みの全データ
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


Cross <- function(pblData){ #pblDataは多重線形性排除済みの全プロジェクトデータ
　　ret <- 0; 

  for(i in 1:length(pblData[,1])){
    Study <- pblData[-(i),];
    Est <- pblData[i, ];
    ret <-  rbind(ret, StepFive(Study, Est));
  }
  ret <- ret[-1,]
  return(ret)
}

CrossCor <- function(pblData){ #pblDataは多重線形性排除済みの全プロジェクトデータ
　　ret <- 0; 

  for(i in 1:length(pblData[,1])){
    Study <- pblData[-(i),];
    Est <- pblData[i, ];
    ret <-  rbind(ret, StepFiveCor(Study, Est));
  }
  ret <- ret[-1,]
  return(ret)
}


CrossMulti <- function(pblData){ #pblDataは多重線形性排除済みの全プロジェクトデータ
　　ret <- 0; 

  for(i in 1:length(pblData[,1])){
    Study <- pblData[-(i),];
    Est <- pblData[i, ];
    ret <-  rbind(ret, StepFiveMulti(Study, Est));
  }
  ret <- ret[-1,]
  return(ret)
}


Choice <- function(pbl){ #多重線形性排除済みの全データ
　　Data <- pbl; Data[is.na(Data)] <- 0 #欠損値を排除
  ret <- 0 #返り値用の変数

  for(i in 5:length(Data)){
    ret <- c(ret, cor((Data$Actual -  Data$Estimated), Data[, i])) #retに工数と要因間の相関係数を代入
  }
  
  ret <- ret[-1]　#初期化の時の0を削除
  names(ret) <- names(pbl[1,5:length(pbl)]) #ラベル付

  sortList <- order(-abs(ret)) #大きい順にソート
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

ChoiceMKtwo <- function(pbl){ #多重線形性排除済みの全データ
　　Data <- pbl; Data[is.na(Data)] <- 0 #欠損値を排除
  ret <- 0 #返り値用の変数

  for(i in 5:length(Data)){
    ret <- c(ret, cor((Data$Actual -  Data$Estimated), Data[, i])) #retに工数と要因間の相関係数を代入
  }
  
  ret <- ret[-1]　#初期化の時の0を削除
  names(ret) <- names(pbl[1,5:length(pbl)]) #ラベル付

  sortList <- order(-abs(ret)) #大きい順にソート

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
  SortedData <- Choice(PBL) #工数変動要因を工数誤差との相関係数の高い順にソート
  ChoiceCoeff <- ChoiceMKtwo(PBL) #各変動要因の相関係数を代入
  
  multico <- Multico(SortedData) #多重線形性を排除する
  
#  write.table(ChoiceCoeff, file = "output.txt", append = TRUE, quote = FALSE);
 write.csv(ChoiceCoeff, file = "/Users/saitotakeru/Documents/Study/workspace/work1/SubData/MetrixInfluence.csv", quote = FALSE, col.names = FALSE);
  write.csv(multico, file = "/Users/saitotakeru/Documents/Study/workspace/work1/SubData/UseMetrix.csv", quote = FALSE, col.names = FALSE)
  MeanRes <- 0
  MedianRes <- 0
  VarRes <- 0
  for(i in 6:length(multico)){ #見積もり工数の誤差をresultに代入
    result <- CrossMulti(multico[, 1:i])[, 3]
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
  multico <- Multico(SortedData) #多重線形性を排除する
#  ret <- StepFiveMulti(multico[-length(PBL[,1]),],multico[length(PBL[,1]),])
  num <- as.integer(names(error))+4
  cat("num = ", num, ", multico = ", length(multico[1,]) )
  Study <- multico[-(i),1:num];
  Est <- multico[i, 1:num];
  ret <-  StepFiveMulti(Study, Est);
  
#  names(ret) <- c("Estimate", "Actual","Error" , "Alfa", "FPTrial", "predict", "Sigma");
#  write.table(ret, file = "output.txt", append = TRUE, quote = FALSE);
#  write.table(ret, file = "SubData/HyperParameters.csv", append = TRUE, quote = FALSE)
  names(ret) <- c(NULL, NULL,NULL , NULL, NULL, NULL, NULL);
  AandB <- 0
  AandB <- rbind(AandB, ret)
  write.table(AandB[-1,], file = "/Users/saitotakeru/Documents/Study/workspace/work1/SubData/HyperParameters.csv", append = TRUE, quote = FALSE)
  
options(scipen=5);return (ret[1]  - error)
#  return(ret[,3] - error)
#  return (length(multico))
#  return (ret[3]- error)
}

#MakeModel(PBLData)
#CalcManHour(PBLData,3)
#system.time(CalcManHour(PBLData))
m1 <- 0
for(i in 1:length(PBLData[,1])){
#  write.table(label[i], file = "output.txt", append = TRUE, quote = FALSE);
  m1 <- rbind(m1 , CalcManHour(PBLData,i))
}
m1 <- m1[-1]

plot(m, xaxt="n", pch = 20, xlab="projects", ylab="Error of scheduled man-hours and actual man-hours");
axis(side = 1, at = 1:length(m), labels = label)
text(1:length(m), m - 50, ceiling(m))