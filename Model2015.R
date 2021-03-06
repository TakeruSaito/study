#作業手順
#１．規模の比例定数αを求める
#αを求めるため、最小二乗法を用いる
#最小二乗法は関数lsfit(x,y)を用いる。　x：説明変数の行列　　y；被説明変数のベクトル
#f <- lsfit(予測工数, 規模)
#f$coefficients で表示される　Intercept  V1がy = ax + bのb,aに該当する。

CalcAlpha <- function(pbl){ #csvの元データそのまま

  TriAlfa <- lsfit(pbl$Estimated, pbl$FPTrial)$coeff
  AppAlfa <- lsfit(pbl$Estimated, pbl$FPApproxi)$coeff 

  return (TriAlfa)
}

#２．各変動要因間の相関性を調べ、多重線形性を排除する
#source("http://aoki2.si.gunma-u.ac.jp/R/src/tolerance.R", encoding="euc-jp")
#をRコンソールにペーストしたのち、tolerance(x)　x：説明変数だけのデータ行列　を使う
#tolerance(x)の出力結果のうち、toleranceの値が大きいものを削除する

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

CalcCoeff <- function(SelectedData, type = "multi"){
  factor <- 0
  data <- SelectedData
  data[is.na(data)] <- 0
  
  if(type == "multi" || is.na(type)){
    
    library("QuantPsyc", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
    factor <- lm.beta(lm((data$Actual - data$Estimated)?., data[-1:-4]))
    
  }else{
    for(i in 5:length(data)){
      if(type == "regression"){
        factor <- rbind(factor, lm((Actual-Estimated)?data[,i], data)$coefficients)
      }else if(type == "cor"){
        factor <- rbind(factor, cor(data$Actual - data$Estimated, data[,i]))
      }else if(type == "t"){
        factor <- rbind(factor, t.test(data$Actual - data$Estimated, data[,i], paired = TRUE)$statistic)
      }
    }
    factor <- factor[-1,]
  }
  factor[is.na(factor)] <- 0
  if(type == "regression"){
    factor <- factor[,2]
    kname <- names(data)[5:length(data)]
    names(factor) <- kname
  }else if(type == "cor" || type == "t"){
    sname <- names(data)[5:length(data)]
    names(factor) <- sname
  }
  return (factor)
}

#４．最小二乗法を用いて、見積もりモデル式を作成する
#現在想定しているモデル式は(工数)　= α　×　FP(規模)　×　(１＋　β　×　Σ{(変動要因)　×　(各変動要因の回帰係数)})
#なので、(工数)- α×FP　=　α　×　FP×　β　×　(変動要因とその回帰係数の積　の総和)と式変形し、
#をlsfit(αFP(変動要因とその回帰係数の積　の総和), (工数-αFP))とし、１.と同様にして、係数βを算出する

#(工数)- α×FP = data$Actual - (TriAlfa[2] * data$FPTrial)

CalcBeta <- function(data, TriAlfa, factor){
  
  FactMulCoeff <- 0
  for(i in 5:length(data)){
    FactMulCoeff <- c(FactMulCoeff, c(data[,i] * factor[i - 4]))
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
  return( lm(leftFormula?rightFormula, data=x)$coeff )
}

#５．クロスバリデーションで推定精度を算出する。
#Rでのクロスバリデーションのメソッドはあるようだが、使い方がわからんので、最悪手打ちで

CalcPriMH <- function(pbl, checkData, type = "multi", mode = "default"){
  # pbl:検査対象プロジェクト削除済みかつ多重線形性排除済みの全データ
  alpha <- CalcAlpha(pbl)
  
  pbl <- Choice(pbl)
  
  factor <- abs(CalcCoeff(pbl, type))
  
  if(mode == "debug"){
    write.table(rownames(checkData), file = "/Users/saitotakeru/Documents/Study/workspace/work1/SubData/Factors.csv", append = TRUE, quote = FALSE,sep = ",")
    write.table(t(factor), file = "/Users/saitotakeru/Documents/Study/workspace/work1/SubData/Factors.csv", append = TRUE, quote = FALSE,sep = ",")
  }
  predict <- CalcBeta(pbl, alpha, factor)
  
  FactMulCoeff <- 0
  
  for(i in 5:length(pbl)){
    FactMulCoeff <- c(FactMulCoeff, c(checkData[i] * factor[i - 4]))
  }

  FactMulCoeff <- FactMulCoeff[-1]
  sigma <- 0
  for(i in 1:length(FactMulCoeff)){
    sigma <- sigma + as.numeric(FactMulCoeff[i])
  }
  
  ret <- alpha [2] * checkData$FPTrial * (1 + predict[2] * sigma)
  ret <- abs(ret)
  ret <- c(ret, checkData$Actual, ret - checkData$Actual, alpha[2], checkData$FPTrial, predict[2], sigma)
  return(ret) 
}

CrossValid <- function(pblData, type = "multi"){
  ret <- 0;
  
  for(i in 1:length(pblData[,1])){
    Study <- pblData[-(i),];
    Est <- pblData[i,];
    ret <- rbind(ret, CalcPriMH(Study, Est, type))
  }
  ret <- ret[-1,]
  return (ret)
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
  return(ret)
}

ChoiceValues <- function(pbl){ #多重線形性排除済みの全データ
　　Data <- pbl; Data[is.na(Data)] <- 0 #欠損値を排除
  ret <- 0 #返り値用の変数

  for(i in 5:length(Data)){
    ret <- c(ret, cor((Data$Actual -  Data$Estimated), Data[, i])) #retに工数と要因間の相関係数を代入
  }
  
  ret <- ret[-1]　#初期化の時の0を削除
  names(ret) <- names(pbl[1,5:length(pbl)]) #ラベル付

  sortList <- order(-abs(ret)) #大きい順にソート

  return(ret)
}

Bunsan <- function(Data, type = "multi"){
  ret <- 0
  for(i in 6:length(Data)){
    ret <- c(ret , sqrt(variance(CrossValid(Data[,1:i])[,3], type)))
  }
  #ret <- ret[-1]
  return (ret)
}

variance <- function(x) var(x)*(length(x)-1)/length(x)

MakeModel <- function(PBL, type = "multi"){
  SortedData <- Choice(PBL) #工数変動要因を工数誤差との相関係数の高い順にソート
  
  ChoiceCoeff <- ChoiceValues(PBL) #各変動要因の相関係数を代入
  
  multico <- Multico(SortedData) #多重線形性を排除する
  
#  write.table(ChoiceCoeff, file = "output.txt", append = TRUE, quote = FALSE);

#  write.csv(ChoiceCoeff, file = "/Users/saitotakeru/Documents/Study/workspace/work1/SubData/MetrixInfluence.csv", quote = FALSE, col.names = FALSE);
#  write.csv(multico, file = "/Users/saitotakeru/Documents/Study/workspace/work1/SubData/UseMetrix.csv", quote = FALSE, col.names = FALSE)
  MeanRes <- 0
  MedianRes <- 0
  VarRes <- 0
  for(i in 6:length(multico)){ #見積もり工数の誤差をresultに代入
    result <- CrossValid(multico[, 1:i], type)[, 3]
    
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

CalcManHour <- function(PBL,i, type = "multi"){
#  error <- MakeModel(PBL[-length(PBL[,1]),])
  
  error <- MakeModel(PBL, type)
  SortedData <- Choice(PBL)
  
  multico <- Multico(SortedData) #多重線形性を排除する
#  ret <- StepFiveMulti(multico[-length(PBL[,1]),],multico[length(PBL[,1]),])
  num <- as.integer(names(error))+4
  
  Study <- multico[-(i),1:num];
  Est <- multico[i, 1:num];
  
  ret <-  CalcPriMH(Study, Est, type, mode = "debug");
  
  names(ret) <- c("Estimate", "Actual","Error" , "Alfa", "FPTrial", "predict", "Sigma");
  AandB <- 0
  AandB <- rbind(AandB, ret)
  #writeLines(paste(i), "/Users/saitotakeru/Documents/Study/workspace/work1/SubData/HyperParameters.csv")
  #write.table(t(AandB[-1,]), file = "/Users/saitotakeru/Documents/Study/workspace/work1/SubData/HyperParameters.csv", append = TRUE, quote = FALSE,sep = ",")
  
  options(scipen=5);return (ret[1]  - error)
}

CalcRelativeError <- function(experiment, calculated){
  return ( ((experiment - calculated) / calculated) * 100 )  
}

CalcTestData <- function(PBLData, type = "multi"){
  ret <- 0
  for(i in 1:length(PBLData[,1])){
    ret <- rbind(ret , CalcManHour(PBLData, i, type) )
    cat( "loop : ", i, "?n")
  }
  ret <- ret[-1]
  return(ret)
}


m1 <- 0
for(i in 1:length(PBLData[,1])){
#  write.table(label[i], file = "output.txt", append = TRUE, quote = FALSE);
  m1 <- rbind(m1 , CalcManHour(PBLData,i))
}
m1 <- m1[-1]

plot(m, xaxt="n", pch = 20, xlab="projects", ylab="Error of scheduled man-hours and actual man-hours");
axis(side = 1, at = 1:length(m), labels = label)
text(1:length(m), m - 50, ceiling(m))

m2 <- 0
for(i in 1:length(PBLData[,1])){
  #  write.table(label[i], file = "output.txt", append = TRUE, quote = FALSE);
  m2 <- rbind(m2 , CalcManHour(PBLData,i, type = "t") )
}
m2 <- m2[-1]