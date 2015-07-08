# Rotina para reproduzir os resultados do artigo 
# "Avaliando as medidas de núcleo da inflação no Brasil"
#  
# revised: Jul 7, 2014
#
# This script requires the following packages:
# FitAR, fUnitRoots, urca, sandwich, xtable and forecast
# make sure you install these packages before you load them.
#
#R Core Team (2014). R: A language and environment for statistical
#computing. R Foundation for Statistical Computing, Vienna, Austria. URL
#http://www.R-project.org/.
#


########################
### importando dados ###
########################

dados.artigo <- read.csv("dados/dados-artigo-atualizado.csv", sep=";")
ipca <- ts(dados.artigo[,c("IPCA")], star=c(1996,1), freq=12)
ms <- ts(dados.artigo[,c("IPCA.MS")], star=c(1996,1), freq=12)
ex2 <- ts(dados.artigo[,c("IPCA.EX2")], star=c(1996,1), freq=12)
dp <- ts(dados.artigo[,c("IPCA.DP")], star=c(1996,1), freq=12)

#############################################
### Figura 2: Núcleos da Inflação e IPCA. ###
#############################################

# imprime o grafico 1
temp <- range(ipca)
pdf(file="ex2.pdf", width=10, height=2.8)
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(1.5,2.5,1,1), cex.main=1, cex.axis=1.5, family="serif", font.main=1)
plot(ipca, lty=2, ylab="", xlab="", ylim=temp,  xaxt="n")
lines(ex2, lty=1, lwd=2)
axis(1,at=c(2000, 2005, 2010), label=c("","",""))
legend("topright", inset=.00, c("IPCA", "IPCA-EX2"), lty=c(2, 1), lwd=c(1, 2), bty="n", cex=1.5)
dev.off()

# imprime o grafico 2
temp <- range(ipca)
pdf(file="ms.pdf", width=10, height=2.8)
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(1.5,2.5,1,1), cex.main=1, cex.axis=1.5, family="serif", font.main=1)
plot(ipca, lty=2, ylab="", xlab="", ylim=temp, xaxt="n")
axis(1,at=c(2000, 2005, 2010), label=c("","",""))
lines(ms, lty=1, lwd=2)
legend("topright", inset=.00, c("IPCA", "IPCA-MS"), lty=c(2, 1), lwd=c(1, 2), bty="n", cex=1.5)
dev.off()

# imprime o grafico 3
temp <- range(ipca)
pdf(file="dp.pdf", width=10, height=3)
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(2.5,2.5,1,1), cex.main=1, cex.axis=1.5, family="serif", font.main=1)
plot(ipca, lty=2, ylab="", xlab="", ylim=temp)
lines(dp, lty=1, lwd=2)
legend("topright", inset=.00, c("IPCA", "IPCA-DP"), lty=c(2, 1), lwd=c(1, 2), bty="n", cex=1.5)
dev.off()

###########################################################
### Tabela 2: Resultados dos testes de estacionariedade ###
###########################################################

#Teste ADF
adf<-function(serie, tipo="c"){ # tipo é "nc", "c", "ct"
  require(FitAR)
  require(fUnitRoots)
  # retira os NAs
  serie<-na.omit(serie) 
  # seleciona o lag de acordo com o criterio de AIC
  lag<-SelectModel(serie, lag.max=12, Criterion="AIC", ARModel="AR", Best=1)
  # aplica o teste ADF aumentado com base em McKinnon(1996)
  teste<-unitrootTest(serie, lags=lag, type=tipo)
  return(teste)
}

#Teste KPSS
kpss<-function(serie,tipo){ # tipo é "mu" ou "tau"
  require(urca)
  require(sandwich)
  serie<-na.omit(serie)
  # seleciona automaticamente o lag com base em Newey & West (1994)
  fit<-lm(serie~1)
  lag<-floor(bwNeweyWest(fit, kernel="Bartlett", prewhite=0))
  # KPSS
  teste<-ur.kpss(serie, use.lag=lag, type=tipo)
  return(teste)
}

#Construcao da tabela
names <- c("IPCA", "IPCA-EX2", "IPCA-MS", "IPCA-DP")
teste.adf <- c(adf(ipca, tipo="c")@test$statistic, adf(ex2, tipo="c")@test$statistic, adf(ms, tipo="c")@test$statistic, adf(dp, tipo="c")@test$statistic )
teste.adf.p <-  c(adf(ipca, tipo="c")@test$p.value[1], adf(ex2, tipo="c")@test$p.value[1], adf(ms, tipo="c")@test$p.value[1], adf(dp, tipo="c")@test$p.value[1])
teste.kpss <- c(kpss(ipca, tipo="mu")@teststat, kpss(ex2, tipo="mu")@teststat, kpss(ms, tipo="mu")@teststat, kpss(dp, tipo="mu")@teststat)
teste.adf.i <- c(adf(window(ipca, end=c(2009,12)), tipo="c")@test$statistic, adf(window(ex2, end=c(2009,12)), tipo="c")@test$statistic, adf(window(ms, end=c(2009,12)), tipo="c")@test$statistic, adf(window(dp, end=c(2009,12)), tipo="c")@test$statistic )
teste.adf.p.i <-  c(adf(window(ipca, end=c(2009,12)), tipo="c")@test$p.value[1], adf(window(ex2, end=c(2009,12)), tipo="c")@test$p.value[1], adf(window(ms, end=c(2009,12)), tipo="c")@test$p.value[1], adf(window(dp, end=c(2009,12)), tipo="c")@test$p.value[1])
teste.kpss.i <- c(kpss(window(ipca, end=c(2009,12)), tipo="mu")@teststat, kpss(window(ex2, end=c(2009,12)), tipo="mu")@teststat, kpss(window(ms, end=c(2009,12)), tipo="mu")@teststat, kpss(window(dp, end=c(2009,12)), tipo="mu")@teststat)
tabela <- cbind(teste.adf.i, teste.adf.p.i,  teste.kpss.i, teste.adf, teste.adf.p,  teste.kpss)
tabela <- cbind(teste.adf.i,  teste.kpss.i, teste.adf,  teste.kpss)
rownames(tabela) <- names
require(xtable)
xtable(tabela)


#################################################################################
### Tabela 3: Resultado do teste F para ausência de viés da medidas de núcleo ###
#################################################################################

# A function to do "nested" F-tests
f.test.lm <- function(R.lm, F.lm) {
  # R.lm: modelo lm restrito
  # F.lm: modelo lm irrestrito
  SSE.R <- sum(resid(R.lm)^2)
  SSE.F <- sum(resid(F.lm)^2)
  df.num <- R.lm$df - F.lm$df
  df.den <- F.lm$df
  F <- ((SSE.R - SSE.F) / df.num) / (SSE.F / df.den)
  p.value <- 1 - pf(F, df.num, df.den)
  return(data.frame(F, df.num, df.den, p.value))
}

## teste de não viesada
tabela <- matrix(NA,3,1)
rest<- lm(ipca ~ 0 + offset(1*ex2))
irest <- lm(ipca ~ ex2)
tabela[1,1] <- f.test.lm(R.lm=rest, F.lm=irest)$p.value

rest<- lm(ipca ~ 0 + offset(1*ms))
irest <- lm(ipca ~ ms)
tabela[2,1] <- f.test.lm(R.lm=rest, F.lm=irest)$p.value

rest<- lm(ipca ~ 0 + offset(1*dp))
irest <- lm(ipca ~ dp)
tabela[3,1] <- f.test.lm(R.lm=rest, F.lm=irest)$p.value

require(xtable)
xtable(tabela, digits=4)


###################################################################################################
### Tabela 4: Resultados do teste t para os coeficientes de ajustamento da inflação e do núlceo ###
###################################################################################################

require(lmtest)
require(sandwich)
tabela<-matrix(NA, 6,8)

for(h in c(3,6,9,12)){
  dados1 <- ts.intersect(lag(ipca,h) - ipca, (ipca- ex2), lag(ipca,-1), lag(ipca,-2), lag(ipca,-3), lag(ipca,-4), lag(ipca,-5), lag(ipca,-6), dframe = TRUE)
  dados2 <- ts.intersect(lag(ex2,h) - ex2, (ipca- ex2), lag(ex2,-1), lag(ex2,-2), lag(ex2,-3), lag(ex2,-4), lag(ex2,-5), lag(ex2,-6), dframe = TRUE)
  v1 <- names(dados1)
  v2 <- names(dados2)
  aic.best1<-aic.best2<-Inf
  for(i in 1:6){
    fit1<-lm(as.formula(paste(v1[1], " ~ ", paste(v1[2:(2+i)], collapse = " + "))), data=dados1, na.action = NULL)
    aic1<-AIC(fit1)
    fit2<-lm(as.formula(paste(v2[1], " ~ ", paste(v2[2:(2+i)], collapse = " + "))), data=dados2, na.action = NULL)
    aic2<-AIC(fit2)
    if(aic1<aic.best1){
      fit.best1<-fit1
      aic.best1<-aic1
    }
    if(aic2<aic.best2){
      fit.best2<-fit2
      aic.best2<-aic2
    }
  }
  tabela[1,2*(h/3)-1] <- coeftest(fit.best1)[2,1]
  tabela[2,2*(h/3)-1] <- coeftest(fit.best1)[2,4]
  tabela[1,2*(h/3)] <- coeftest(fit.best2)[2,1]
  tabela[2,2*(h/3)] <- coeftest(fit.best2)[2,4]
}

for(h in c(3,6,9,12)){
  dados1 <- ts.intersect(lag(ipca,h) - ipca, (ipca- ms), lag(ipca,-1), lag(ipca,-2), lag(ipca,-3), lag(ipca,-4), lag(ipca,-5), lag(ipca,-6), dframe = TRUE)
  dados2 <- ts.intersect(lag(ms,h) - ms, (ipca- ms), lag(ms,-1), lag(ms,-2), lag(ms,-3), lag(ms,-4), lag(ms,-5), lag(ms,-6), dframe = TRUE)
  v1 <- names(dados1)
  v2 <- names(dados2)
  aic.best1<-aic.best2<-Inf
  for(i in 1:6){
    fit1<-lm(as.formula(paste(v1[1], " ~ ", paste(v1[2:(2+i)], collapse = " + "))), data=dados1, na.action = NULL)
    aic1<-AIC(fit1)
    fit2<-lm(as.formula(paste(v2[1], " ~ ", paste(v2[2:(2+i)], collapse = " + "))), data=dados2, na.action = NULL)
    aic2<-AIC(fit2)
    if(aic1<aic.best1){
      fit.best1<-fit1
      aic.best1<-aic1
    }
    if(aic2<aic.best2){
      fit.best2<-fit2
      aic.best2<-aic2
    }
  }
  tabela[3,2*(h/3)-1] <- coeftest(fit.best1)[2,1]
  tabela[4,2*(h/3)-1] <- coeftest(fit.best1)[2,4]
  tabela[3,2*(h/3)] <- coeftest(fit.best2)[2,1]
  tabela[4,2*(h/3)] <- coeftest(fit.best2)[2,4]
}

for(h in c(3,6,9,12)){
  dados1 <- ts.intersect(lag(ipca,h) - ipca, (ipca - dp), lag(ipca,-1), lag(ipca,-2), lag(ipca,-3), lag(ipca,-4), lag(ipca,-5), lag(ipca,-6), dframe = TRUE)
  dados2 <- ts.intersect(lag(dp,h) - dp, (ipca - dp), lag(dp,-1), lag(dp,-2), lag(dp,-3), lag(dp,-4), lag(dp,-5), lag(dp,-6), dframe = TRUE)
  v1 <- names(dados1)
  v2 <- names(dados2)
  aic.best1<-aic.best2<-Inf
  for(i in 1:6){
    fit1<-lm(as.formula(paste(v1[1], " ~ ", paste(v1[2:(2+i)], collapse = " + "))), data=dados1, na.action = NULL)
    aic1<-AIC(fit1)
    fit2<-lm(as.formula(paste(v2[1], " ~ ", paste(v2[2:(2+i)], collapse = " + "))), data=dados2, na.action = NULL)
    aic2<-AIC(fit2)
    if(aic1<aic.best1){
      fit.best1<-fit1
      aic.best1<-aic1
    }
    if(aic2<aic.best2){
      fit.best2<-fit2
      aic.best2<-aic2
    }
  }
  tabela[5,2*(h/3)-1] <- coeftest(fit.best1)[2,1]
  tabela[6,2*(h/3)-1] <- coeftest(fit.best1)[2,4]
  tabela[5,2*(h/3)] <- coeftest(fit.best1)[2,1]
  tabela[6,2*(h/3)] <- coeftest(fit.best1)[2,4]
}

require(xtable)
xtable(tabela, digits=3)


###################################################################################
### Figura 3: Variação percentual em 12 meses do IPCA e dos núcleos da inflação ###
###################################################################################

# funcao acumala ultimos 12 meses
acum<-function(x){
  date<-time(x)
  x12<- vector()
  for (i in 1:(length(date)-11)){
    x12[i] <- round((prod(1+(window(x, start=date[i], end=date[i+11]))/100)-1)*100,2)
  }
  x12<-ts(x12, start=date[12], freq=12)
  return(x12)
}

# transforma para acumalado dos 12 meses
ipca12<-acum(ipca)
ex212<-acum(ex2)
ms12<-acum(ms)
dp12<-acum(dp)

# imprime o grafico
temp <- range(ipca12)
pdf(file="nucleos.pdf", width=8, height=2.5)
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(2.5,2.5,2,1), cex.main=1, cex.axis=1.2, family="serif", font.main=1)
plot(ipca12, lwd=2, ylab="", xlab="", ylim=temp+c(0,0.2*temp[2]))
lines(ms12, lty=1, lwd=1)
lines(ex212, lty=2, lwd=1)
lines(dp12, lty=3, lwd=2)
legend("topright", inset=.00, c("IPCA", "IPCA-EX2", "IPCA-MS", "IPCA-DP"), lty=c(1, 1, 2, 3), lwd=c(2, 1, 1, 2), bty="n", cex=1)
dev.off()

################################################################################
### Tabela 5: REQM da previsão fora da amostra para o IPCA com o modelo ARMA ###
################################################################################

## Seleciona a ordem do ARMA pelo criterio AIC

order.arma<-function(serie){
  aic.best<-Inf
  for(i in 1:3){
    for(j in 1:3){
      fit <- arima(serie, order=c(i,0,j), method="CSS", optim.control=list(maxit=2000))
      raizMA <- sum(abs(fit$coef[(i+1):(i+j)]))
      raizAR <- sum(abs(fit$coef[1:i]))
      n<-length(serie)
      aic <- -2*fit$loglik + (2*(j+i+1)*n)/(n-j-i-2)
      if(aic<aic.best & raizMA<1 & raizAR<1){
        fit.best<-fit
        aic.best<-aic
        order<-c(i,0,j)
      }
    }
  }
  for(i in 1:3){
    fit <- arima(serie, order=c(i,0,0), method="CSS", optim.control=list(maxit=2000))
    aic <- -2*fit$loglik + (2*(i+1)*n)/(n-i-2)
    raizAR <- sum(abs(fit$coef[1:i]))
    if(aic<aic.best & raizAR<1){
      fit.best<-fit
      aic.best<-aic
      order<-c(i,0,0)
    }
  }
  for(j in 1:3){
    fit <- arima(serie, order=c(0,0,j), method="CSS", optim.control=list(maxit=2000))
    raizMA <- sum(abs(fit$coef[1:j]))
    aic <- -2*fit$loglik + (2*(j+1)*n)/(n-j-2)
    if(aic<aic.best & raizMA<1){
      fit.best <-fit
      aic.best<-aic
      order<-c(0,0,j)
    }
  }
  return(list(fit=fit.best, order=order))
}

## função que normaliza a matriz de previsao e erro para que se tenha o mesmo numero de previsoes para cada horizonte
trans <- function(erro,k,h){
  # erro: matriz de erros de previsao fora da amostra
  # h: horizonte das previsoes fora da amostra
  # k: numero de previsoes fora da amostraa
  # transforma a matriz para que cada linha corresponda ao mesmo periodo de previsao
  #
  ans<-erro[h:(k+h-1),1]
  for(i in 1:(h-1)){
    ans<-cbind(ans,erro[(h-i):(k+h-1-i),i+1])
  }
  return(ans)
}

## Prever a série com ARMA e compara com os valores da série ipca
ARMAoutsample <- function(serie, serie.ipca, k=36,  h=12){
  # serie: ts que será prevista com ARMA
  # ipca: ts usada para calcular os erros de previsão
  # h: horizonte das previsoes fora da amostra
  # k: numero de previsoes fora da amostraa
  #
  T  <- length(serie)
  date <- time(serie)
  erro <- prev <- matrix(NA, (k+h-1), h)
  j<-0
  ordem.arma <- matrix(NA,2,3)
  # função simples, realiza previsões com a amostra restrita
  for(i in 1:(k+h-1)){
    serie.amostra <- window(serie, end=date[T-k-h+i])
    model <- order.arma(serie.amostra)
    fcast <- predict(model$fit, n.ahead = 12)
    prev[i,] <- fcast$pred
    erro[i,] <- c(fcast$pred-serie.ipca,rep(NA,12-length(fcast$pred-serie.ipca)))
    if(i==1|i==(k+h-1)){
      # registra o arma escolhido no início e no fim da amostra
      j<-j+1
      ordem.arma[j,] <- model$order
    }
  }
  prev <- ts(trans(prev, h=h, k=k), start=date[T-k+1], freq=12)
  erro <- ts(trans(erro, h=h, k=k), start=date[T-k+1], freq=12)
  eqm <- apply(erro^2, 2, mean)
  eqm <- round(eqm[c(3,6,9,12)],3)
  reqm <- round(sqrt(eqm),3)
  return(list(erro=erro[,c(3,6,9,12)], reqm=reqm, eqm=eqm, prev=prev[,c(3,6,9,12)], ordem=ordem.arma))
}


ARMAipca24 <- ARMAoutsample(ipca12, ipca12, k=24)
ARMAex224 <- ARMAoutsample(ex212, ipca12, k=24)
ARMAms24 <- ARMAoutsample(ms12, ipca12, k=24)
ARMAdp24 <- ARMAoutsample(dp12, ipca12, k=24)

ARMAipca36 <- ARMAoutsample(ipca12, ipca12, k=36)
ARMAex236 <- ARMAoutsample(ex212, ipca12, k=36)
ARMAms36 <- ARMAoutsample(ms12, ipca12, k=36)
ARMAdp36 <- ARMAoutsample(dp12, ipca12, k=36)

ARMAipca48 <- ARMAoutsample(ipca12, ipca12, k=48)
ARMAex248 <- ARMAoutsample(ex212, ipca12, k=48)
ARMAms48 <- ARMAoutsample(ms12, ipca12, k=48)
ARMAdp48 <- ARMAoutsample(dp12, ipca12, k=48)

tabela <- rbind(rbind(ARMAipca24$reqm, ARMAex224$reqm, ARMAms24$reqm, ARMAdp24$reqm),
                rbind(ARMAipca36$reqm, ARMAex236$reqm, ARMAms36$reqm, ARMAdp36$reqm),
                rbind(ARMAipca48$reqm, ARMAex248$reqm, ARMAms48$reqm, ARMAdp48$reqm))

tabela <- rbind(rbind(ARMAipca24$eqm/ARMAipca24$eqm, ARMAex224$eqm/ARMAipca24$eqm, ARMAms24$eqm/ARMAipca24$eqm, ARMAdp24$eqm/ARMAipca24$eqm),
                rbind(ARMAipca36$eqm/ARMAipca36$eqm, ARMAex236$eqm/ARMAipca36$eqm, ARMAms36$eqm/ARMAipca36$eqm, ARMAdp36$eqm/ARMAipca36$eqm),
                rbind(ARMAipca48$eqm/ARMAipca48$eqm, ARMAex248$eqm/ARMAipca48$eqm, ARMAms48$eqm/ARMAipca48$eqm, ARMAdp48$eqm/ARMAipca48$eqm))

tabela <- rbind(rbind(ARMAipca24$reqm/ARMAipca24$reqm, ARMAex224$reqm/ARMAipca24$reqm, ARMAms24$reqm/ARMAipca24$reqm, ARMAdp24$reqm/ARMAipca24$reqm),
                rbind(ARMAipca36$reqm/ARMAipca36$reqm, ARMAex236$reqm/ARMAipca36$reqm, ARMAms36$reqm/ARMAipca36$reqm, ARMAdp36$reqm/ARMAipca36$reqm),
                rbind(ARMAipca48$reqm/ARMAipca48$reqm, ARMAex248$reqm/ARMAipca48$reqm, ARMAms48$reqm/ARMAipca48$reqm, ARMAdp48$reqm/ARMAipca48$reqm))

rownames(tabela) <-rep(c("IPCA 24", "IPCA-EX2 24", "IPCA-MS 24", "IPCA-DP 24"))
require(xtable)
xtable(tabela)



# calcula o p-valor do teste DM modificado

require(forecast)
tabela <- matrix(NA,9,4)
for(i in 1:4){
  tabela[1,i] <- dm.test(ARMAipca24$erro[,i], ARMAex224$erro[,i], h=3*i, alternative="two.sided", power=2)$p.value
  tabela[2,i] <- dm.test(ARMAipca24$erro[,i], ARMAms24$erro[,i], h=3*i, alternative="two.sided", power=2)$p.value
  tabela[3,i] <- dm.test(ARMAipca24$erro[,i], ARMAdp24$erro[,i], h=3*i, alternative="two.sided", power=2)$p.value
  tabela[4,i] <- dm.test(ARMAipca36$erro[,i], ARMAex236$erro[,i], h=3*i, alternative="two.sided", power=2)$p.value
  tabela[5,i] <- dm.test(ARMAipca36$erro[,i], ARMAms36$erro[,i], h=3*i, alternative="two.sided", power=2)$p.value
  tabela[6,i] <- dm.test(ARMAipca36$erro[,i], ARMAdp36$erro[,i], h=3*i, alternative="two.sided", power=2)$p.value
  tabela[7,i] <- dm.test(ARMAipca48$erro[,i], ARMAex248$erro[,i], h=3*i, alternative="two.sided", power=2)$p.value
  tabela[8,i] <- dm.test(ARMAipca48$erro[,i], ARMAms48$erro[,i], h=3*i, alternative="two.sided", power=2)$p.value
  tabela[9,i] <- dm.test(ARMAipca48$erro[,i], ARMAdp48$erro[,i], h=3*i, alternative="two.sided", power=2)$p.value
}

require(xtable)
xtable(tabela, digits=4)
