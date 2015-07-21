# Codigo em R para reproduzir os resultados do artigo 
# "Avaliando as medidas de núcleo da inflação no Brasil"
#  
# revised: Jul 7, 2014
#
# This script requires the following packages:
# FitAR, fUnitRoots, urca, sandwich, lmtest, xtable, dyn, zoo and car
# make sure you install these packages before you load them.
#



########################
### importa os dados ###
########################

dados.artigo <- read.csv("dados/dados-artigo-atualizado.csv", sep=";")
ipca <- ts(dados.artigo[,c("IPCA")], star=c(1996,1), freq=12)
ex <- ts(dados.artigo[,c("IPCA.EX")], star=c(1996,1), freq=12)
ex2 <- ts(dados.artigo[,c("IPCA.EX2")], star=c(1996,1), freq=12)
ms <- ts(dados.artigo[,c("IPCA.MS")], star=c(1996,1), freq=12)
dp <- ts(dados.artigo[,c("IPCA.DP")], star=c(1996,1), freq=12)
# inclui o nucleo 


#############################################
### Figura 2: Núcleos da Inflação e IPCA. ###
#############################################

# imprime o grafico 1
temp <- range(ipca)
pdf(file="ex.pdf", width=10, height=2.8)
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(1.5,2.5,1,1), cex.main=1, cex.axis=1.5, family="serif", font.main=1)
plot(ipca, lty=2, ylab="", xlab="", ylim=temp,  xaxt="n", bty="n")
abline(h=-0.65, v=c(1995.25), lwd=2.5)
lines(ex, lty=1, lwd=2)
axis(1,at=c(2000, 2005, 2010), label=c("","",""))
legend("topright", inset=.00, c("IPCA", "IPCA-EX"), lty=c(2, 1), lwd=c(1, 2), bty="n", cex=1.5)
dev.off()

# imprime o grafico 2
temp <- range(ipca)
pdf(file="ex2.pdf", width=10, height=2.8)
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(1.5,2.5,1,1), cex.main=1, cex.axis=1.5, family="serif", font.main=1)
plot(ipca, lty=2, ylab="", xlab="", ylim=temp,  xaxt="n",  bty="n")
abline(h=-0.65, v=c(1995.25), lwd=2.5)
lines(ex2, lty=1, lwd=2)
axis(1,at=c(2000, 2005, 2010), label=c("","",""))
legend("topright", inset=.00, c("IPCA", "IPCA-EX2"), lty=c(2, 1), lwd=c(1, 2), bty="n", cex=1.5)
dev.off()

# imprime o grafico 3
temp <- range(ipca)
pdf(file="ms.pdf", width=10, height=2.8)
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(1.5,2.5,1,1), cex.main=1, cex.axis=1.5, family="serif", font.main=1)
plot(ipca, lty=2, ylab="", xlab="", ylim=temp, xaxt="n",  bty="n")
abline(h=-0.65, v=c(1995.25), lwd=2.5)
axis(1,at=c(2000, 2005, 2010), label=c("","",""))
lines(ms, lty=1, lwd=2)
legend("topright", inset=.00, c("IPCA", "IPCA-MS"), lty=c(2, 1), lwd=c(1, 2), bty="n", cex=1.5)
dev.off()

# imprime o grafico 4
temp <- range(ipca)
pdf(file="dp.pdf", width=10, height=3)
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(2.5,2.5,1,1), cex.main=1, cex.axis=1.5, family="serif", font.main=1)
plot(ipca, lty=2, ylab="", xlab="", ylim=temp, bty="n")
abline(h=-0.65, v=c(1995.25), lwd=2.5)
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
names <- c("IPCA", "IPCA-EX", "IPCA-EX2", "IPCA-MS", "IPCA-DP")
teste.adf <- c(adf(ipca, tipo="c")@test$statistic, adf(ex, tipo="c")@test$statistic, adf(ex2, tipo="c")@test$statistic, adf(ms, tipo="c")@test$statistic, adf(dp, tipo="c")@test$statistic )
teste.adf.p <-  c(adf(ipca, tipo="c")@test$p.value[1], adf(ex, tipo="c")@test$p.value[1], adf(ex2, tipo="c")@test$p.value[1], adf(ms, tipo="c")@test$p.value[1], adf(dp, tipo="c")@test$p.value[1])
teste.kpss <- c(kpss(ipca, tipo="mu")@teststat, kpss(ex, tipo="mu")@teststat, kpss(ex2, tipo="mu")@teststat, kpss(ms, tipo="mu")@teststat, kpss(dp, tipo="mu")@teststat)
teste.adf.i <- c(adf(window(ipca, end=c(2009,12)), tipo="c")@test$statistic, adf(window(ex, end=c(2009,12)), tipo="c")@test$statistic, adf(window(ex2, end=c(2009,12)), tipo="c")@test$statistic, adf(window(ms, end=c(2009,12)), tipo="c")@test$statistic, adf(window(dp, end=c(2009,12)), tipo="c")@test$statistic )
teste.adf.p.i <-  c(adf(window(ipca, end=c(2009,12)), tipo="c")@test$p.value[1], adf(window(ex, end=c(2009,12)), tipo="c")@test$p.value[1],  adf(window(ex2, end=c(2009,12)), tipo="c")@test$p.value[1], adf(window(ms, end=c(2009,12)), tipo="c")@test$p.value[1], adf(window(dp, end=c(2009,12)), tipo="c")@test$p.value[1])
teste.kpss.i <- c(kpss(window(ipca, end=c(2009,12)), tipo="mu")@teststat, kpss(window(ex, end=c(2009,12)), tipo="mu")@teststat, kpss(window(ex2, end=c(2009,12)), tipo="mu")@teststat, kpss(window(ms, end=c(2009,12)), tipo="mu")@teststat, kpss(window(dp, end=c(2009,12)), tipo="mu")@teststat)
tabela <- cbind(teste.adf.i, teste.adf.p.i,  teste.kpss.i, teste.adf, teste.adf.p,  teste.kpss)
tabela <- cbind(teste.adf.i,  teste.kpss.i, teste.adf,  teste.kpss)
rownames(tabela) <- names

options(OutDec= ",") 

require(xtable)
xtable(tabela)


#################################################################################
### Tabela 3: Resultado do teste F para ausência de viés da medidas de núcleo ###
#################################################################################



## teste de não viesada
Ftest <- function(x){
  require(car)
  require(sandwich)
  hm <- rbind(c(1,0),c(0,1))
  rhs <- c(0,1)
  test <- linearHypothesis(model=x,hypothesis.matrix=hm, rhs=rhs, vcov=NeweyWest(x))
  p.value <- test[[4]][2]
  return(p.value)
}

tabela <- matrix(NA,4,3)

model <- lm(ipca ~ ex)
tabela[1,] <- c(coef(model),Ftest(model))

model <- lm(ipca ~ ex2)
tabela[2,] <- c(coef(model),Ftest(model))

model <- lm(ipca ~ ms)
tabela[3,] <- c(coef(model),Ftest(model))

model <- lm(ipca ~ dp)
tabela[4,] <- c(coef(model),Ftest(model))

rownames(tabela) <- c("ex", "ex2", "ms", "dp")
colnames(tabela) <- c("alpha", "beta", "valor p")
  
require(xtable)
xtable(tabela, digits=4)


###################################################################################################
### Tabela 4: Resultados do teste t para os coeficientes de ajustamento da inflação e do núlceo ###
###################################################################################################


tab.ajust <- function(x,ipca){
  #input: 
  # x: (ts) serie do  nucleo
  # ipca: (ts) serie do ipca
  #output: 
  # tabela: (matrix) tabela com coeficiente e valor p
  
  require(lmtest)
  require(sandwich)
  tabela<-matrix(NA, 2,8)
    for(h in c(3,6,9,12)){
      dados1 <- ts.intersect(lag(ipca,h) - ipca, (ipca- x), lag(ipca,-1), lag(ipca,-2), lag(ipca,-3), lag(ipca,-4), lag(ipca,-5), lag(ipca,-6), dframe = TRUE)
      dados2 <- ts.intersect(lag(x,h) - x, (ipca- x), lag(x,-1), lag(x,-2), lag(x,-3), lag(x,-4), lag(x,-5), lag(x,-6), dframe = TRUE)
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
      tabela[1,2*(h/3)-1] <- coeftest(fit.best1, vcov=vcovHC(fit.best1, type="HC0"))[2,1]
      tabela[2,2*(h/3)-1] <- coeftest(fit.best1, vcov=vcovHC(fit.best1, type="HC0"))[2,4]
      tabela[1,2*(h/3)] <- coeftest(fit.best2, vcov=vcovHC(fit.best1, type="HC0"))[2,1]
      tabela[2,2*(h/3)] <- coeftest(fit.best2, vcov=vcovHC(fit.best1, type="HC0"))[2,4]
    }
  return(tabela)
}


tabela <- rbind(tab.ajust(ex, ipca), tab.ajust(ex2, ipca), tab.ajust(ms, ipca), tab.ajust(dp, ipca))
tab <- tabela

for(i in 1:4){
  tab[2*i-1, ] <- format(round(tabela[2*i-1, ], 3), nsmall = 3)
  tab[2*i, ] <- paste("(",format(round(tabela[2*i, ], 3), nsmall = 3) , ")", sep="")
}


#Insere colunas em data.frame
insert.dt <- function(dt, x,i){
  #input:
  # dt(df):data.frame
  # x(vt): vector com tamnaho igual ao numero de colunas de dt
  # i(int): numero da linha na qual o vector sera incluido
  #output: 
  # dt(df): data.frame com o vector x incluido da linha i
  
  n <- length(dt[1,])
  if(i==1){
    dt<- cbind(as.list(x),dt)
  }else{
    if(i==n){
      dt<- cbind(dt, as.list(x))
    }else{
      dt <- cbind(dt[,1:(i-1)],  as.list(x),dt[,i:n])
    }
  }
  return(dt)
}

tab <- insert.dt(tab,rep("",8), 3)
tab <- insert.dt(tab,rep("",8), 6)
tab <- insert.dt(tab,rep("",8), 9)

col_name <- c("IPCA-EX", "", "IPCA-EX2", "", "IPCA-MS", "", "IPCA-DP", "")
tab <- cbind(col_name, tab)

require(xtable)
xtable(tab, digits=3)


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
ex12<-acum(ex)
ex212<-acum(ex2)
ms12<-acum(ms)
dp12<-acum(dp)

# imprime o grafico
temp <- range(ipca12)
pdf(file="nucleos.pdf", width=8, height=2.5)
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(2.5,2.5,1,1), cex.main=1, cex.axis=1, family="serif", font.main=1)
plot(ipca12, ylab="", xlab="", ylim=temp+c(-0.5*temp[1],0.2*temp[2]), bty="n")
abline(h = -0.01, v=c(1996.17), lwd=c(2.5, 2.5))
lines(ex12, lty=2, lwd=1)
lines(ex212, lty=3, lwd=1)
lines(ms12, lty=4, lwd=1)
lines(dp12, lty=5, lwd=1)
legend("topright", inset=.00, c("IPCA","IPCA-EX", "IPCA-EX2", "IPCA-MS", "IPCA-DP"), lty=c(1,2,3,4,5), lwd=c(1,1,1,1,1), bty="n", cex=0.8)
dev.off()


# imprime o grafico 4
temp <- range(ipca)
pdf(file="dp.pdf", width=10, height=3)
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(2.5,2.5,1,1), cex.main=1, cex.axis=1.5, family="serif", font.main=1)
plot(ipca, lty=2, ylab="", xlab="", ylim=temp, bty="n")
abline(h=-0.65, v=c(1995.25))
lines(dp, lty=1, lwd=2)
legend("topright", inset=.00, c("IPCA", "IPCA-DP"), lty=c(2, 1), lwd=c(1, 2), bty="n", cex=1.5)
dev.off()

######################################################################################
### Tabela 5: REQM obtido com os núcleos relativo ao REQM do modelo de referência  ###
######################################################################################


## Constroi dummies sazonais 
SeasDummy <- function(wts, s, t0, type){
  #input: 
  # wts: (ts) serie temporal
  # s: (int) frequencia da serie
  # t0: (vector) data inicial da serie
  # type: (character) dominio do tempo "alg" ou dominio da frequencia "trg"
  #output: 
  # VFE: (ts) serie com dummies sazonais
  
  
  N <- length(wts)
  if(type == "alg"){        # Empleada en Barsky & Miron (1989)
    auxD <- matrix(0,nrow=N, ncol=s)
    sq   <- seq(1,N,s)
    k    <- 0
    
    for(j in 1:s){
      n <- ifelse(sq[length(sq)] + k > N, length(sq)-1,  N)
      for(i in 1:n)
        auxD[sq[i]+k,j] <- 1
      k <- k+1
    }
    VFE <- auxD
    if(t0[2] != 1){
      VFE <- matrix(nrow=N, ncol=s)
      VFE[,1:(t0[2]-1)] <- auxD[,(s-t0[2]+2):s]; VFE[,t0[2]:s] <- auxD[,1:(s-t0[2]+1)]
    }
    if(t0[2] == 1){ VFE <- auxD }
  }
  
  if(type == "trg"){        # Empleada en Granger & Newbold (1986)
    qq  <- s/2
    VFE <- matrix(nrow=N, ncol=(s-1))
    
    sq1 <- seq(1,qq*2,2)
    sq2 <- seq(2,qq*2,2)
    j   <- c(1:(qq-1))
    
    for(i in 1:N){
      for(k in 1:(s-qq-1)){
        tmp <- i * j[k] * pi / qq
        VFE[i,sq1[k]] <- cos(tmp)
        VFE[i,sq2[k]] <- sin(tmp)
      }
      VFE[i,(s-1)] <- (-1)^i
    }
  }
  VFE
}



## Regressao com series temporais
tsreg <- function(y, z, p, x=NULL, q, h, s=FALSE){
  #input: 
  # y: (ts) variavel a ser prevista  
  # z: (ts) variavel independente 1
  # x: (ts) variavel independente 2
  # p: (int) numero de defasagens de z
  # q: (int) numero de defasagens de x
  # h: (int) horizonte de previsao 
  # s: (logical) TRUE insere dummies sazonais  
  #output: 
  # best.fit: o modelo selecionado (dyn)
  
  require(dyn)
  require(zoo)
  aux <- ts(rep(NA,h),start=tsp(y)[2]+1/12 , freq=tsp(y)[3])
  y<-ts.union(y,aux)[,1]
  L <- function(x, k = 1) lag(x, -k)
  SDum <- SeasDummy(y, s=frequency(y), t0=start(y), type="alg")
  colnames(SDum)<-as.character(seq(1:length(SDum[1,])))
  if(is.null(x)){ 
    dados <- as.zoo(cbind(y,z,s=SDum))
    v <- colnames(dados)
    if(s==TRUE){
      fit <- dyn$lm(as.formula(paste("y ~ 0 + L(z, (h):(h+p-1)) +", paste(v[-(1:2)], collapse=" + ", sep=""), sep="")) , data=dados)
    }else{
      fit <- dyn$lm(y ~ L(z, (h):(h+p-1)), data=dados)
    } 
  }else{
    x<-ts.union(x,aux)[,1]
    dados <- as.zoo(cbind(y, z, x, v=SDum))
    v <- colnames(dados)
    if(s==TRUE){
      fit <- dyn$lm(as.formula(paste("y ~ 0 + L(z, (h):(h+p-1)) +  L(x, (h):(h+q-1)) + ", paste(v[-(1:3)], collapse=" + ", sep=""), sep="")) , data=dados)
    }else{
      fit <- dyn$lm(y ~ L(z, (h):(h+p-1)) +  L(x, (h):(h+q-1)), data=dados)
    } 
  }
  # get prediction for h 
  fcast <- tail(na.omit(predict(fit, dados)), 1)
  return(list(fcast=fcast, fit=fit))
}


## Seleciona e estima o modelo tsreg com menor BIC
tsreg.select <- function(y, z, p, x=NULL, q, h, s=FALSE){
  #input: 
  # y: (ts) variavel a ser prevista
  # z: (ts) variavel independente 1
  # x: (ts) variavel independente 2
  # p: (int) numero máximo de defasagens de z
  # q: (int) numero máximo de defasagens de x
  # h: (int) horizonte de previsao 
  # s: (logical) TRUE insere dummies sazonais  
  #output: 
  # fit: (dyn) o modelo selecionado
  # order: (vector) ordem selecionada
  # fcast: (zoo) previsao selecionada
  
  bic.best<-Inf
  if(is.null(x)){
    for(i in 1:p){
      model <- tsreg(y=y, z=z, p=i, h=h, s=s)
      bic <- BIC(model$fit)
      if(bic<bic.best){
        bic.best <- bic
        order.best <- c(i,0)
        fit.best <- model$fit
        fcast.best <- model$fcast
      }
    }
  }else{
    for(i in 1:p){
      for(j in 1:q){
        model <- tsreg(y=y, z=z, p=i, x=x, q=j, h=h, s=s)
        bic <- BIC(model$fit)
        if(bic<bic.best){
          bic.best <- bic
          order.best <- c(i,j)
          fit.best <- model$fit
          fcast.best <- model$fcast
        }
      }
    }
  }
  return(list(fit=fit.best, order=order.best, fcast=fcast.best))
}


# Realiza as previsões fora da amostra
outsample <- function(y, z, x=NULL, p,  q,  k, s=FALSE){
  #input: 
  # y: (ts) variavel a ser prevista
  # z: (ts) variavel independente 1
  # x: (ts) variavel independente 2
  # p: (int) numero máximo de defasagens de z
  # q: (int) numero máximo de defasagens de x
  # h: (int) horizonte de previsao 
  # s: (logical) TRUE insere dummies sazonais  
  #output: 
  # erro(zoo): erros de previsao 
  # prev (zoo): valores previstos
  # eqm(num): erro quadratiico medio
  # reqm(num): raiz do erro quadratiico medio
  # model (dyn): ultimo modelo selecionado
  
  T  <- length(y)
  date <- time(y)
  erro <- prev <- matrix(NA, k, 4)
  # função simples, realiza previsões com a amostra restrita
  list.model<-vector("list")
  for(h in c(3,6,9,12)){
    for(i in 1:k){
      y.amostra <- window(y, end=date[T-k-h+i])
      z.amostra <- window(z, end=date[T-k-h+i])
      x.amostra <- x
      if(!is.null(x)){
        x.amostra <- window(x, end=date[T-k-h+i])
      }
      # estima o modelo VAR por OLS
      model <- tsreg.select(y=y.amostra, z=z.amostra, x=x.amostra, p=p, q=q, h=h, s=s)
      # prever as variaveis do VAR h meses a frente
      prev[i,(h/3)] <- model$fcast
      erro[i,(h/3)] <- model$fcast-as.zoo(y)
      if(i==k){
        list.model[[paste(h)]] <- model
      }
    }
  }
  prev <- ts(prev, end=index(model$fcast), freq=12)
  erro <- ts(erro, end=index(model$fcast), freq=12)
  eqm <- apply(erro^2, 2, mean)
  reqm <- sqrt(eqm)
  return(list(erro=erro, reqm=reqm, eqm=eqm, prev=prev, model=list.model))
}


# Realiza as simulacoes de previsao fora da amosta
output_ex <- output_ex2 <- output_dp <- output_ms <- output_ipca <- vector("list")
for(i in c(24, 36, 48)){
  output_ipca[[paste(i)]] <- outsample(y=ipca12, z=ipca12, p=4, q=0, k=i)
  output_ex[[paste(i)]] <- outsample(y=ipca12, z=ex12,  x=ipca12, p=4, q=4, k=i)
  output_ex2[[paste(i)]] <- outsample(y=ipca12, z=ex212, x=ipca12, p=4, q=4, k=i)
  output_ms[[paste(i)]] <- outsample(y=ipca12, z=ms12, x=ipca12, p=4, q=4, k=i)
  output_dp[[paste(i)]] <- outsample(y=ipca12, z=dp12, x=ipca12, p=4, q=4, k=i)
}




# Tabela com o REQM relativo
tabela1 <- matrix(NA, 12, 4)
id <- c(1,5,9)
for(i in 1:3){
  tabela1[id[i],] <- output_ex[[i]]$reqm / output_ipca[[i]]$reqm
  tabela1[id[i]+1,] <- output_ex2[[i]]$reqm / output_ipca[[i]]$reqm
  tabela1[id[i]+2,] <- output_ms[[i]]$reqm / output_ipca[[i]]$reqm
  tabela1[id[i]+3,] <- output_dp[[i]]$reqm / output_ipca[[i]]$reqm
}



# muda o separador decimal de . para ,
options(OutDec= ".") 

# Formata a tabela
tabela <- matrix(NA, 12, 4)
for(i in 1:12){
  tabela[i,] <- format(round(tabela1[i,],digits = 2), nsmall = 2)
}


# nomeando as colunas e linhas da tabela
id_row <- rep(c("IPCA-EX", "IPCA-EX2","IPCA-MS", "IPCA-DP"),3)
tabela <- cbind(id_row, tabela)
id_col <- c("modelo", "h=3", "h=6", "h=9", "h=12")
tabela <- rbind(id_col,tabela)

tab <- data.frame(tabela)
require(xtable)
xtable(tab)


#######################################################################################
### Apendice - Tabela 6: Especificação dos modelos de previsão no final da amostra  ###
#######################################################################################


tabela2 <- matrix(NA, 4, 5)
for(i in 1:4){
  tabela2[i,1] <- paste("p=",output_ipca[[1]]$model[[i]]$order[1], ", q=", output_ipca[[1]]$model[[i]]$order[2], sep="")
  tabela2[i,2] <- paste("p=",output_ex[[1]]$model[[i]]$order[1], ", q=", output_ex[[1]]$model[[i]]$order[2], sep="")
  tabela2[i,3] <- paste("p=",output_ex2[[1]]$model[[i]]$order[1], ", q=", output_ex2[[1]]$model[[i]]$order[2], sep="")
  tabela2[i,4] <- paste("p=",output_ms[[1]]$model[[i]]$order[1], ", q=", output_ms[[1]]$model[[i]]$order[2], sep="")
  tabela2[i,5] <- paste("p=",output_dp[[1]]$model[[i]]$order[1], ", q=", output_dp[[1]]$model[[i]]$order[2], sep="")
}

id_col <- c("IPCA","EX", "EX2","MS", "DP") 
tabela2 <- rbind(id_col,tabela2)

id_row <- c("serie", "h=3", "h=6", "h=9", "h=12")
tabela2 <- cbind(id_row, tabela2)

require(xtable)
xtable(tabela2)


