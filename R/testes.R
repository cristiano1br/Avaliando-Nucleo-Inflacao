# Project: Avaliando Nucleo Inflacao
# Testes 
# Script usado para testar as funções deste projeto
#
#

ls()

dados.artigo <- read.csv("dados/dados-artigo-atualizado.csv", sep=";")
ipca <- ts(dados.artigo[,c("IPCA")], star=c(1996,1), freq=12)
ms <- ts(dados.artigo[,c("IPCA.MS")], star=c(1996,1), freq=12)
ex2 <- ts(dados.artigo[,c("IPCA.EX2")], star=c(1996,1), freq=12)
dp <- ts(dados.artigo[,c("IPCA.DP")], star=c(1996,1), freq=12)

# funcao acumala ultimos 12 meses
acum<-function(x){
  date<-time(x)
  x12<- vector()
  for (i in 1:(length(date)-11)){
    x12[i] <- (prod(1+(window(x, start=date[i], end=date[i+11]))/100)-1)*100
  }
  x12<-ts(x12, start=date[12], freq=12)
  return(x12)
}

# transforma para acumalado dos 12 meses
ipca12<-acum(ipca)
ex212<-acum(ex2)
ms12<-acum(ms)
dp12<-acum(dp)

y=ipca12
x=ms12

plot(y)
plot(x)


source('C:/Users/Cristiano/OneDrive/RStudio/generalfunctions.R')
teste <- tsreg(y=y, p=3, h=12, s=TRUE)
teste$fcast
summary(teste$fit)


teste <- tsreg(y=y, p=3, x=x,q=3, h=12, s=TRUE)
teste$fcast
summary(teste$fit)
BIC(teste$fit)


teste <- tsreg.select(y=y, p=3, x=x,q=3, h=12, s=TRUE)

system.time(teste1 <- outsample(y=y, p=3, h=12, m=36))
system.time(teste2 <- outsample(y=y, p=3, x=x, q=3, h=12, m=36))


VARfit<-VAR(cbind(ipca,ms), p=6, type="const", season=12)
