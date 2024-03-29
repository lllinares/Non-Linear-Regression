library("rstudioapi")
getwd()
setwd(dirname(getActiveDocumentContext()$path))       # Set working directory to source file location
getwd()
scalebiaised<-function(x){ #centrer et r�duire
  ce<-x-mean(x)
  re<-sqrt(var(x)*(length(x)-1)/length(x))
  return(ce/re)
}

Predict<-function(y,Z,lambda,M,x){ #fait la r�gression lin�aire des moindres carr� p�nalis�e ou non et pr�dit un nouveau point
  beta_hat<-solve(1/length(y)*t(Z)%*%Z+lambda*M)%*%t(Z)%*%cbind(y)*1/length(y) #l'estimation des beta ridge
  pred=rbind(x)%*%cbind(beta_hat) #pr�diction d'un y pour un vecteur x donn?
  return(list("beta_hat"=beta_hat,"predi"=pred)) #on renvoie les deux
}

#mod�le lin�aire
#########################################################################################################################
loose<-function(Z,y,lambda,M,afhist){ #fait le leave one out pour la r�gression lin�aire Ridge
  erreurquadr=c() #initialisation des erreurs quadratiques
  beta_isuc=list() #les betas successives
  ytilde_isuc=list() #les estimations des y_i successives
  for (i in 1:length(y)){ #pour chaque �l�ment de la r�ponse y
    y_i=y[-i] #on enl�ve le i�me
    Z_i=Z[-i,] #on enl�ve la i?me ligne des variables explicatives correspondante
    pred=Predict(y_i,Z_i,lambda,M,Z[i,]) #utilise la fonction Predict
    beta_isuc=append(beta_isuc,list(pred$beta_hat)) #on rajoute les beta � la liste
    ytilde_isuc=append(ytilde_isuc,list(pred$predi)) #on rajoute les y_i � la liste
    erreurquadr_i=(y[i]-pred$predi)^2 # on calcule l'�cart entre le vrai y_i et la pr�diction
    erreurquadr[i]=erreurquadr_i  # on rajoute l'�cart � la liste
  }
  erreurquadrmoy=mean(erreurquadr) #on fait la moyenne de tout ces �carts
  if (!missing(afhist)){
    hist(erreurquadr,main=NULL,ylab=NULL,xlab=NULL) #on affiche ou non l'histogramme
  }
  return(list("beta_isuc"=beta_isuc,"ytilde_isuc"=ytilde_isuc,"erreurquadr"=erreurquadr,"erreurquadrmoy"=erreurquadrmoy)) #on renvoie tout ces ?l?ments
}

set.seed(20210127) #simulation d'un mod?le lin�aire
n <- 100
x_test1 <- scalebiaised(runif(n, min = -2, max = 2)) 
x_test2 <- scalebiaised(runif(n, min = -2, max = 2))
beta_0 <- -2
beta_1 <- 3
beta_2<-1.5
epsilon_test <- rnorm(n, mean = 0, sd = 1)
Z_test<-cbind(rep(1,n),x_test1,x_test2)
y_test <- beta_0 + beta_1 * x_test1 +beta_2*x_test2+epsilon_test
b<-loose(Z_test,y_test,0,diag(1,3),afhist="la") #on obtient un histogramme des erreurs quadratiques
b 
title(main="Erreurs quadratiques de pr�diction ",sub="pour une r�gression lin�aire sur donn�es simul�es",
      ylab="Fr�quence",xlab="erreurs quadratiques")
#########################################################################################################################


Zind<-function(x,N){  #repartit les nombres de x dans les intervalles et sort la matrice des indicatrices
  inter=c()
  if (!(length(x)%%N==0)) #si la longueur de x n'est pas divisible par N le nombre de tranches on a une mauvaise valeur
  {
    print("mauvaise valeur")
    return()
  }
  Tr=length(x)/N #le nombre d'�l�ment dans une tranche
  xs=x#on fait une copie de x
  for (i in  1:N){ #pour chaque num�ro de tranche
    for (j in 1:Tr){ #pour chaque �l�ment de la tranche
      ind=which.min(xs) #on cherche le minimum et son indice
      xs[ind]=max(xs)+1 #on remplace cet �l�ment par un autre beaucoup trop grand
      inter[ind]=i #on affecte cet emplacement au num�ro de la tranche
    }
  }
  d<-sort(x) #on range les x
  tuple=list() #on initialise une liste de couple
  tuple=append(tuple,list(rbind(d[1],(d[Tr]+d[Tr+1])/2)))# le premier couple c'est la premi�re valeur et la demi-somme du  bout de la tranche avec la suivante
  for (i in 1:(N-2)){
    tuple=append(tuple,list(rbind((d[i*Tr]+d[i*Tr+1])/2,(d[(i+1)*Tr]+d[(i+1)*Tr+1])/2))) #on recalcule les demi-sommes pour qu'il y ait correspondance entre l'intervalle et la classe
  }
  tuple=append(tuple,list(rbind((d[Tr*(N-1)]+d[Tr*(N-1)+1])/2,d[length(x)])))
  
  return(list("matrice"=model.matrix(~as.factor(inter)-1),"codage"=cbind(x,inter),"intervalle"=tuple)) #on renvoie une matrice d'indicatrice suivant les modalit?s du num?ro de la tranche les intervalles et le codage des x
  
}
RidgeM<-function(y,Z,lambda,M){ #r�sout la r�gression lin�aire p�nalis� et sort les estimations des param�tres 
  beta_hat<-solve((1/length(y)*t(Z)%*%Z+lambda*M))%*%t(Z)%*%cbind(y)*1/length(y)  #l'estimation des beta ridge
  y_hat<-Z%*%cbind(beta_hat) #l'estimation des y_i
  return(list("beta_hat"=beta_hat,"y_hat"=y_hat)) #on renvoie ces deux vecteurs
}


div<-function(n){ #les diviseurs sauf 1 et le nombre lui-m?me
  d=c()
  for (i in 2:(n-1)){
    if (n%%i==00){
      d=c(d,i)
    }}
  return(d)
}
ind <- function(x,a,b)  ifelse(x >= a & x <= b, 1,0) #une indicatrice

estimbetadim<-function(X,y,N,lambda,M,x){# permet de prendre en compte le nombre de pr?dicteurs
  Z<-c()
  for (i in 1:dim(X)[2]){ #on regarde chaque colonne de X et on en fait une matrice d'indicatrice avec le partage en N tranche
    Z_i<-Zind(cbind(X)[,i],N)$matrice
    Z<-cbind(Z,Z_i) #on les colle ensemble
  }
  Z<-as.matrix(Z)
  
  return(list("beta_hat"=RidgeM(y,Z,lambda,M)$beta_hat,"y_hat"=RidgeM(y,Z,lambda,M)$y_hat)) #on utilise Ridge avec cete matrice
}

affichereg<-function(x,y,lambda,M,N){ #permet d'afficher la r�gression p�nalis�e avec indicatrice
  Z<-Zind(x,N)$matrice #on fait le partage de x
  Ma=metriqueIetD(N) #on calcule les m�triques pour un partage en 50
  
  if(M=="I"){Ma=Ma$I}
  if(M=="D"){Ma=Ma$D}
  if(M=="E"){Ma=Ma$E}
  beta<-RidgeM(y,Z,lambda,Ma)$beta_hat #on extrait le beta_hat
  plot(x,y,sub=paste(" R�gression sur indicatrices avec la m�trique",M,"pour un coefficient de p�nalit� de",as.character(lambda))) #on affiche les points (x,y) 
  alpha=1/N
  ii<-Zind(x,N)$intervalle #on extrait les intervalles
  lia=rep(0,N*2)
  for (j in 1:N){ #permet de construire la fonction en escalier
    lia[2*j-1]<-ii[[j]][1] #on prend les extr�mit�s des intervalles
    lia[2*j]<-ii[[j]][2]
  }
  lio=c()
  for (j in 1:N){
    lio<-c(lio,rep(beta[j],2)) # on prend les valeurs beta de l'escalier
  }
  lines(lia,lio,col='blue') #on l'affiche en lignee
}
metriqueIetD<-function(N,d=1){ #renvoie les m�triques d'int�r�ts
  N1=N*d #multiplie par la dimension
  if (missing(d)){
    N1=N
  }
  M1=diag(1,N1) #l'identit?
  M2demi<-diag(1,N1)
  M2demi[row(M2demi)-col(M2demi)==-1]<--1
  M2demi=rbind(M2demi[-N1,]) #m�trique des �carts successifs
  M2<-t(M2demi)%*%M2demi
  M3demi<-diag(1,N1) 
  M3demi[row(M3demi)-col(M3demi)==-1]<--2
  M3demi[row(M3demi)-col(M3demi)==-2]<-1
  M3demi<-M3demi[-N*d,]
  M3demi<-M3demi[-(N*d-1),]
  M3<-t(M3demi)%*%M3demi #m�trique des �carts des �carts
  
  
  return(list("I"=M1,"D"=M2,"E"=M3))
}

par(mfrow=c(1,1))
par(mar=c(5.1, 4.1, 4.1, 2.1))


library(readODS)
dataset <- data.matrix(read_ods("Datasim1.ods")) 
x_test<-scalebiaised(dataset[,2])
y_test<-dataset[,1]
plot(x_test, y_test)

N=50
affichereg(x_test,y_test,0.001,"I",N) #(sur-apprentissage) (effet limit� de la p�nalit�)
title(main="R�gression de y sur x avec une tr�s faible p�nalit�")
affichereg(x_test,y_test,0.01,"I",N) #ne fait qu'�craser la courbe vers 0 (sous-apprentissage)
title(main="R�gression de y sur x avec une faible p�nalit�")
affichereg(x_test,y_test,1,"I",N) #ne fait qu'�craser la courbe vers 0 (sous-apprentissage)
title(main="R�gression de y sur x avec une forte p�nalit�")



affichereg(x_test,y_test,0.001,"D",N) #(sur-apprentissage) (effet limit� de la p�nalit�)
title(main="R�gression de y sur x avec une tr�s faible p�nalit�")
affichereg(x_test,y_test,0.01,"D",N)
title(main="R�gression de y sur x avec une faible p�nalit�")
affichereg(x_test,y_test,1,"D",N)#affine la courbe en l'?crasant vers le centre (sous-apprentissage)
title(main="R�gression de y sur x avec une forte p�nalit�")


affichereg(x_test,y_test,0.001,"E",N) #(sur-apprentissage) (effet limit� de la p�nalit�)
title(main="R�gression de y sur x avec une tr�s faible p�nalit�")
affichereg(x_test,y_test,0.01,"E",N)#pas de d�faut majeur lisse la fonction en escalier
title(main="R�gression de y sur x avec une faible p�nalit�")
affichereg(x_test,y_test,1,"E",N)#pas de d?faut majeur
title(main="R�gression de y sur x avec une forte p�nalit�")


loola<-function(x,y,lambda,choix,alt,afhist){ #permet d'avoir les erreur moyennes pour plusieurs alpha et lambda possibles
  alpha=1/div(length(x)) #on prend les diviseurs de la longueur de x
  if (!missing(alt)){
    alpha=alt #on a le choix de s�lectionner un alpha
  }
  errormoy<-list()
  for(alp in alpha){ #pour chaque partage de x
    for(lam in lambda){ #pour chaque p�nalit�
      if (choix=="I"){ #je choisis la m�trique
        error<-loose(Zind(x,1/alp)$matrice,y,lam,metriqueIetD(1/alp)$I,afhist)$erreurquadrmoy #je fais le leave one out apr�s avoir partag� les x
      }
      if(choix=="D"){
        
        error<-loose(Zind(x,1/alp)$matrice,y,lam,metriqueIetD(1/alp)$D,afhist)$erreurquadrmoy
      }
      if (choix=="E"){
        error<-loose(Zind(x,1/alp)$matrice,y,lam,metriqueIetD(1/alp)$E,afhist)$erreurquadrmoy
        
      }
      errormoy<-append(errormoy,list(c(error,lam,alp)))}}
  return(errormoy)} 
alt=1/50
lambda=seq(0,0.1,0.01)
loola(x_test,y_test,lambda,"I",alt) #on cherche pour la m�trique I (erreur,p�nalit�,alpha)
print("lambda=(0),  erreur = 0.00616")
affichereg(x_test,y_test,0,"I",N) #sans p�nalit� 
title(main="R�gression de y sur x sans p�nalit�")


lambda=seq(0,0.1,0.01)
loola(x_test,y_test,lambda,"I",alt)
print("lambda=(0) erreur 0.00616")


lambda=seq(0,0.001,0.0001)
loola(x_test,y_test,lambda,"I",alt)
print("lambda=(0.0002) erreur 0.006120794")


affichereg(x_test,y_test,0.0002,"I",N) #p�nalit� 0.02 obtenu apr�s t�tonnement
title(main="R�gression de y sur x avec une p�nalit� optimale") #erreur de LOO :|0.006106342|

loofit<-function(x,y){
  error=c()
  if (!(is.null(dim(x)[2]))){ #si x correspond � plusieurs variables explicatives
    for( i in 1:dim(x)[1]){
      fit=lm(y[-i]~x[-i,]) #on enl�ve l'individu i
      error[i]=(fit$coefficients[1]+rbind(fit$coefficients[-1])%*%cbind(x[i,])-y[i])^2 #on calcule l'erreur entre la pr?diction et la vraie valeur
    }
    erreurquadrmoy=mean(error)
    return(erreurquadrmoy)}
  
  for(i in 1:length(x)){ #pour chaque �l�ment de x
    fit=lm(y[-i]~x[-i]) #on fait la r�gression avec un �l�ment en moins
    error[i]=(fit$coefficients[1]+fit$coefficients[2]*x[i]-y[i])^2 #on calcule le l'erreur entre la vraie valeur et la pr�diction
  }
  erreurquadrmoy=mean(error)
  return(erreurquadrmoy)
  
}
fit=lm(y_test~x_test)
loofit(x_test,y_test) #erreur de LOO : |0.00627|
abline(fit$coefficients[1],fit$coefficients[2],col="red") #comparaison avec une r�gression lin�aire simple
legend(-1,0.6,legend=c("indicatrices","lin�aire"),col=c("blue","red"),lty=1:2,cex=0.8)
loola(x_test,y_test,0.0002,"I",alt,T)
title(main="Erreurs quadratiques de pr�diction",sub="Pour une r�gression sur indicatrices avec la m�trique I et de coefficient de p�nalit� 0.0002",ylab="Fr�quence",xlab="erreurs quadratiques")


lambda=seq(0,5,1)
loola(x_test,y_test,lambda,"D",alt) #on cherche avec la m�trique D
print("lambda=1, erreur =0.004483")

lambda=seq(0.08,0.12,0.01)
loola(x_test,y_test,lambda,"D",alt)
print("lambda=(0.1), erreur =0.003102166")

lambda=seq(0.095,0.105,0.001)
loola(x_test,y_test,lambda,"D",alt)
print("lambda=(0.097), erreur =0.00310205")


affichereg(x_test,y_test,0.97,"D",N) #avec p�nalit� 9.6
title(main="R�gression de y sur x avec une p�nalit� optimale")#erreur LOO: |0.0031020|
abline(fit$coefficients[1],fit$coefficients[2],col="red")  #erreurLoo: |0.00627|
legend(-1,0.6,legend=c("indicatrices","lin�aire"),col=c("blue","red"),lty=1:2,cex=0.8)
loola(x_test,y_test,0.97,"D",alt,T)
title(main="Erreurs quadratiques de pr�diction",sub="Pour une r�gression sur indicatrices avec la m�trique D et de coefficient de p�nalit� 0.97",ylab="Fr�quence",xlab="erreurs quadratiques")



lambda=seq(0,100,10)
loola(x_test,y_test,lambda,"E",alt) #on cherche avec la m�trique E
print("lambda=(10),aplha=1/50 erreur =0.002922605")
affichereg(x_test,y_test,10,"E",N) 
lambda=seq(9,12,0.1)
lambda=11.7
loola(x_test,y_test,lambda,"E",alt)
print("lambda=(11.7),aplha=1/50 erreur =0.002922241")
affichereg(x_test,y_test,11.7,"E",N) #avec p�nalit� 11.7
title(main="R�gression de y sur x avec une p�nalit� optimale")#erreur LOO: |0.002922241|
abline(fit$coefficients[1],fit$coefficients[2],col="red")  ##erreurLoo: |0.00627|
legend(-1,0.6,legend=c("indicatrices","lin�aire"),col=c("blue","red"),lty=1:2,cex=0.8)
loola(x_test,y_test,11.7,"E",alt,T)
title(main="Erreurs quadratiques de pr�diction",sub="Pour une r�gression sur indicatrices avec la m�trique E et de coefficient de p�nalit� 11.7",ylab="Fr�quence",xlab="erreurs quadratiques")



library(readODS)
dataset <- data.matrix(read_ods("Datasim1.ods",sheet=2)) #on teste avec la simulation avec deux variables explicatives
print(dataset)
x_test1<-scalebiaised(dataset[,2])
x_test2<-scalebiaised(dataset[,3])
y_test<-dataset[,1]
N=25
M=metriqueIetD(N,2)
beta<-estimbetadim(cbind(x_test1,x_test2),y_test,N,1,M$D)$beta_hat #attention aux syst�mes singuliers
alpha=1/N
ii<-Zind(x_test1,N)$intervalle
lia=rep(0,N*2)
for (j in 1:N){ #permet de construire la fonction en escalier
  lia[2*j-1]<-ii[[j]][1]
  lia[2*j]<-ii[[j]][2]
}
lio=c()
for (j in 1:N){
  lio<-c(lio,rep(beta[j],2))
}
plot(lia,lio,col='blue',type='l',xlab="x1",ylab="beta1",main="Coefficients des indicatrices sur les tranches de la premi�re variable",sub="Pour la m�trique D et valeur de la p�nalit� 1")

ii<-Zind(x_test2,N)$intervalle
lia=rep(0,N*2)
for (j in 1:N){ #permet de construire la fonction en escalier
  lia[2*j-1]<-ii[[j]][1]
  lia[2*j]<-ii[[j]][2]
}
lio=c()
for (j in 1:N){
  lio<-c(lio,rep(beta[j+N],2))
}
plot(lia,lio,col='blue',type='l',xlab="x2",ylab="beta2",main="Coefficients des indicatrices sur les tranches de la seconde variable",sub="Pour la m�trique D et valeur de la p�nalit� 1")

plot(y_test,estimbetadim(cbind(x_test1,x_test2),y_test,N,1,M$D)$y_hat,ylab="y_hat",main="Correspondance entre les y observ�s 
     et estim�s pour la m�trique D et p�nalit� de coefficient 1") #correspondance y_hat  et y_test ?
abline(0,1,col="blue")
legend(1,1.6,legend="bissectrice",col="blue",lty=1:2,cex=0.8)
plot(y_test,estimbetadim(cbind(x_test1,x_test2),y_test,N,1,M$I)$y_hat,ylab="y_hat",main="Correspondance entre les y observ�s 
     et estim�s pour la m�trique I et p�nalit� de coefficient 1") #correspondance y_hat  et y_test ?
abline(0,1,col="blue")
legend(1,1.6,legend="bissectrice",col="blue",lty=1:2,cex=0.8)
plot(y_test,estimbetadim(cbind(x_test1,x_test2),y_test,N,1,M$E)$y_hat,ylab="y_hat",main="Correspondance entre les y observ�s 
     et estim�s pour la m�trique E et p�nalit� de coefficient 1") #correspondance y_hat  et y_test ?
abline(0,1,col="blue")
legend(1,1.6,legend="bissectrice",col="blue",lty=1:2,cex=0.8)



loolm<-function(x,y,lambda,alp){ #compare les m�triques
  li<-list()
  li<-append(li,list("I"=loola(x,y,lambda,"I",alp)))
  li<-append(li,list("D"=loola(x,y,lambda,"D",alp)))
  li<-append(li,list("E"=loola(x,y,lambda,"E",alp)))
  
  return(li)
}
dataset <- data.matrix(read_ods("Datasim1.ods"))
x_test<-scalebiaised(dataset[,2])
y_test<-dataset[,1]

loolm(x_test,y_test,0.01,0.02) #avec la p�nalit� 1 , il semble que E soit beaucoup mieux

dataset <- as.matrix(read.csv("Food_Supply_Quantity_kg_Data.csv",sep=',',header=T,row.names=1))
dim(dataset)
dataset
N=83
obesity=as.double(dataset[,24]) #pourcentage d'ob�sit� par pays
obesity=obesity[-c(110,53,148)]
meat=scalebiaised(as.double((dataset[,9])))#quantit� de viande
meat=meat[-c(110,53,148)]
meat=meat[-c(1)]
obesity=obesity[-c(1)]
affichereg(meat,obesity,127,"E",N)
title(main="R�gression de la pr�valence de l'ob�sit� sur la consommation de viande",xlab="Consommation de viande",ylab="Pourcentage d'ob�sit�")
h=seq(100,200,1)
h=127
loola(meat,obesity,h,"E",1/83)
print("erreur=62.41674489 lambda= (127)")
loola(meat,obesity,0,"E",1/83,T)
title(main="Erreurs quadratiques de pr�diction",sub="Pour une r�gression sur indicatrices avec la m�trique E et de coefficient de p�nalit� 127",ylab="Fr�quence",xlab="erreurs quadratiques")



affichereg(meat,obesity,0.6,"D",N)
title(main="R�gression de la pr�valence de l'ob�sit� sur la consommation de viande",xlab="Consommation de viande",ylab="Pourcentage d'ob�sit�")
h=seq(0,1,0.1)
loola(meat,obesity,h,"D",1/83)
print('erreur=63.55188252 lambda=(0.6)')
loola(meat,obesity,0.6,"D",1/83,T)
title(main="Erreurs quadratiques de pr�diction",sub="Pour une r�gression sur indicatrices avec la m�trique D et de coefficient de p�nalit� 0.6",ylab="Fr�quence",xlab="erreurs quadratiques")


cereal=as.double((dataset[,5]))#quantit� de c�r�ales 
cereal=scalebiaised(cereal[-c(110,53,148)])
cereal=cereal[-c(1)]
affichereg(cereal,obesity,597,"E",N)
title(main="R�gression de la pr�valence de l'ob�sit� sur la consommation de c�r�ales",xlab="Consommation de c�r�ales",ylab="Pourcentage d'ob�sit�")
h=seq(580,600,1)
loola(cereal,obesity,h,"E",1/83)
print("erreur=67.20206989 lambda= 597")
loola(cereal,obesity,597,"E",1/83,T)
title(main="Erreurs quadratiques de pr�diction",sub="Pour une r�gression sur indicatrices avec la m�trique E et de coefficient de p�nalit� 597",ylab="Fr�quence",xlab="erreurs quadratiques")


affichereg(cereal,obesity,0.8,"D",N)
title(main="R�gression de la pr�valence de l'ob�sit� sur la consommation de c�r�ales",xlab="Consommation de c�r�ales",ylab="Pourcentage d'ob�sit�")
h=seq(0,1,0.1)
loola(cereal,obesity,h,"D",1/83)
print("erreur=67.887 lambda= 0.8")
loola(cereal,obesity,0.8,"D",1/83,T)
title(main="Erreurs quadratiques de pr�diction",sub="Pour une r�gression sur indicatrices avec la m�trique D et de coefficient de p�nalit� 0.8",ylab="Fr�quence",xlab="erreurs quadratiques")


loofit(cereal,obesity) #erreur 66.47949
loofit(meat,obesity) #erreur 68.82583

#Attention on a enlev� un �l�ment car 167 est premier.

kernelgaussian<-function(x,h,p){  #noyau gaussien
  if (missing(p)){
    return(1/(h*sqrt(2*pi))*exp(-(x**2)/(2*h**2))) #loi normale unidimensionnel  ?l?ment par ?l?ment
  }
  return(1/((h*sqrt(2*pi))**p)*exp(-(t(x)%*%diag(rep(1/h**2,length(x)))%*%x)/2)) #loi normale multidimensionnelle pour un vecteur
  
}
kernelgaussian(cbind(c(1,2,0)),2,T)
kernelgaussian(1,3)
nadarayawatson<-function(z,x,y,h){#nadaraya-watson
  num<-0
  den<-0
  for (i in 1:length(x)){
    num<-num + y[i]* kernelgaussian(x[i]-z,h) #num�rateur somme pond�r�
    den<-den +kernelgaussian(x[i]-z,h) #divis� par la pond�ration �l�ment par �lement de z
  }
  return(num/den)
}
library(readODS)
dataset <- data.matrix(read_ods("Datasim1.ods"))
x_test<-scalebiaised(dataset[,2])
y_test<-dataset[,1]
titre="R�gression de nadaraya-watson pour un noyau d'�cart-type "
h=0.05 
plot(x_test,y_test,sub=paste(titre,h),main="R�gression de y sur x ")
z=seq(-2,2,0.01)
y_hat<-nadarayawatson(z,x_test,y_test,h)
a<-data.frame(z,y_hat)
a<-a[order(z),]
lines(z,a$y_hat,col="blue") #surajustement

h=0.1
plot(x_test,y_test,sub=paste(titre,h),main="R�gression de y sur x ")
z=seq(-2,2,0.01)
y_hat<-nadarayawatson(z,x_test,y_test,h)
a<-data.frame(z,y_hat)
a<-a[order(z),]
lines(z,a$y_hat,col="blue")

h=0.3
plot(x_test,y_test,sub=paste(titre,h),main="R�gression de y sur x ")
z=seq(-2,2,0.01)
y_hat<-nadarayawatson(z,x_test,y_test,h)
a<-data.frame(z,y_hat)
a<-a[order(z),]
lines(z,a$y_hat,col="blue")

h=3
plot(x_test,y_test,sub=paste(titre,h),main="R�gression de y sur x ")
z=seq(-2,2,0.01)
y_hat<-nadarayawatson(z,x_test,y_test,h)
a<-data.frame(z,y_hat)
a<-a[order(z),]
lines(z,a$y_hat,col="blue") #devient une droite



loonw<-function(x,y,h,afhist){ #fait le leave one out 
  erreurquadr=c()
  pred_suc=list()
  ytilde_isuc=list()
  for (i in 1:length(y)){
    y_i=y[-i]
    x_i=x[-i]
    pred=nadarayawatson(x[i],x_i,y_i,h) #on pr?dit sur l'?l?ment enlev?
    pred_suc=append(pred_suc,list(pred))
    erreurquadr_i=(y[i]-pred)^2
    erreurquadr[i]=erreurquadr_i  # calcule les erreurs successives de la pr?diction
  }
  erreurquadrmoy=mean(erreurquadr)
  if (!missing(afhist)){
    hist(erreurquadr,main=NULL,ylab=NULL,xlab=NULL) #on affiche ou non l'histogramme
  }
  return(list("pred_suc"=pred_suc,"erreurquadr"=erreurquadr,"erreurquadrmoy"=erreurquadrmoy))
}

loonw(x_test,y_test,0.1)

loonwch<-function(x,y,hach,afhist){ #permet d'avoir les erreur moyennes pour plusieurs h possibles
  errormoy<-list()
  for(h in hach){
    error<-loonw(x,y,h,afhist)$erreurquadrmoy
    errormoy<-append(errormoy,list(c(error,h)))}
  
  return(errormoy)} 
h=seq(0.1,5,0.1)
loonwch(x_test,y_test,h)
print("h=(0.2) erreur= 0.003114466")
h=seq(0.1,0.3,0.01)
loonwch(x_test,y_test,h)
print("h=(0.17) erreur= 0.003100911")
h=seq(0.1,5,0.1)
loonwch(x_test,y_test,h)
print("h=0.174 erreur= 0.003100559")

h=0.174
plot(x_test,y_test,sub=paste(titre,h),main="R�gression de y sur x optimale")
z=seq(-2,2,0.01)
y_hat<-nadarayawatson(z,x_test,y_test,h)
a<-data.frame(z,y_hat)
a<-a[order(z),]
lines(z,a$y_hat,col="blue")
fit=lm(y_test~x_test)
loofit(x_test,y_test) #erreur 0.00627 
abline(fit$coefficients[1],fit$coefficients[2],col="red")
legend(-1,0.6,legend=c("nadaraya-watson","r�gression lin�aire"),col=c("blue","red"),lty=1:2,cex=0.8)
loonwch(x_test,y_test,0.174,T)
title(main="Erreurs quadratiques de pr�diction ",sub="pour une r�gression � noyau d'�cart-type 0.174 ",
      ylab="Fr�quence",xlab="erreurs quadratiques")


matecpuis<-function(X,x,p) { #matrice pour approximation locale polynomiale
  z=c()
  for (j in 0:p){
    z=cbind(z,(X-x)**j) # X c'est un vecteur d'observations d'une variables explicatives et x une valeur d'int?r?t p degr? du polyn?me
  }
  
  return(z)
}

x_test<-scalebiaised(dataset[,2])
x<-2
matecpuis(x_test,x,3)

regpolynoy<-function(z,x,y,h,p){ 
  y_hat=c()
  for(i in 1:length(z)){ #un vecteur de z ? pr?dire
    X=matecpuis(x,z[i],p)
    beta<-solve(t(X)%*%diag(kernelgaussian(x-z[i],h))%*%X)%*%t(X)%*%diag(kernelgaussian(x-z[i],h))%*%y #c'est le beta solution du programme
    y_hat[i]=beta[1] #l'estimation de y c'est le coefficient constant 
  }
  return(y_hat)
  
}

h=0.174
p=2
titre=paste("R�gression polynomiale � noyau de degr�",p,"et d'�cart-type",h)
plot(x_test,y_test,sub=paste(titre),main="R�gression de y sur x ")
y_hat<-regpolynoy(z,x_test,y_test,0.174,2)
lines(z,y_hat,col="blue")

h=0.174
p=5
titre=paste("R�gression polynomiale � noyau de degr�",p,"et d'�cart-type",h)
plot(x_test,y_test,sub=paste(titre,h),main="R�gression de y sur x ")
y_hat<-regpolynoy(z,x_test,y_test,0.174,5)
lines(z,y_hat,col="blue")

h=1
p=2
titre=paste("R�gression polynomiale � noyau de degr�",p,"et d'�cart-type",h)
plot(x_test,y_test,sub=paste(titre,h),main="R�gression de y sur x ")
y_hat<-regpolynoy(z,x_test,y_test,1,2)
lines(z,y_hat,col="blue")


loonwp<-function(x,y,h,afhist){ #fait le leave one out 
  erreurquadr=c()
  pred_suc=list()
  ytilde_isuc=list()
  for (i in 1:length(y)){
    y_i=y[-i]
    x_i=x[-i]
    pred=regpolynoy(x[i],x_i,y_i,h,2)
    pred_suc=append(pred_suc,list(pred))
    erreurquadr_i=(y[i]-pred)^2
    erreurquadr[i]=erreurquadr_i  # calcule les erreurs successives de la pr?diction
  }
  erreurquadrmoy=mean(erreurquadr)
  if (!missing(afhist)){
    hist(erreurquadr,main=NULL,ylab=NULL,xlab=NULL) #on affiche ou non l'histogramme
  }
  return(list("pred_suc"=pred_suc,"erreurquadr"=erreurquadr,"erreurquadrmoy"=erreurquadrmoy))
}

loonwpch<-function(x,y,hach,afhist){ #permet d'avoir les erreur moyennes pour plusieurs h possibles
  errormoy<-list()
  for(h in hach){
    error<-loonwp(x,y,h,afhist)$erreurquadrmoy
    errormoy<-append(errormoy,list(c(error,h)))}
  
  return(errormoy)} 

h=seq(0.7,1,0.01)
loonwpch(x_test,y_test,h)
print("h=0.75, erreur= 0.002731755")

h=0.75
p=2
titre=paste("R�gression polynomiale � noyau de degr�",p,"et d'�cart-type",h)
plot(x_test,y_test,sub=paste(titre),main="R�gression de y sur x optimale")
y_hat<-regpolynoy(z,x_test,y_test,h,p)
lines(z,y_hat,col="blue")
abline(fit$coefficients[1],fit$coefficients[2],col="red")#erreur 0.00627 vs 0.002731755
legend(-1,0.6,legend=c("noyau","lin�aire"),col=c("blue","red"),lty=1:2,cex=0.8)
loonwpch(x_test,y_test,0.75,T)
title(main="Erreurs quadratiques de pr�diction ",sub="pour une r�gression polynomiale � noyau d'�cart-type 0.75 de degr� 2 ",
      ylab="Fr�quence",xlab="erreurs quadratiques")


dataset <- as.matrix(read.csv("Food_Supply_Quantity_kg_Data.csv",sep=',',header=T,row.names=1))
dim(dataset)
dataset
obesity=as.double(dataset[,24]) #pourcentage d'ob�sit� par pays
obesity=obesity[-c(110,53,148)]
meat=scalebiaised(as.double((dataset[,9])))#quantit� de viande
meat=meat[-c(110,53,148)]
z=seq(-2,3,0.1)

h=0.79
p=2
titre=paste("R�gression polynomiale � noyau de degr�",p,"et d'�cart-type",h)
plot(meat,obesity,main="R�gression de la pr�valence de l'ob�sit� sur la consommation de viande",xlab="Consommation de viande",ylab="Pourcentage d'ob�sit�",sub=titre)
y_hat<-regpolynoy(z,meat,obesity,0.79,2)
lines(z,y_hat,col="blue")
fit=lm(obesity~meat)

loofit(meat,obesity) #erreur 68.80 vs 61.9229
abline(fit$coefficients[1],fit$coefficients[2],col="red") #comparaison avec une r�gression lin�aire simple
legend(-1,40.6,legend=c("noyau","linÃ©aire"),col=c("blue","red"),lty=1:2,cex=0.8)


h=seq(0.1,3,0.01)
h=0.79
loonwpch(meat,obesity,h)
print("h=0.79 erreur= 61.9229")
loonwpch(meat,obesity,0.79,T)
title(main="Erreurs quadratiques de pr�diction ",sub="pour une r�gression polynomiale � noyau d'�cart-type 0.79 de degr� 2 ",
      ylab="Fr�quence",xlab="erreurs quadratiques")


cereal=scalebiaised(as.double((dataset[,5])))#quantit� de c�r�ales 5
dataset
cereal=cereal[-c(110,53,148)]
z=seq(-2,4,0.1)
h=5
p=2
titre=paste("R�gression polynomiale � noyau de degr�",p,"et d'�cart-type",h)
plot(cereal,obesity,main="R�gression de la pr�valence de l'ob�sit� sur la consommation de c�r�ales",xlab="Consommation de c�r�ales",ylab="Pourcentage d'ob�sit�",sub=titre)

y_hat<-regpolynoy(z,cereal,obesity,5,2)
lines(z,y_hat,col="blue")

h=seq(5,50,1)
h=5
loonwpch(cereal,obesity,h)
print("h=5 erreur= 67.00844")
fit=lm(obesity~cereal)
loofit(cereal,obesity) #erreur 66.10071 
abline(fit$coefficients[1],fit$coefficients[2],col="red") #comparaison avec une r�gression lin�aire simple
legend(-1,40,legend=c("noyau","lin�aire"),col=c("blue","red"),lty=1:2,cex=0.8)
loonwpch(cereal,obesity,5,T)
title(main="Erreurs quadratiques de pr�diction ",sub="pour une r�gression polynomiale � noyau d'�cart-type 5 de degr� 2 ",
      ylab="Fr�quence",xlab="erreurs quadratiques")


matecpuismulti<-function(X,x,d) { #la matice pour mutidimensionnel et polynomial
  z=cbind(rep(1,dim(X)[1])) #lignes = individus
  for (j in 1:d){
    r=c()
    for(i in 1:dim(X)[1]){
      r<-rbind(r,(X[i,]-x)**j)#  x c'est un vecteur 
    }
    z=cbind(z,r)
    
  }
  
  return(z)
}


regpolynoymulti<-function(Z,X,y,h,p){
  beta_hat=list()
  y_hat=c()
  pol=c()
  pol=matrix(rep(0,dim(Z)[1]*dim(Z)[2]),ncol=dim(Z)[1],nrow=dim(Z)[2]) #on pr?pare le polyn?me o? les lignes correspondent au variables
  for(k in 1:dim(Z)[1]){#pour chaque nouvelle ligne de Z
    H=matecpuismulti(X,Z[k,],p)# pour chaque vecteur ligne on calcule la matrice degr? p 
    beta<-solve(t(H)%*%diag(apply(X,1,function(l) kernelgaussian(l-Z[k,],h,dim(X)[2])))%*%H)%*%t(H)%*%diag(apply(X,1,function(l)
      kernelgaussian(l-Z[k,],h,dim(X)[2])))%*%y #on applique ? chaque vecteur ligne de X
    beta_hat=append(beta_hat,list(beta)) #on a r?solu le programme 
    b1=beta[1] #on prend le premier coefficient
    beta=beta[-1]
    
    y_hat[k]=b1
  }
  
  return(y_hat)
}

loonwpmulti<-function(X,y,h,afhist){ #fait le leave one out 
  erreurquadr=c()
  pred_suc=list()
  ytilde_isuc=list()
  for (i in 1:length(y)){
    y_i=y[-i]
    x_i=X[-i,]
    pred=regpolynoymulti(rbind(X[i,]),cbind(x_i),cbind(y_i),h,2)[[1]]
    pred_suc=append(pred_suc,list(pred))
    erreurquadr_i=(y[i]-pred)^2
    erreurquadr[i]=erreurquadr_i  # calcule les erreurs successives de la pr?diction
  }
  erreurquadrmoy=mean(erreurquadr)
  if (!missing(afhist)){
    hist(erreurquadr,main=NULL,ylab=NULL,xlab=NULL) #on affiche ou non l'histogramme
  }
  return(list("pred_suc"=pred_suc,"erreurquadr"=erreurquadr,"erreurquadrmoy"=erreurquadrmoy))
}

loonwpmultich<-function(x,y,hach,afhist){ #permet d'avoir les erreur moyennes pour plusieurs h possibles
  errormoy<-list()
  for(h in hach){
    error<-loonwpmulti(X,y,h,afhist)$erreurquadrmoy
    errormoy<-append(errormoy,list(c(error,h)))}
  
  return(errormoy)} 

library(readODS)
dataset <- data.matrix(read_ods("Datasim1.ods",sheet=2))
x_test1<-scalebiaised(dataset[,2])
x_test2<-scalebiaised(dataset[,3])
y_test<-dataset[,1]
y_hat=regpolynoymulti(cbind(seq(-1.5,1.5,0.1),rep(0,length(seq(-1.5,1.5,0.1)))),cbind(x_test1,x_test2),y_test,12,2)
#plot(seq(-1.5,1.5,0.1),y_hat,type='l',main="Estimation de la rÃ©gression polynomiale Ã  noyau sur un segment",
#     sub="Sur des donnÃ©es rÃ©elles h=12 er p=2") #estimation de y sur un segment
y_hat=regpolynoymulti(cbind(x_test1,x_test2),cbind(x_test1,x_test2),y_test,50,2)
plot(y_test,y_hat) #correspondance y yhat
h=seq(1,10,0.1)
h=10
loonwpmultich(cbind(x_test1,x_test2),y_test,h) #erreur 0.0028905 pour h=10
X=cbind(x_test1,x_test2)
y_i=y_test[-1]
x_i=X[-1,]
h=1
pred=regpolynoymulti(rbind(X[1,]),cbind(x_i),cbind(y_i),h,2)[[1]]


print("h=10 erreur =0.002899354")
y_hat=regpolynoymulti(cbind(x_test1,x_test2),cbind(x_test1,x_test2),y_test,10,2)
plot(y_test,y_hat,main="Correspondance entre y observ�s et estim�s",
     sub="dans une r�gression multidimensionnelle � noyau , h=10 et p=2")
loofit(cbind(x_test1,x_test2),y_test)#erreur 0.006454682 
abline(0,1,col='blue') 
legend(1,1.6,legend="bissectrice",col="blue",lty=1:2,cex=0.8)
loonwpmultich(cbind(x_test1,x_test2),y_test,10,T)
title(main="Erreurs quadratiques de pr�diction ",sub="pour une r�gression polynomiale � noyau multidimensionnel d'�cart-type 10 de degr� 2 ",
      ylab="Fr�quence",xlab="erreurs quadratiques")


dataset <- as.matrix(read.csv("Food_Supply_Quantity_kg_Data.csv",sep=',',header=T,row.names=1))
dim(dataset)
dataset
obesity=as.double(dataset[,24]) #pourcentage d'ob�sit� par pays
obesity=obesity[-c(110,53,148)]
meat=as.double((dataset[,9]))#quantit� de viande
meat=scalebiaised(meat[-c(110,53,148)])
cereal=as.double((dataset[,5]))#quantit� de c�r�ales 3
cereal=scalebiaised(cereal[-c(110,53,148)])


y_hat=regpolynoymulti(cbind(meat,cereal),cbind(meat,cereal),obesity,1.2,2)
plot(obesity,y_hat,main="Correspondance entre l'ob�sit� observ� et estim�",
     sub="dans une r�gression multidimensionnelle � noyau d'�cart-type 1.2 et de degr� 2",ylab="ob�sit�_estim�",xlab="pourcentage d'ob�sit�")
abline(0,1,col='blue')
legend(30,15,legend="bissectrice",col="blue",lty=1:2,cex=0.8)
loonwpmultich(cbind(meat,cereal),obesity,1.2,T)
title(main="Erreurs quadratiques de pr�diction ",sub="pour une r�gression polynomiale � noyau d'�cart-type 1.2 de degr� 2 ",
      ylab="Fr�quence",xlab="erreurs quadratiques")

fit=lm(obesity~cbind(meat,cereal))
plot(obesity,predict(fit),main="Correspondance entre y observ�s et estim�s"
     ,sub="dans une r�gression lin�aire",ylab="ob�sit�_estim�",xlab="pourcentage d'ob�sit�")
abline(0,1,col='blue')
legend(30,15,legend="bissectrice",col="blue",lty=1:2,cex=0.8)


h=seq(1,2,0.1)
h=1.2
loonwpmultich(cbind(meat,cereal),obesity,h)
X=cbind(meat,cereal)
y_i=obesity[-1]
x_i=X[-1,]
h=1
pred=regpolynoymulti(rbind(X[1,]),cbind(x_i),cbind(y_i),h,2)[[1]]
print('h=1.2 error=52.57091')

loofit(cbind(meat,cereal),obesity) #57.20668 


X=cbind(meat,cereal)
y_i=obesity[-1]
x_i=X[-1,]
h=1
pred=regpolynoymulti(rbind(X[1,]),cbind(x_i),cbind(y_i),h,2)[[1]]


# les splines
loosplin<-function(x,y,hach,afhist){
  errormoy<-list()
  for (h in hach){
    error=c()
    for(i in 1:length(x)){
      x_i=x[-i]
      y_i=y[-i]
      s=smooth.spline(x_i,y_i,cv=T,lambda=h)
      error[i]=(predict(s,x=x[i])$y-y[i])^2
    }
    errorquadr=mean(error)
    if (!missing(afhist)){
      hist(error,main=NULL,ylab=NULL,xlab=NULL) #on affiche ou non l'histogramme
    }
    
    errormoy<-append(errormoy,list(c(errorquadr,h)))
  }
  return(errormoy)
}


z=seq(-3,3,0.1)
plot(meat,obesity,main="R�gression de la pr�valence de l'ob�sit� sur la consommation de viande",sub="Comparaison spline et noyau avec p�nalit� de coeficient 0.00859",xlab="Consommation de viande",ylab="Pourcentage d'ob�sit�")
y_hat<-regpolynoy(z,meat,obesity,0.79,2)
s1<-smooth.spline(meat,obesity,cv=T)
lines(z,y_hat,col="blue")
lines(z,predict(s1,x=z)$y,col="red") 
legend(-1,40,legend=c("noyau","spline"),col=c("blue","red"),lty=1:2,cex=0.8)
s1 #lambda=0.00859
loosplin(meat,obesity,0.00859,T) #ereur 62.13156 vs 61,9229
title(main="Erreurs quadratiques de pr�diction ",sub="pour une r�gression spline p�nalis� de coefficient 0.00859 ",
      ylab="Fr�quence",xlab="erreurs quadratiques")

h=seq(1,2,0.01)
#loonwpch(meat,obesity,h)
print("h=0.79 erreur= 61.9229")

z=seq(-3,3,0.1)
plot(meat,obesity,main="R�gression de la pr�valence de l'ob�sit� sur la consommation de c�r�ales",sub="Comparaison spline et noyau avec p�nalit� de coeficient 3.891",xlab="Consommation de c�r�ales",ylab="Pourcentage d'ob�sit�")
y_hat<-regpolynoy(z,cereal,obesity,15,2)
lines(z,y_hat,col="blue")
s2<-smooth.spline(cereal,obesity,cv=T)
lines(z,predict(s2,x=z)$y,col="red") 
legend(-1,40,legend=c("noyau","spline"),col=c("blue","red"),lty=1:2,cex=0.8)
s2 #lambda= 3.891786
loosplin(cereal,obesity,3.891789,T) #ereur 66.142971 vs 67.00844
title(main="Erreurs quadratiques de pr�diction ",sub="pour une r�gression spline p�nalis� de coefficient 3.891789 ",
      ylab="Fr�quence",xlab="erreurs quadratiques")

print("h=5 erreur= 67.00844")

dataset <- data.matrix(read_ods("Datasim1.ods"))
x_test<-scalebiaised(dataset[,2])
y_test<-dataset[,1]
plot(x_test,y_test,sub="Comparaison spline et noyau avec p�nalit� de coefficient 0.005877111",main="R�gression de y sur x ")
y_hat<-regpolynoy(z,x_test,y_test,0.75,2)
lines(z,y_hat,col="blue")
s3<-smooth.spline(x_test,y_test,cv=T)
lines(z,predict(s3,x=z)$y,col="red") 
legend(-1,0.6,legend=c("noyau","spline"),col=c("blue","red"),lty=1:2,cex=0.8)
s3 #lambda=0.005877111
loosplin(x_test,y_test,0.005877111,T) #ereur 0.0027799500
title(main="Erreurs quadratiques de pr�diction ",sub="pour une r�gression spline p�nalis� de coefficient 0.005877111 ",
      ylab="Fr�quence",xlab="erreurs quadratiques")

print("h=0.75, erreur= 0.002731755")



