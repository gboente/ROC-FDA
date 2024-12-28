##############################
# Data Metabolic Syndrome    #
##############################

rm(list = ls())


library(fda)
library(scales)

library(lattice)
library('robustbase') # lmrob
library('MASS')       # rlm

source('Functions.R')
source('Functions_different_ROC.R')


cardio<- read.csv("BC_cardiotox_functional_variable.csv",header=T,sep=";",dec = ",")



fix(cardio)

attach(cardio)
dim(cardio)

names(cardio)

###################################################################
## Data from A cardiotoxicity dataset for breast cancer patients   #
## Beatriz Piñeiro-Lamas et al. (2024)                            #
## www.nature.com/scientificdata                                  #
###################################################################

oldpar <- par(mfrow = c(1,1))
n1 <- sum(CTRCD)
n0 <- dim(cardio)[1]-n1
grilla.tes <- seq(0,1,length=1001)


nombre="matplot.pdf"
pdf(nombre, bg='transparent')
par(mfrow=c(1,1))
par(mar=c(2,2,3,3))

matplot(grilla.tes,t(cardio[CTRCD==0,2:1002]),col="black",type = "l",xlab = "Cycle", ylab = " ",lty=1,lwd=2)
matplot(grilla.tes,t(cardio[CTRCD==1,2:1002]),col="magenta",type = "l",xlab = "Cycle", ylab = " ",lty=1,lwd=2,add=TRUE)


media0=apply(cardio[CTRCD==0,2:1002], 2, mean)
media1=apply(cardio[CTRCD==1,2:1002], 2, mean)

lines(grilla.tes, media0, col="green", lwd=3)
lines(grilla.tes, media1, col="gold",lwd=3)

legend(0.6,15,col=c("black","magenta","green","gold"),
       legend=c("CTRCD=0", "CTCRD=1", "Mean CTRCD=0","Mean CTRCD=1"), 
       text.col=c("black","magenta","green","gold"), bty="n",cex=0.8, lty=c(1,1))

dev.off()



############################################
# Definition of XD and XH                    #
############################################

XD <- as.matrix(cardio[CTRCD==1,2:1002]) # "CTRCD=1"
XH <- as.matrix(cardio[CTRCD==0,2:1002]) # "CTRCD=0"

nD <- dim(XD)[1]
nH <- dim(XH)[1]


matplot(grilla.tes,t(XH),col="black",type = "l",xlab = "Cycle", ylab = "Velocity (cm/s)",lty=1,lwd=2)
matplot(grilla.tes,t(XD),col="magenta",type = "l",xlab = "Cycle", ylab = "Velocity (cm/s)",lty=1,lwd=2,add=TRUE)

legend(0.8,14,col=c("black","magenta"),
       legend=c("CTRCD=0", "CTCRD=1"), 
       text.col=c("black","magenta"), bty="n",cex=0.8, lty=c(1,1))


#### Different colours

XD <- as.matrix(cardio[CTRCD==1,2:1002]) # "CTRCD=1"
XH <- as.matrix(cardio[CTRCD==0,2:1002]) # "CTRCD=0"

nD <- dim(XD)[1]
nH <- dim(XH)[1]


matplot(grilla.tes,t(XH),col="black",type = "l",xlab = "Cycle", ylab = "Velocity (cm/s)",lty=1,lwd=2)
matplot(grilla.tes,t(XD),col="mediumspringgreen",type = "l",xlab = "Cycle", ylab = "Velocity (cm/s)",lty=1,lwd=2,add=TRUE)

legend(0.8,14,col=c("black","mediumspringgreen"),
       legend=c("CTRCD=0", "CTCRD=1"), 
       text.col=c("black","mediumspringgreen"), bty="n",cex=0.8, lty=c(1,1))

####

nombre="matplot_bis.pdf"
pdf(nombre, bg='transparent')
par(mfrow=c(1,1))
par(mar=c(4,4.5,3,3))

matplot(grilla.tes,t(cardio[CTRCD==0,2:1002]),col="gray40",type = "l",xlab = "Cycle", ylab = "Velocity (cm/s)",lty=1,lwd=2)
matplot(grilla.tes,t(cardio[CTRCD==1,2:1002]),col="mediumspringgreen",type = "l",xlab = "Cycle", ylab = "Velocity (cm/s)",lty=1,lwd=2,add=TRUE)


media0=apply(cardio[CTRCD==0,2:1002], 2, mean)
media1=apply(cardio[CTRCD==1,2:1002], 2, mean)

lines(grilla.tes, media0, col="black", lwd=4)
lines(grilla.tes, media1, col="gold",lwd=4)

# legend(0.6,15,col=c("black","mediumspringgreen","magenta","gold"),
#        legend=c("CTRCD=0", "CTCRD=1", "Mean CTRCD=0","Mean CTRCD=1"), 
#        text.col=c("black","mediumspringgreen","magenta","gold"), bty="n",cex=0.8, lty=c(1,1))

dev.off()


nombre="matplot_bis_sinmedias.pdf"
pdf(nombre, bg='transparent')
par(mfrow=c(1,1))
par(mar=c(4,4.5,3,3))

matplot(grilla.tes,t(cardio[CTRCD==0,2:1002]),col="gray40",type = "l",xlab = "Cycle", ylab = "Velocity (cm/s)",lty=1,lwd=2)
matplot(grilla.tes,t(cardio[CTRCD==1,2:1002]),col="mediumspringgreen",type = "l",xlab = "Cycle", ylab = "Velocity (cm/s)",lty=1,lwd=2,add=TRUE)

dev.off()



############################################
# Functional Boxplots                      #
############################################
library(fda)

nombre="fbplot_CTRCD_0.pdf"
pdf(nombre, bg='transparent')

par(mfrow=c(1,1))
par(mar=c(4,4.5,3,3))

fbplot(t(XH),x=grilla.tes,xlim=c(0,1),ylim=c(-20,15),ylab="Velocity (cm/s)",xlab=" ")

dev.off()


nombre="fbplot_CTRCD_1.pdf"
pdf(nombre, bg='transparent')

par(mar=c(4,4.5,3,3))

fbplot(t(XD),x=grilla.tes,xlim=c(0,1),ylim=c(-20,15),ylab="Velocity (cm/s)",xlab=" ")

dev.off()

############################################################
#  Call the function that computes all the ROC's           #
############################################################

varianza.exp<-0.95 #explained percentage
ROC.all<- VARIAS.ROC(XD , XH, pe=seq(0,1,length=101), varianza.exp, grilla.tes,  base="PC")




ROC.all$auc.min
#   0.6770309

ROC.all$auc.max
#   0.4546563

ROC.all$auc.mediat
#   0.5328456

ROC.all$auc.dif.mediana
#   0.6867856

ROC.all$auc.dif.medias
#   0.6819082

ROC.all$AUC.PC.lineal
#   0.7033989

ROC.all$AUC.PHI.PC.lineal
#   0.6948398

ROC.all$AUC.PC.CUAD
#   0.8876696

ROC.all$kn
#  11

ROC.all$porcentaje
# 0.9541613


ROC.all.2<- VARIAS.ROC(XD , XH, pe=seq(0,1,length=101), varianza.exp, grilla.tes,  base="PC", pooled="NO")

 ROC.all.2$auc.min
#   0.6770309

ROC.all.2$auc.max
#   0.4546563

ROC.all.2$auc.mediat
#   0.5328456

ROC.all.2$auc.dif.mediana
#   0.6867856

ROC.all.2$auc.dif.medias
#   0.6819082

ROC.all.2$AUC.PC.lineal
#   0.7093431


ROC.all.2$AUC.PHI.PC.lineal
#     0.7072274

ROC.all.2$AUC.PC.CUAD
#   0.8876696

ROC.all.2$kn
#  11

ROC.all.2$porcentaje
#0.9541613

###############################################
# Plot  the ROC's                             #
###############################################

# > names(ROC.all)
# [1] "roc.min"           "auc.min"           "roc.max"           "auc.max"           "roc.mediat"       
# [6] "auc.mediat"        "roc.dif.mediana"   "auc.dif.mediana"   "roc.dif.medias"    "auc.dif.medias"   
# [11] "ROC.PC.lineal"     "AUC.PC.lineal"     "ROC.PHI.PC.lineal" "AUC.PHI.PC.lineal" "kn"               
# [16] "porcentaje"        "ROC.PC.CUAD"       "AUC.PC.CUAD"      

nombre="ROC-with-PC-linear-quadratic_min.pdf"
pdf(nombre, bg='transparent')

colores=c("black","magenta", "green", "gold","blue")
par(mfrow=c(1,1))
par(mar=c(4,4.5,3,3))

pe=seq(0,1,length=101)
plot(pe, ROC.all$ROC.PC.CUAD, type="l", col=colores[1],ylim=c(0,1),lwd=3,
     xlab="p",ylab=expression(hat(ROC)), 
     main=paste("k= ",ROC.all$kn,"(%=",varianza.exp,")"))
lines(pe, ROC.all$ROC.PC.lineal, type="l", col=colores[2], lwd=3) 

lines(pe, ROC.all$roc.dif.medias, type="l", col=colores[3], lwd=3)
lines(pe, ROC.all$roc.mediat, type="l", col=colores[4], lwd=3) 
lines(pe, ROC.all$roc.min, type="l", col=colores[5], lwd=3) 

lines(pe,pe,lty=2)


legend("bottomright",col=colores, lty=c(1,1,1,1,1), lwd=2,
       legend=c("QUADRATIC", "LINEAR", expression(mu[D]-mu[H]), 
                expression(paste(beta, "=1")),"MIN" ), text.col=colores, bty="n",cex=0.8)


dev.off()


nombre="ROC-with-PC-linear-quadratic_min_POOLED_NO.pdf"
pdf(nombre, bg='transparent')

colores=c("black","magenta", "green", "gold","blue")
par(mfrow=c(1,1))
par(mar=c(4,4.5,3,3))

pe=seq(0,1,length=101)
plot(pe, ROC.all.2$ROC.PC.CUAD, type="l", col=colores[1],ylim=c(0,1),lwd=3,
     xlab="p",ylab=expression(hat(ROC)), 
     main=paste("k= ",ROC.all$kn,"(%=",varianza.exp,")"))
lines(pe, ROC.all.2$ROC.PC.lineal, type="l", col=colores[2], lwd=3) 

lines(pe, ROC.all.2$roc.dif.medias, type="l", col=colores[3], lwd=3)
lines(pe, ROC.all.2$roc.mediat, type="l", col=colores[4], lwd=3) 
lines(pe, ROC.all.2$roc.min, type="l", col=colores[5], lwd=3) 

lines(pe,pe,lty=2)


legend("bottomright",col=colores, lty=c(1,1,1,1,1), lwd=2,
       legend=c("QUADRATIC", "LINEAR", expression(mu[D]-mu[H]), 
                expression(paste(beta, "=1")),"MIN" ), text.col=colores, bty="n",cex=0.8)


dev.off()



########################
# Another plot         #
########################

nombre="ROC-with-PC-linear-quadratic_min_dif.pdf"
pdf(nombre, bg='transparent')

colores=c("black","magenta", "green", "gold","blue")
par(mfrow=c(1,1))
par(mar=c(4,4.5,3,3))

pe=seq(0,1,length=101)
plot(pe, ROC.all$ROC.PC.CUAD, type="l", col=colores[1],ylim=c(0,1),lwd=3,
     xlab="p",ylab=expression(hat(ROC)), 
     main=paste("k= ",ROC.all$kn,"(%=",varianza.exp,")"))
lines(pe, ROC.all$ROC.PC.lineal, type="l", col=colores[2], lwd=3) 

lines(pe, ROC.all$roc.dif.medias, type="l", col=colores[3], lwd=3)
lines(pe, ROC.all$roc.dif.mediana, type="l", col=colores[4], lwd=3) 
lines(pe, ROC.all$roc.min, type="l", col=colores[5], lwd=3) 

lines(pe,pe,lty=2)


legend("bottomright",col=colores, lty=c(1,1,1,1,1), lwd=2,
       legend=c("QUADRATIC", "LINEAR", expression(mu[D]-mu[H]), expression(med[D]-med[H])
                ,"MIN" ), text.col=colores, bty="n",cex=0.8)


dev.off()

########################
# Another plot: only 4 #
########################

nombre="ROC-with-PC-linear-quadratic_min_difmean.pdf"
pdf(nombre, bg='transparent')

colores=c("black","magenta", "green","blue")
par(mfrow=c(1,1))
par(mar=c(4,4.5,3,3))

pe=seq(0,1,length=101)
plot(pe, ROC.all$ROC.PC.CUAD, type="l", col=colores[1],ylim=c(0,1),lwd=3,
     xlab="p",ylab=expression(hat(ROC))
#     ,main=paste("k= ",ROC.all$kn,"(%=",varianza.exp,")")
     )
lines(pe, ROC.all$ROC.PC.lineal, type="l", col=colores[2], lwd=3) 

lines(pe, ROC.all$roc.dif.medias, type="l", col=colores[3], lwd=3)
lines(pe, ROC.all$roc.min, type="l", col=colores[4], lwd=3) 

lines(pe,pe,lty=2)


# legend("bottomright",col=colores, lty=c(1,1,1,1), lwd=2,
#        legend=c("QUADRATIC", "LINEAR", expression(mu[D]-mu[H]), "MIN" ), text.col=colores, bty="n",cex=0.8)


dev.off()

nombre="ROC-with-PC-linear-quadratic_min_difmean_POOLED_NO.pdf"
pdf(nombre, bg='transparent')

colores=c("black","magenta", "green","blue")
par(mfrow=c(1,1))
par(mar=c(4,4.5,3,3))

pe=seq(0,1,length=101)
plot(pe, ROC.all.2$ROC.PC.CUAD, type="l", col=colores[1],ylim=c(0,1),lwd=3,
     xlab="p",ylab=expression(hat(ROC))
#  , main=paste("k= ",ROC.all$kn,"(%=",varianza.exp,")")
	)
lines(pe, ROC.all.2$ROC.PC.lineal, type="l", col=colores[2], lwd=3) 

lines(pe, ROC.all.2$roc.dif.medias, type="l", col=colores[3], lwd=3) 
lines(pe, ROC.all.2$roc.min, type="l", col=colores[4], lwd=3) 

lines(pe,pe,lty=2)


# legend("bottomright",col=colores, lty=c(1,1,1,1), lwd=2,
#        legend=c("QUADRATIC", "LINEAR", expression(mu[D]-mu[H]), "MIN" ), text.col=colores, bty="n",cex=0.8)

dev.off()
