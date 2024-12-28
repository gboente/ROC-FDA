source('Functions.R')


VARIAS.ROC<- function(XD , XH, pe=seq(0,1,length=101), varianza.exp, grilla.tes,  base="PC", pooled="SI"){

	nD <- dim(XD)[1]
	nH <- dim(XH)[1]

	if(dim(XD)[2]!=dim(XH)[2]){
		print('ERROR: SAMPLES RECORDED IN DIFFERENT GRID POINTS')
	}

	if(dim(XD)[2]!=length(grilla.tes)){
		print('ERROR: THE GRID POINTS DO NOT MATCH THE SAMPLE RECORDS')
	}

	minimo<-min(min(XD),min(XH))
	maximo<-max(max(XD),max(XH))



####################################################### 
# Empirical ROC curve and AUC based on the maximum    #  
# over t of the curve of each individual              #
#######################################################

	yH.max <- apply(XH,1,max)  # "HEALTHY"
	yD.max <- apply(XD,1,max)  # "DISEASED"


	roc.max<-ROC.clasica.dim1(yD.max,yH.max)
########################################################
# The grid is p=seq(0,1,length=101)                    #
# hence roc.max[1]==0                                  #
########################################################

	roc.max[1]= 0
	auc.max<-AUC(yD.max,yH.max)

####################################################### 
# Empirical ROC curve and AUC based on the minimum    #  
# over t of the curve of each individual              #
#######################################################

	yH.min <- apply(XH,1,min)  # "HEALTHY"


	yD.min <- apply(XD,1,min)  # "DISEASED"


	roc.min<-ROC.clasica.dim1(yD.min,yH.min)

########################################################
# The grid is p=seq(0,1,length=101)                    #
# hence roc.min[1]==0                                  #
########################################################

	roc.min[1]= 0

	auc.min<-AUC(yD.min,yH.min)




####################################################### 
# Empirical ROC curve and AUC based on the mean value #
# over t of the curve of each individual              #
# i.e., taking as the direction beta(t)=1             #
#######################################################

	yH <- rowMeans(XH)  # "HEALTHY"
	yD <- rowMeans(XD)  # "DISEASED"


	roc.mediat<-ROC.clasica.dim1(yD,yH)

	auc.mediat<- AUC(yD,yH)

#########################################################
# The grid is p=seq(0,1,length=101)                     #
# hence roc.mediat[1]==0                                #
#########################################################

	roc.mediat[1]= 0


##################################################
# Based on the difference of medians             #
##################################################

	medianaD<-l1.median(XD)
	medianaH<-l1.median(XH)
	direccion.mediana<- medianaD-medianaH


	mediaD=apply(XD, 2, mean)
	mediaH=apply(XH, 2, mean)

	yD.max.mediana <- as.matrix(XD)%*%direccion.mediana
	yH.max.mediana <- as.matrix(XH)%*%direccion.mediana

	roc.dif.mediana <- ROC.clasica.dim1(yD.max.mediana,yH.max.mediana)
	auc.dif.mediana<- AUC(yD.max.mediana,yH.max.mediana)

#########################################################
# The grid is p=seq(0,1,length=101)                     #
# hence roc.dif.mediana[1]==0                           #
#########################################################

	roc.dif.mediana[1]= 0

#########################################################
# Based on the difference of means                      #
#########################################################

	direccion.AUC <- mediaD-mediaH

	yD.max.AUC <- XD%*%direccion.AUC
	yH.max.AUC <- XH%*%direccion.AUC

	roc.dif.medias <- ROC.clasica.dim1(yD.max.AUC,yH.max.AUC)


	auc.dif.medias<- AUC(yD.max.AUC,yH.max.AUC)

#########################################################
# The grid is p=seq(0,1,length=101)                     #
# hence roc.dif.medias[1]==0                            #
######################################################### 

	roc.dif.medias[1]= 0


#########################################################
# Using  % of the explained variance: varianza.exp       #
#########################################################

 	RC.PC<- ROC.CL.funcional.alfa(grilla.tes, XD, XH,pe=pe,alfa=NULL, varexp=varianza.exp, base=base, pooled=pooled)

	kn.PC.lineal <- RC.PC$kn

	direccion.PC.lineal <- RC.PC$beta
	porcentaje.lineal <-RC.PC$porcentaje

	ROC.PC.lineal.proyectada  <- RC.PC$ROC.emp
	ROC.PHI.PC.lineal.proyectada <- RC.PC$ROC.phi

	AUC.PC.lineal.proyectada  <- RC.PC$AUC
	AUC.PHI.PC.lineal.proyectada <- mean(ROC.PHI.PC.lineal.proyectada )


######################################################## 
# Using PC and quadratic rule                          #
########################################################

	pe=seq(0,1,length=101)

	RC.PC.CUAD<- ROC.CUAD.funcional(grilla.tes, XD, XH,pe=pe,varexp=varianza.exp, base=base)

	kn.PC.CUAD<- RC.PC.CUAD$kn


	porcentaje.CUAD <-RC.PC.CUAD$porcentaje

	ROC.PC.CUAD.proyectada  <- RC.PC.CUAD$ROC.emp 

	AUC.PC.CUAD.proyectada <- RC.PC.CUAD$AUC 




return(list(roc.min=roc.min,auc.min=auc.min,roc.max=roc.max,auc.max=auc.max, 
	roc.mediat=roc.mediat, auc.mediat=auc.mediat,	
	roc.dif.mediana=roc.dif.mediana,auc.dif.mediana=auc.dif.mediana,
	roc.dif.medias=roc.dif.medias,auc.dif.medias=auc.dif.medias,
	ROC.PC.lineal=ROC.PC.lineal.proyectada, AUC.PC.lineal=AUC.PC.lineal.proyectada,
	ROC.PHI.PC.lineal=ROC.PHI.PC.lineal.proyectada, AUC.PHI.PC.lineal=AUC.PHI.PC.lineal.proyectada,
	kn=kn.PC.lineal, porcentaje=porcentaje.CUAD,
	ROC.PC.CUAD=ROC.PC.CUAD.proyectada, AUC.PC.CUAD=AUC.PC.CUAD.proyectada))

}
