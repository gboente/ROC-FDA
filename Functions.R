library(rrcov)
library(mvtnorm)
library(rrcov)
require('fda')        # splines


########################
# Classical ROC        #
########################

########################
# DIMENSION 1          #
########################

ROC.clasica.dim1<- function(xD,xH,p=seq(0,1,length=101)) {
  	xD.emp <- ecdf(xD)
	roc.cl = 1- xD.emp(quantile(xH,1-p))
	return(roc.cl)
}


# Empirical ROC curve

# Empirical distribution function
emp.df=function(n,x,point){sum(x<=point)/n}

# Empirical quantile function
emp.quantile=function(n,x,p){ifelse(n*p==floor(n*p),sort(x)[floor(n*p)],sort(x)[floor(n*p)+1])}


ROC.emp <- function(yD,yH,p=seq(0,1,length=101)){
  nD <-length(yD)
  nH <-length(yH)
  lp <- length(p)
  ROCp <- rep(0,lp)
  for (i in 1:lp){
    ROCp[i] <- 1-emp.df(nD,yD,emp.quantile(nH,yH,1-p[i]))
  }
  plot(p,ROCp,xlim=c(0,1),ylim=c(0,1),xlab="p",ylab="ROC(p)",type="s")
  abline(0,1,lty=3)
  AUC <- mean(outer(yD,yH,"-")>0) 
  return(list(p,ROCp,AUC))
}

AUC <- function(yD,yH){ mean(outer(yD,yH,"-")>0) }

ROC.CL.alfa<-function(yD, yH,pe,alfa=NULL){

	nD=dim(yD)[1]
	nH=dim(yH)[2]

################################################
# Classical estimators of location and scale   #
################################################

	mediaD  = colMeans(yD)
 	covaD   = cov(yD)

 	mediaH   =colMeans(yH)
 	covaH    = cov(yH)

	Sigma = (nD*covaD+nH*covaH)/(nD+nH)


 	if(is.null(alfa)){
		alfa= solve(Sigma)%*%(mediaD-mediaH)
	}


################################
# Binormal Case                #
################################

	norma= sqrt(mahalanobis(mediaD, mediaH, Sigma))


	zetap=qnorm(1-pe)

	ROC.phi= 1-pnorm(zetap-norma)


###################################################################
# Classical empirical function, projected over the direction alfa #
###################################################################

	xD = yD %*%alfa 
	xH = yH %*%alfa 

	largo=length(pe)

      ROC.emp=rep(NA,length=largo)
  	for ( i in 1:largo){
            punto = pe[i]
  		roc.cl = ROC.clasica.dim1(xD=xD,xH=xH,p=punto) 
    		ROC.emp[i]<- roc.cl
	  }	

	return(list(ROC.phi=ROC.phi, ROC.emp=ROC.emp, alfa=alfa))
 }

########################
# Optimal Beta AUC     #
########################
funcionese <- function(u){ 
		ese=1/(1+exp(-u))
		return(ese)
	}

funcionele <- function(u){
	ele= funcionese(u)*(1-funcionese(u))
	return(ele)
}


beta.auc.cl <- function(yD, yH,alfaini, ache){

	nD=dim(yD)[1]
	nH=dim(yH)[1]

	d=length(alfaini)
      d1 = d-1
	theta.ini=alfaini[2:d]/alfaini[1]

	yD.star= yD[, 2:d]
	yH.star= yH[, 2:d]

	if(d1==1){
		yD.star=as.matrix(yD.star,nrow=nD,ncol=1)
		yH.star=as.matrix(yH.star,nrow=nH,ncol=1)
	}

	z= rep(NA,nD*nH)
	Z.star <- matrix(rep(NA,d1*nD*nH),nrow=nD*nH,ncol=d1)

	indice=0
	for(i in 1:nD){ 
		for(j in 1:nH){
			indice=indice+1
			z[indice]= yD[i,1]-yH[j,1]
			Z.star[indice,] = yD.star[i,]- yH.star[j,]

		}
	}


	funcionobj <- function(theta){
		thetavec=theta
		arg=(z+ Z.star%*%thetavec)/ache
		objetivo=mean(funcionese(arg))
		return(objetivo)
	}

	derfuncionobj <- function(theta){
		thetavec=theta
		arg=(z+ Z.star%*%thetavec)/ache
		derobjetivo=mean(as.vector(funcionele(arg))*Z.star)
		return(derobjetivo)
	}	

	#theta.opt=optimize(funcionobj,interval=c(-7,-5),maximum=TRUE)$max  

	if(d1==1){
		maximizacion= optim(theta.ini, fn=funcionobj, gr=derfuncionobj,method="Brent",lower=-20,upper=20,control = list(fnscale=-1 ))
		theta.opt=maximizacion$par 
		}

	if(d1!=1){
		theta.opt=optim(theta.ini, fn=funcionobj, gr=derfuncionobj,control = list(fnscale=-1 ))$par 
		}


	return(theta.opt)
}


###############################################################
# The grid where the trajectories are observed                #
# needs to be the same                                        #
###############################################################
###############################################################
# Compute the norm of a data function (as list) with a        #
#  discretization of size h                                   #
###############################################################
L2.norm.h <- function(dato,h)
{
  return(sqrt(sum(dato*dato)*h))
}

#################################################################################
# Given a data collection (a matrix with data by row), it returns the  mean     #
#################################################################################

media.puntual <- function(data)
{
   apply(data,FUN=mean,MARGIN=2)
}

is.odd <- function(x) x %% 2 != 0


#########################################################################
# Using basis                                                           #
#########################################################################
# kn     : number of splines for beta, if base="splines"
#	 : number of the elements of the fourier basis if basis="fourier", in which case it must be odd
#	 : number of common principal components if basis="PC" 
# 	 : number of element of the true basis of principal direction if basis="verdad"
# varexp    : Percentage of variance explained, if varexp=1 the value of kn is taken as the number of PC
# pooled : Use the pooled covariance to compute the direction (option "SI")
#        : Use the average  between covariance to compute the direction (option "NO")

ROC.CL.funcional.alfa<-function(grilla.tes, XD, XH,pe=seq(0,1,length=100),alfa=NULL, kn=10,varexp=1, norder=4, base="splines", base.verdad, pooled="SI"){

	nD=dim(XD)[1]
	nH=dim(XH)[1]

	tt=length(grilla.tes)

	dt <- 1 / (tt - 1)
	porcentaje=NA

	if(base=="splines"){
#####################################################
# Generation of the splines basis                   #
#####################################################

	 	nodos.spl   <- seq(min(grilla.tes), max(grilla.tes), length = kn- norder + 2)
       		base.spl   <- create.bspline.basis(rangeval = c(min(grilla.tes), max(grilla.tes)), norder = norder, breaks = nodos.spl)
       		spl.en.grilla <- getbasismatrix(grilla.tes, base.spl)

       		base_proyec <-  spl.en.grilla
	}

	if(base=="PC"){

		mediaD.func <- media.puntual(XD)
		mediaH.func <- media.puntual(XH)

		XD.center <-  t(t(XD) - mediaD.func)
		XH.center <-  t(t(XH) - mediaH.func)

		xc <- rbind(XD.center, XH.center)
		ntot<- nD+nH
		cov.mat <- crossprod(xc) / ntot

		cov.dec <- eigen(cov.mat, symmetric = TRUE)

		autovalores <- cov.dec$values

		if(varexp!=1){
			porcen_acum <- cumsum(autovalores)/sum(autovalores)
 
#########################################
# Those with >= % fixed variance        #
#########################################

  			cuales <- 1*(varexp <= porcen_acum)
  			names(cuales) <- 1: length(autovalores)
  			distintoscero <- as.numeric(names(cuales[cuales!=0]))

  			kn= min(distintoscero)
		}


		porcentaje <- sum(autovalores[1:kn])/sum(autovalores) 

		print(porcentaje)

    		pc.CL  <- cov.dec$vectors[, 1:kn]
    		normas   <- sqrt(colSums(pc.CL^2) * dt)  #### To get norm L2 1
    		base_proyec <- pc.CL %*% diag(1 / normas)

	}

	if(base=="verdad"){

		pc.CL  <- base.verdad[, 1:kn]
    		normas   <- sqrt(colSums(pc.CL^2) * dt)  #### To get norm L2 1
    		base_proyec <- pc.CL %*% diag(1 / normas)

	}

	if(base=="fourier"){

		pepe<- is.odd(kn)
		if(pepe==TRUE){
			base.fourier   <- create.fourier.basis(rangeval=c(min(grilla.tes), max(grilla.tes)), nbasis=kn, period=diff(rangeval))

			base.fourier.grilla <- getbasismatrix(grilla.tes, base.fourier)

			normas   <- sqrt(colSums(base.fourier.grilla^2) * dt)  #### To get norm L2 1

			base_proyec <-   base.fourier.grilla %*% diag(1 / normas)
		}
	}

        ## Estimated   coefficients in the basis   (by row)

      yD <- XD %*% base_proyec * dt
	yH <- XH %*% base_proyec * dt


################################################
# Classical estimators of locations and scale  #
################################################

	mediaD  = colMeans(yD)
 	covaD   = cov(yD)

 	mediaH   =colMeans(yH)
 	covaH    = cov(yH)

	if(pooled=="SI"){
		Sigma = (nD*covaD+nH*covaH)/(nD+nH)
	}
	if(pooled!="SI"){
		Sigma =  (covaD+ covaH)/2
	}

 	if(is.null(alfa)){
		alfa= solve(Sigma)%*%(mediaD-mediaH)
	}

#####################################################
# Reconstruction of the functional Beta             #
#####################################################

	beta <-   base_proyec %*% alfa

	posicionD <-  base_proyec %*% mediaD
	posicionH <-  base_proyec %*% mediaH

################################
# Binormal Case                #
################################

	norma= sqrt(mahalanobis(mediaD, mediaH, Sigma))


	zetap=qnorm(1-pe)

	ROC.phi= 1-pnorm(zetap-norma)


######################################################################
# Classical empirical function, projected over the direction alfa    #
######################################################################

	xD = yD %*%alfa 
	xH = yH %*%alfa 

	largo=length(pe)

      ROC.emp=rep(NA,length=largo)
  	for ( i in 1:largo){
            punto = pe[i]

  		roc.cl = ROC.clasica.dim1(xD=xD,xH=xH,p=punto) 
    		ROC.emp[i]<- roc.cl

		if(punto==0){ROC.emp[i]=0}
		if(punto==1){ROC.emp[i]=1}
	  }	


	AREA.PROY<- AUC(xD,xH) 

	return(list(ROC.phi=ROC.phi, ROC.emp=ROC.emp, alfa=alfa, beta=beta,kn=kn, AUC=AREA.PROY,mediaD=posicionD,mediaH=posicionH, porcentaje=porcentaje, xD=xD, xH=xH))
 }


ROC.CUAD.funcional<-function(grilla.tes, XD, XH,pe=seq(0,1,length=100), kn=10,varexp=1, base="PC", base.verdad) {

		nD=dim(XD)[1]
		nH=dim(XH)[1]

		tt=length(grilla.tes)

		dt <- 1 / (tt - 1)

		if(base=="PC"){

			mediaD.func <- media.puntual(XD)
			mediaH.func <- media.puntual(XH)

			XD.center <-  t(t(XD) - mediaD.func)
			XH.center <-  t(t(XH) - mediaH.func)

			xc <- rbind(XD.center, XH.center)
			ntot<- nD+nH
			cov.mat <- crossprod(xc) / ntot

			cov.dec <- eigen(cov.mat, symmetric = TRUE)

			autovalores <- cov.dec$values

			if(varexp!=1){
				porcen_acum <- cumsum(autovalores)/sum(autovalores)

 #####################################
 # Those with >= % fixed variance    #
 #####################################

  				cuales <- 1*(varexp <= porcen_acum)
  				names(cuales) <- 1: length(autovalores)
  				distintoscero <- as.numeric(names(cuales[cuales!=0]))

  				kn= min(distintoscero)
			}

			porcentaje <- sum(autovalores[1:kn])/sum(autovalores) 

			print(porcentaje)
    			pc.CL  <- cov.dec$vectors[, 1:kn]
    			normas   <- sqrt(colSums(pc.CL^2) * dt)  #### To get norm L2 1
    			base_proyec <- pc.CL %*% diag(1 / normas)

	}

	if(base=="verdad"){

		pc.CL  <- base.verdad[, 1:kn]
    		normas   <- sqrt(colSums(pc.CL^2) * dt)  #### To get norm L2 1
    		base_proyec <- pc.CL %*% diag(1 / normas)

	}


        ## Estimated coefficients in the basis (by row)

        yD <- XD %*% base_proyec * dt
	yH <- XH %*% base_proyec * dt



################################################
# Classical estimators of location and scale   #
################################################

	mediaD  = colMeans(yD)
 	covaD   = cov(yD)

 	mediaH   =colMeans(yH)
 	covaH    = cov(yH)


	#cte <- log(det(covaH)/det(covaD)) - (mahalanobis(mediaD, rep(0,kn), covaD) - mahalanobis(mediaH, rep(0,kn), covaH))
	cuad<- (solve(covaD)-solve(covaH)) 
	alfa<-  2*(solve(covaD)%*%mediaD-solve(covaH)%*%mediaH)


	posicionD <-  base_proyec %*% mediaD
	posicionH <-  base_proyec %*% mediaH
	
######################################################################
# Classical empirical function, projected over the direction alfa    #
######################################################################

	xD =   - diag(yD %*%cuad%*% t(yD)) +yD%*%alfa
	xH =   - diag(yH %*%cuad%*% t(yH)) +yH%*%alfa

	largo=length(pe)

      ROC.emp=rep(NA,length=largo)
  	for ( i in 1:largo){
            punto = pe[i]
  		roc.cl = ROC.clasica.dim1(xD=xD,xH=xH,p=punto) 
    		ROC.emp[i]<- roc.cl


		if(punto==0){ROC.emp[i]=0}
		if(punto==1){ROC.emp[i]=1}

	  }	


	AREA.PROY<- AUC(xD,xH) 

	return(list(ROC.emp=ROC.emp, alfa=alfa, matriz=cuad, kn=kn, AUC=AREA.PROY, mediaD=posicionD, mediaH=posicionH, porcentaje=porcentaje, xD=xD, xH=xH))
 }


##############
# MEDIAN L1  #
###################################################################################
# Given a data collection (a matrix with data by row), it returns the median  L1  #
###################################################################################

l1.median <- function(data)
{
  if(ncol(data)==1){
    return(median(data[,1]))
  }
  return(pcaPP::l1median(data))

}

