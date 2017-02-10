# data(sysdata, envir=environment())

.demog.nat <-
function(yrborn, sex, demographics.df=NULL){

  i<-age.m<-age.d<-deg.1.demog<-NULL  #initialize
  if(is.null(demographics.df)){
    print("No demographics given.  Using USDemographics.df")
    demographics.df=USDemographics.df
  }

  #inyr<-c(1800, 1900, 1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000) #entry year
  #out.year<-demographics.df$in.year+c(100, rep(10, each=10), 200)  #exit year extend to 2050.
  #f.age.mar<-c(22.0, 21.09, 21.6, 21.2, 21.3, 21.5, 20.3, 20.8, 22.0, 23.9, 25.1, 26.1)  #average marriage age
  #m.age.mar<-c(26.1, 26.1, 25.1, 24.6, 24.3, 24.3, 22.8, 23.2, 24.7, 26.1, 26.8, 28.2)
  #national data below 1800 assumed to be the same as 1890
  # f.age.death<-c(46.3, 48.3, 51.8, 54.6, 61.6, 65.2, 71.1, 74.7, 77.4, 78.8, 79.7, 81.1) #average life expectancy including infant mortality(what would be better is average life span excluding infant mortality, we know all people lived to reproductive age)
  #m.age.death<-c(44.3, 46.3, 48.4, 53.6, 58.1, 60.8, 65.6, 67.1, 70.0, 71.8, 74.3, 76.2) # average life expectancy, including infant mortality
  #f.age.death<-c(60.2, 62.03, 63.77, 64.88, 66.46, 68.52, 71.38, 74.56, 76.29, 77.244, 79.44, 80.3) # based on http://www.infoplease.com/ipa/A0005140.html life expectancy at 20
  #m.age.death<-c(60.1, 60.66, 62.19, 62.71, 65.6, 66.02, 67.76, 69.52, 70.25, 70.22, 72.45, 74) # based on http://www.infoplease.com/ipa/A0005140.html life expectancy at 20

  #aveh<-c(4.93, 4.76, 4.54, 4.34, 4.11, 3.77, 3.52, 3.14, 2.76, 2.63, 2.62, 2.61) #average household
  #offsp<-c((2.66*0.67), (2.62*.77), (2.50*.8), (2.09*.83), (1.97*.92), (2.87*.94), (3.29*.962), (2.07*.97), (1.67*.977), (1.87*.985), (1.92*.89), (1.90*.991)) #average children/woman * (1-child mortality)
  #offsp<-c((4.7*0.67), (3.8*.77), (3.6*.8), (3.3*.83), (2.4*.92), (2.2*.94), (3*.962), (3.7*.97), (2.5*.977), (1.8*.985), (2.1*.89), (2.1*.991)) #average children/woman from gapminder * (1-child mortality)
   ####
  if(sex==1){  #age by sex
    age.m<-demographics.df$female.age.marriage
    age.d<-demographics.df$female.age.death
  }else{
    age.m<-demographics.df$male.age.marriage
    age.d<-demographics.df$male.age.death
  }

  age.mar<-offs<-NULL   #marriage age, average children
  for (i in 1:length(demographics.df$in.year)){     ## generate marriage year from demographics
    if (round(as.numeric(yrborn) + age.m[i]) >= demographics.df$in.year[i]-9 & round(as.numeric(yrborn) + age.m[i]) < demographics.df$out.year[i]+9){
      #not monotonic so add fudge factor
      age.mar<-age.m[i] # for marriage age and average children per woman
      offs<-demographics.df$offspring[i]
     }
  }
  age.death<-NULL  #death age
  for (i in 1:length(demographics.df$in.year)){   ## generage death age
    if (yrborn >= demographics.df$in.year[i] & yrborn  < demographics.df$out.year[i]){
      age.death<-age.d[i] # for marriage age and average children per woman
    }
  }
  ##create demographics using tables from population demographic and distributions roughly skewed roughly equal to reality
  deg.1.demog<-list("age.mar"=ifelse(is.null(age.mar)==FALSE,round(rsnorm(1, age.mar, sd = 2.5, xi = 1.5),3),NA ), "age.death"=ifelse(is.null(age.mar)==FALSE,round(rsnorm(1, age.death, sd = (age.death/6), xi = .8),3),NA ), "offs"=ifelse(is.null(offs)==FALSE,rpois(1, lambda=abs(offs)),NA))
  return(deg.1.demog)
}

.demog <-
function(offnum, aveage, sdage){
  n<-NULL
  if(sum(n)==0){
  n<-rpois(1, lambda=abs(offnum))
  }#generate num of offsprings
  print(n)
  age<-round(rnorm(n, aveage, sdage),3)  #generage age
  female<-sample(0:1, n, replace=T) #randomly assigne geneder
  id<-c(1:n)#within family id
  ped<-list("id"=id, "age"=age, "female"=female)
return(ped)
}

.demog.2 <-
function(offnum, y.birth, demographics.df=NULL){
  age<-age.temp<-dead<-y.born<-NULL

  if(is.null(offnum)){
    offnum <-1
  } ## make sure there are generations connecting proband's generation and ancestor.
  if(is.na(offnum)){
    offnum <-1
  }
  if(offnum <= 0){
    offnum <-1
  }
  #  repeat{
  #    offn<-rpois(1, lambda=abs(offnum))
  #    offn<-ifelse(is.na(offn)==TRUE, 0, offn)
  #    if(offn>0) break()
  #  }
  offn <- offnum

  #initializing the arrays
  age=rep(0,offn)
  age.temp=age
  dead=age
  y.born=age

  female<-sample(0:1, offn, replace=T) #randomly assigne geneder
  id<-c(1:offn)#within family id

  for (i in 1:offn){
    y.born[i] <-ceiling(round(rnorm(1, y.birth, 5),3)) #initialize from y.birth of generation
    age.temp[i]<-round(2010-y.born[i])
    age[i]<-.demog.nat(y.born[i],female[i],demographics.df)$age.death
    #age[i]<-.demog.nat.china(y.born[i],female[i])$age.death
    dead[i]<-1
    if (age.temp[i]< round(rnorm(1, 81.1, 5),3)){ #took the older age
      age[i]<-age.temp[i]
      dead[i]<-0
    }
  }
  ped<-list("id"=id, "age"=age, "female"=female, "y.born"=y.born, "dead"=dead)
  #print(ped)
return(ped)
}

.genotype <-
function(mom, dad, ped){
  fromdad<-frommom<-NULL

  if(all(is.numeric(ped$id))==FALSE&all(is.numeric(ped$age))==FALSE&all(is.numeric(ped$female))==FALSE&all(is.numeric(ped$y.born))==FALSE){
    ped$geno<-NA
  }else{
  # allele from dad
 if(dad==1) {
   fromdad <- sample(0:1, length(ped$id), replace=TRUE)
   while(sum(fromdad)==0){
     fromdad <- sample(0:1, length(ped$id), replace=TRUE)
   }
  } else if(dad==0){
   fromdad <- sample(0, length(ped$id), replace=TRUE)
  }
  # allele from mom
 if(mom==1) {
     frommom <- sample(0:1, length(ped$id), replace=TRUE)
     while(sum(frommom)==0){
      frommom <- sample(0:1, length(ped$id), replace=TRUE)
     }
  }else if(mom==0){
    frommom <- sample(0, length(ped$id), replace=TRUE)
  }
  # return kids' genotypes
  ped$geno<-fromdad + frommom
  }
return(ped)
}

.grow.p <-
function(prev.dat,demographics.df=NULL){
  deg.num <- max(prev.dat$degree)
  next.deg <- deg.num+1
  tcur.deg <- subset(prev.dat, prev.dat$degree==deg.num)
  i<- toffs<- deg<-fem<-mal<-dad<-mom<-momid<-dadid<-geno<-female<-id<-dead <-y.born<-age<- age.temp<-ids<-temp<-next.dat<- NULL
  tids <- tcur.deg$id
  kids <- 0

  #initializing arrays
  toffs=rep(0,nrow(tcur.deg))

  while(kids <1){   ## while loop to make sure there is at least one person with offspring
    for (i in 1:nrow(tcur.deg)){
      toffs[i] <- .demog.nat(tcur.deg[i,]$y.born, 1,demographics.df)$offs  #check for offspring #use offspring regardless of gender for now#, if none do not proceed
      #toffs[i] <- .demog.nat.china(tcur.deg[i,]$y.born, 1)$offs  #check for offspring #use offspring regardless of gender for now#, if none do not proceed
    }
    #print(toffs)
    kids <- sum(toffs,na.rm = TRUE)
  }
  cur.deg.1 <- cbind(tcur.deg, toffs)
  cur.deg <- subset(cur.deg.1, cur.deg.1$toffs > 0) ## do not generate spouse if no offspring
  toffs <- cur.deg$toffs
  #print(toffs)
  cur.deg$toffs <- NULL
  ids <- cur.deg$id  ## get id for individuals with offspring

  #initializing more arrays
  geno=female=y.born=age=id=age.temp=dead=momid=dadid=fem=mal=dad=mom=rep(0,nrow(cur.deg))

  for (i in 1:nrow(cur.deg)){
    ###generate spouse
    geno[i] <- 0 # these should all be 0 no carriers marry in # abs(prev.dat[prev.dat$id==car[i],]$geno-1)# 0 since choose the carrier
    female[i]<-abs(cur.deg[cur.deg$id==ids[i],]$female-1) # gender opposite
    if(female[i] == 0){
      y.born[i]<-round(rnorm(1, (cur.deg[cur.deg$id==ids[i],]$y.born)-2, 3),3)  # choose age, men tend to be older than the women they marry
     age[i] <- round(rnorm(1, (cur.deg[cur.deg$id==ids[i],]$age)-2, 3),3)  #choose age at death women usually outlive their spouses
    }
    if(female[i] == 1){
      y.born[i]<-round(rnorm(1, (cur.deg[cur.deg$id==ids[i],]$y.born)+2, 3),3)  # choose age, men tend to be older than the women they marry
      age[i] <- round(rnorm(1, (cur.deg[cur.deg$id==ids[i],]$age)+2, 3),3)  #choose age at death women usually outlive their spouses
    }
    id[i]<-cur.deg[cur.deg$id==ids[i],]$id+0.1  ## id of spouse is id of carrier-mate +0.1
    #dead[i]<-1
    #   if (round(2010-y.born[i]) < round(rnorm(1, 81.1, 5),3)){ dead[i]<-0}
    # code to fix spouse death status and correct age to current age or age at death
    age.temp[i]<-round(2010-y.born[i])
    age[i]<-.demog.nat(y.born[i],female[i],demographics.df)$age.death
    #age[i]<-.demog.nat.china(y.born[i],female[i])$age.death
    dead[i]<-1
    if (age.temp[i]< round(rnorm(1, 81.1, 5),3)){ #took the older age
      age[i]<-age.temp[i]
      dead[i]<-0
    }

    in.2 <- data.frame(list("degree"=deg.num, "momid"=NA, "dadid"=NA, "id"=id, "age"=age, "female"=female, "y.born"=y.born, "dead"=dead, "geno"=geno))

    ##### generate offspring
    in.1 <- cur.deg
    temp<-list("degree"=next.deg, "momid"=NA, "dadid"=NA, "id"=NA, "age"=NA, "female"=NA,"y.born"=NA, "dead"=NA, "geno"=NA)
    if (in.1[in.1$id==ids[i],]$female==1){  # if female parent from the initial pedigree
      momid[i]<-in.1[in.1$id==ids[i],]$id
      dadid[i]<-in.2$id[i]
      fem[i]<-match(1, in.1[in.1$id==ids[i],]$female)  #determine mother and father genotypes
      mal[i]<-match(0, in.2$female[i])
      dad[i] <- in.1[in.1$id==ids[i],]$geno[fem[i]]
      mom[i] <- in.2$geno[mal[i]]
      toffyb <-  (in.1[in.1$id==ids[i],]$y.born) + (.demog.nat(in.1[in.1$id==ids[i],]$y.born, in.1[in.1$id==ids[i],]$female,demographics.df)$age.mar) + 5 ###children clustered around 5 years after age of marriage
      #toffyb <-  (in.1[in.1$id==ids[i],]$y.born) + (.demog.nat.china(in.1[in.1$id==ids[i],]$y.born, in.1[in.1$id==ids[i],]$female)$age.mar) + 5 ###children clustered around 5 years after age of marriage
      temp1<-.demog.2(toffs[i], toffyb, demographics.df) #generate age, gender, no. of offsprings with list
      temp2<-.genotype(mom[i], dad[i], temp1)#generate genotypes of offsprings with list
      temp$momid<-momid[i]
      temp$dadid<-dadid[i]
      temp[c("id", "age", "female", "y.born", "dead", "geno")]<-temp2
      temp<-as.data.frame(temp)
    }else if(in.1[in.1$id==ids[i],]$female==0){
      momid[i]<-in.2$id[i]
      dadid[i]<-in.1[in.1$id==ids[i],]$id
      fem[i]<-match(1, in.2$female[i])  #determine mother father genotypes
      mal[i]<-match(0, in.1[in.1$id==ids[i],]$female)
      dad[i] <- in.2$geno[mal[i]]
      mom[i] <- in.1[in.1$id==ids[i],]$geno[fem[i]]
      toffyb <-  (in.2$y.born[i]) + (.demog.nat(in.2$y.born[i],1,demographics.df)$age.mar) + 5 ###children clustered around 5 years after age of marriage
      #toffyb <-  (in.2$y.born[i]) + (.demog.nat.china(in.2$y.born[i],1)$age.mar) + 5 ###children clustered around 5 years after age of marriage
      temp1<-.demog.2(toffs[i], toffyb, demographics.df) #generate age, gender, no. of offsprings with list
      temp2<-.genotype(mom[i], dad[i], temp1)#generate genotypes of offsprings with list
      temp$momid<-momid[i]
      temp$dadid<-dadid[i]
      temp[c("id", "age", "female", "y.born", "dead", "geno")]<-temp2
      temp<-as.data.frame(temp)
    }

    next.dat<- rbind(next.dat,temp)

  }
### can we remove the NA from the
next.dat$id <- 1:length(next.dat$id)
next.dat$id <- next.dat$id+round(max(prev.dat$id),0)
next.dat.1<-rbind(prev.dat,in.2,next.dat)
return(next.dat.1)
}


.risktoinci <-
function(x,y, nzero = 20, spli = 3){    #x is age (in year)time points, y is cumulative risk at those points
	fit <- lm( y~ns(x, spli) )
	xx <- seq(0,100, length.out=100)
	xxx <- seq(1,101, length.out=100)
	annual <- predict(fit, data.frame(x=xxx)) - predict(fit, data.frame(x=xx))
	plot(x,y, xlim = c(0,100), ylim = c(0,0.9))
	lines(xx, predict(fit, data.frame(x=xx)), col='orange')
	annual <- c(rep(0,nzero),annual[(nzero+1):length(annual)])
	for (i in 1:length(annual)){
		if (annual[i] < 0){
			annual[i] <- 0
			}
		}
	return(annual)
}

.crisk <-
function(age, sex, geno, frequencies.df){
  female=carrier=cancer.type=NULL#this line is here to appease R CMD Check
  cancer.names=unique(frequencies.df$cancer.type)
  number.cancers=length(cancer.names)
  number.ages=length(unique(frequencies.df$age))
  freqs <- matrix(data=0,nrow = number.ages, ncol = number.cancers)  ### mtx of freqs annual incidence imputed from lifetime risk studies.
  aff.result <- rep(0,number.cancers)
  aoo.result <- rep(NA,number.cancers)

  if(age >= 20){ #people below 20 don't have the cancer we are looking at
    sub.frequencies.df=subset(frequencies.df,female==sex & carrier=={geno>=1})
    for(i in 1:number.cancers){
      temp=subset(sub.frequencies.df,cancer.type==cancer.names[i])
      temp=temp[order(temp$age),] #arranges the ages to be ascending
      freqs[,i]=temp$frequencies
    }

    ageT <- round(age)
    if(ageT >=100){
      ageT <- 99
    } ### truncate age at 100 (risks at this age are rough estimates anyway)

    samp <- matrix(nrow = number.ages, ncol = number.ages)

    #here we sample to see whether each individual is affected
    for(h in 1:number.cancers){
      rbfreqs <-   rbinom(number.ages,100,freqs[,h])
      for(i in 1:number.ages){
        samp[i,] <- c(rep(0,(100-rbfreqs[i])),rep(1,rbfreqs[i]))
      }
      for(j in 1:ageT){
        if(sample(x=samp[j,], size = 1) > 0){
          aff.status <- 1
          aff.result[h] <- aff.status
          if(is.na(aoo.result[h] == TRUE)){
            age.onset <- j
            aoo.result[h] <- age.onset
          }
        }
      }
    }
  }
  return(c(aff.result,aoo.result))
}












#This is R code for OriGen made by John Michael O. Ranola ranolaj@uw.edu
#if the function is not to be accessed by users, start with a period(.)

.is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol




.add.parent.cols=function(ped){
	#this function finds all parent row numbers and saves it to the pedigree
	number.people=length(ped$id)
	ped$momrow=0
	ped$dadrow=0
	for(i in 1:number.people){
		temp=which(ped$id==ped$momid[i],arr.ind=TRUE)
		if(length(temp)>0){
			ped$momrow[i]=temp
			ped$dadrow[i]=which(ped$id==ped$dadid[i],arr.ind=TRUE)
		}
	}
	return(ped)
}


.add.pedigree.degree=function(ped){
	#this function modifies the pedigree so that it has the relevant degree information
	#this function uses ped$female vec which is 1 if the individual is female and 0 otherwise
	#this function adds degree information to the pedigree
	if(length(ped$degree)>0){
		#print("Degree is already present. Overwriting degree information")
	}

	#start with first row individual.  Assign degree 1.  Assign all offspring degree 2 and parents degree 0 if any.  Move to all newly assigned individuals and assign degrees that are either +-1 of the current degree of the individual based off parent/offspring relationship.  In the end, add a constant to the degree vector to make the lowest degree 1.
	ped=.add.parent.cols(ped)
	number.people=length(ped$id)
	ped$degree=NA
	ped$degree[1]=1
	individual.done=array(FALSE,dim=c(number.people))
	while(sum(individual.done)<number.people){
		for(i in 1:number.people){
			if(!individual.done[i] & !is.na(ped$degree[i])){
				#sets parents to degree -1
				if(ped$momrow[i]!=0){
					ped$degree[ped$momrow[i]]=ped$degree[i]-1
					ped$degree[ped$dadrow[i]]=ped$degree[i]-1
				}
				#sets offspring to degree +1
				if(ped$female[i]){
					ped$degree[ped$momrow==i]=ped$degree[i]+1
				}else{
					ped$degree[ped$dadrow==i]=ped$degree[i]+1
				}

				individual.done[i]=TRUE
			}
		}
	}
	ped$degree=ped$degree+abs(min(ped$degree))+1
	return(ped)
}



.penetrance.prob=function(genotype,affected.vector,age,gender,gene){
	#fMutant=c(52.3,13.89,0.821),mMutant=c(63.48,12.24,0.021), #John's BRCA1 Estimate
  #fMutant=c(59.83,11.82,0.7396),mMutant=c(61.31,11.96,0.5851), #John's MLH1 Estimate
	#fMutant=c(53,16.5,0.96),mMutant=c(94.5,20,0.0025), #BRCA1 Mohammadi
	#fMutant=c(53.9,16.5,0.96),mMutant=c(94.5,20,0.0025), #BRCA1 Jonker
	#fMutant=c(58.5,13.8,1),mMutant=c(58.5,13.8,0.15), #BRCA2 Mohammadi
	#fNorm=c(64.02,10.38,0.091),mNorm=c(67.29,9.73,0.0015)){ #John's BRCA1 Estimate
  #fNorm=c(65.39,11.54,0.1115),mNorm=c(67.36,10.36,0.1014)){ #John's MLH1 Estimate
	#fNorm=c(72,20,0.15),mNorm=c(94.5,20,0.0025)){ #Mohammadi listed
	#fNorm=c(66.3,14.9,0.08),mNorm=c(94.5,20,0.0025)){ #Jonker model 1
	#fNorm=c(72,16.5,0.10),mNorm=c(94.5,20,0.0025)){ #Jonker model 2
	#mBRCA1 is set to mNorm

	#this function returns the probability that an individual has or doesn't have cancer given their genotype, age, and gender.  Note penetrance is given as (\mu,\sigma,r) for the normal distribution
  fMutant=c(52.3,13.89,0.821)
  mMutant=c(63.48,12.24,0.021) #John's BRCA1 Estimate
  fNorm=c(64.02,10.38,0.091)
  mNorm=c(67.29,9.73,0.0015) #John's BRCA1 Estimate
  if(gene=="MLH1"){
    fMutant=c(59.83,11.82,0.7396)
    mMutant=c(61.31,11.96,0.5851) #John's MLH1 Estimate
    fNorm=c(65.39,11.54,0.1115)
    mNorm=c(67.36,10.36,0.1014) #John's MLH1 Estimate
  }else if(gene=="BRCA2"){
    fMutant=c(54.07,13.30,0.679)
    mMutant=c(57.27,13.43,0.090) #John's BRCA2 Estimate
    fNorm=c(64.02,10.38,0.091)
    mNorm=c(67.29,9.73,0.0015) #John's SEER estimate (same as BRCA1)
  }else if(gene=="BRCA1"){
    fMutant=c(52.3,13.89,0.821)
    mMutant=c(63.48,12.24,0.021) #John's BRCA1 Estimate
    fNorm=c(64.02,10.38,0.091)
    mNorm=c(67.29,9.73,0.0015) #John's BRCA1 Estimate
  }else{
    print("No gene type given or identified. Using BRCA1 for penetrance. ")
  }

	number.people=length(genotype)
	prob=genotype*0
	if(length(affected.vector)!=number.people | length(age)!=number.people | length(gender)!=number.people){
		print("Vectors are different lengths")
		return(0)
	}
	for(i in 1:number.people){
		if(affected.vector[i]==1){#individual is affected
			#for some reason the paper says it is given by the derivative if individual has cancer...
			if(genotype[i]==1){
				if(gender[i]==1){#male
					temp=mMutant[3]*dnorm(age[i],mMutant[1],mMutant[2])/mMutant[2]
				} else {
					temp=fMutant[3]*dnorm(age[i],fMutant[1],fMutant[2])/fMutant[2]
				}
			} else {
				if(gender[i]==1){#male
					temp=mNorm[3]*dnorm(age[i],mNorm[1],mNorm[2])/mNorm[2]
				} else {
					temp=fNorm[3]*dnorm(age[i],fNorm[1],fNorm[2])/mNorm[2]
				}
			}
			prob[i]=temp
		} else if(affected.vector[i]==0){ #individual is unaffected
			if(genotype[i]==1){
				if(gender[i]==1){#male
					temp=mMutant[3]*pnorm(age[i],mMutant[1],mMutant[2])
				} else {
					temp=fMutant[3]*pnorm(age[i],fMutant[1],fMutant[2])
				}
			} else {
				if(gender[i]==1){#male
					temp=mNorm[3]*pnorm(age[i],mNorm[1],mNorm[2])
				} else {
					temp=fNorm[3]*pnorm(age[i],fNorm[1],fNorm[2])
				}
			}
			prob[i]=1-temp
		} else{ #individual has unknown phentoype
      if(genotype[i]==1){
				if(gender[i]==1){#male
          temp1=mMutant[3]*dnorm(age[i],mMutant[1],mMutant[2])/mMutant[2]
					temp2=mMutant[3]*pnorm(age[i],mMutant[1],mMutant[2])
				} else {
          temp1=fMutant[3]*dnorm(age[i],fMutant[1],fMutant[2])/fMutant[2]
					temp2=fMutant[3]*pnorm(age[i],fMutant[1],fMutant[2])
				}
			} else {
				if(gender[i]==1){#male
          temp1=mNorm[3]*dnorm(age[i],mNorm[1],mNorm[2])/mNorm[2]
					temp2=mNorm[3]*pnorm(age[i],mNorm[1],mNorm[2])
				} else {
          temp1=fNorm[3]*dnorm(age[i],fNorm[1],fNorm[2])/mNorm[2]
					temp2=fNorm[3]*pnorm(age[i],fNorm[1],fNorm[2])
				}
			}
			prob[i]=0.5*{temp1}+0.5*{1-temp2}
    }
	}
	return(prob)
}



.RemoveUnconnectedIndividuals=function(ped){
  #this function removes people from the pedigree who have no kids or parents in the paedigree
  #only founders will not have any parents in the pedigree so look for momid=NA
  number.people=length(ped$id)
  is.parent=ped$id*0

  for(i in 1:number.people){
  	if((sum(ped$momid==ped$id[i],na.rm=TRUE)>0) | (sum(ped$dadid==ped$id[i],na.rm=TRUE)>0)){
  		is.parent[i]=1
  	}else{
  		is.parent[i]=0
  	}
  }
  #print(is.parent)

  for(i in number.people:1){ #we go backwards because it is removing the row numbers...
  	if(is.na(ped$momid[i])){#no parents
  		if(is.parent[i]==0){ #no kids
  			ped=ped[-i,]
  		}
  	}
  }

  row.names(ped)=1:length(ped$id)

  return(ped)
}





CalculateLikelihoodRatio=
function(ped,affected.vector,gene="BRCA1"){
  #ped should have id, momid, dadid, age, y.born, female, geno and or genotype,
  #In this function we calculate the likelihood ratio "on the fly", meaning that we don't save any possible genotype information.  This is done so we could potentially increase the number of genotype we can process.  Currently we can do 25 non-founders because the number of possible genotype then would have a maximum of 2^25 (possible is about 5% that).  This is the limiting array in terms of storage.  If we do away with it then we will be able to do much more though we will now be limited by computing time.

	#Saving needed variables
	number.people=length(ped$id)
	lr.numerator=0
	lr.denominator=0
	check.num.genotype.probability=0
	check.den.genotype.probability=0
	number.genotypes.found=0

	#first we check if the pedigrees have the cols we need
	ped=.add.parent.cols(ped)

  if(length(ped$genotype)==0){
		ped$genotype=ped$geno
    # stop("Error, pedigree has no genotype information.")
	}

	observed.vector=array(FALSE,dim=c(number.people))
	observed.vector[ped$genotype!=2]=TRUE

	#We save this here for meioses counting
	original.observed.vector=observed.vector

	#Here we modify the observed vector and the genotype vector so that people who are married to carriers are known non-carriers (since there can be only one founder in this model.)
	#as a first pass, we go through each individual and make sure that if one of their parents is a carrier, the other is a non-carrier
	for(i in 1:number.people){
		if(ped$momrow[i]>0){ #person is not a founder
			if(observed.vector[ped$momrow[i]]+observed.vector[ped$dadrow[i]]==1){ #This means that exactly one parent is observed.
				if(ped$genotype[ped$momrow[i]]==1){
					ped$genotype[ped$dadrow[i]]=0
					observed.vector[ped$dadrow[i]]=TRUE
				} else if(ped$genotype[ped$dadrow[i]]==1){
					ped$genotype[ped$momrow[i]]=0
					observed.vector[ped$momrow[i]]=TRUE
				}
			}
		}
	}

	#Here we modify the observed vector and the genotype vector so that if a carrier has a non-carrier parent then the other parent is a carrier.
	for(i in 1:number.people){
		if(ped$genotype[i]==1 & ped$momrow[i]>0){ #non-founder carrier
			if(observed.vector[ped$momrow[i]]+observed.vector[ped$dadrow[i]]==1){ #This means that exactly one parent is observed.
				if(ped$genotype[ped$momrow[i]]==1){
					ped$genotype[ped$dadrow[i]]=0
					observed.vector[ped$dadrow[i]]=TRUE
				}else if(ped$genotype[ped$momrow[i]]==0){
					ped$genotype[ped$dadrow[i]]=1
					observed.vector[ped$dadrow[i]]=TRUE
				}else if(ped$genotype[ped$dadrow[i]]==1){
					ped$genotype[ped$momrow[i]]=0
					observed.vector[ped$momrow[i]]=TRUE
				}else { #(ped$genotype[ped$dadrow[i]]==0)
					ped$genotype[ped$momrow[i]]=1
					observed.vector[ped$momrow[i]]=TRUE
				}
			}
		}
	}


	#add degree information regardless...
	ped=.add.pedigree.degree(ped)

	#Here, we find all ancestors and descendents of each individual
	#Note that ancestor.descendent.array[i,j]=TRUE if j is an ancestor of i or i is a descendent of j
	#Also note that for convenience in the code, ancestor.descendent.array[i,i]=TRUE
	ancestor.descendent.array=array(FALSE,dim=c(number.people,number.people))
	pedigree.founders={ped$momrow==0}
	for(i in 1:number.people){
		ancestor.vec=array(FALSE,dim=c(number.people))
		future.vec=ancestor.vec
		ancestor.vec[i]=TRUE
		current.vec=ancestor.vec
		changes=TRUE
		while(changes){
			for(j in 1:number.people){
				if(current.vec[j] & !pedigree.founders[j]){
					future.vec[ped$dadrow[j]]=TRUE
					future.vec[ped$momrow[j]]=TRUE
				}
			}
			if(sum(future.vec)>0){
				ancestor.vec=ancestor.vec|future.vec
				current.vec=future.vec
				future.vec[]=FALSE
			} else {
				changes=FALSE
			}
		}
		ancestor.descendent.array[i,]=ancestor.vec[]
	}
	#print(ancestor.descendent.array)

	#The methods below depend on the ordering of the pedigree.  The ancestors of the pedigree should have a lower number than the descendents of the pedigree.  Here, we check to see if that's the case and reorder if not.
	counter=1 #note:counter starts at 2 because 1 is trivially true
	attempt.number=1
	order.vec=0
	while(counter<={number.people-2}){
		counter=counter+1
		#here we check if the number of
		if(sum(ancestor.descendent.array[counter,{counter+1}:number.people])>0){
			#then we need to reorder the pedigree
			print("Pedigree misordered.  Reordering...")
			old.ped=ped
			if(attempt.number==1){
				if(length(ped$y.born)>0){
					order.vec=order(old.ped$y.born,old.ped$degree)
				}else{
					order.vec=order(-old.ped$age,old.ped$degree)
				}
				attempt.number=2
			}else if(attempt.number==2){
				order.vec=order(old.ped$degree,-old.ped$age)
				attempt.number=3
			}else{
				print("Error reordering pedigree... Are there loops in the pedigree?  Please contact maintainer.")
				return()
			}
			ped=old.ped[order.vec,]

			#we also need to reorder the ancestor.descendent.array, observed.vector, affected.vector, momrow, dadrow, and pedigree.founders
			print(order.vec)
			ancestor.descendent.array=ancestor.descendent.array[order.vec,order.vec]
			observed.vector=observed.vector[order.vec]
			pedigree.founders=pedigree.founders[order.vec]
			affected.vector=affected.vector[order.vec]
			ped=.add.parent.cols(ped)

			rm(old.ped)
			counter=1#here we set the counter back at 1 to recheck the ordering
		}
	}

	#now, we enumerate all possible phenotype probabilities (2 per individual)
	#note that 1 is not carrier and 2 is carrier
	all.phenotype.probabilities=array(0,dim=c(2,number.people))
	for(i in 1:2){
		all.phenotype.probabilities[i,]=.penetrance.prob(genotype=rep({i-1},times=number.people),affected.vector=affected.vector,age=ped$age,gender={ped$female+1},gene=gene)
	}

	# #We try to speed up the algorithm by saving the ratios of the phenotype probabilities for use later
	# phenotype.probabilities.ratios=all.phenotype.probabilities[1,]/all.phenotype.probabilities[2,]

	#Here we store how many immediate offspring each individual has for the calculation of genotype probabilities.
	number.offspring=0*ped$id
	for(i in 1:number.people){
		if(ped$momrow[i]!=0){
			temp.int=ped$momrow[i]
			number.offspring[temp.int]=number.offspring[temp.int]+1
			temp.int=ped$dadrow[i]
			number.offspring[temp.int]=number.offspring[temp.int]+1
		}
	}

	#Here we take note of the row numbers of the possible founders
	proband.ancestors=ancestor.descendent.array[ped$proband==1,]
	temp.vec={proband.ancestors & pedigree.founders}
	number.proband.founders=sum(temp.vec)
	founder.cols=which(temp.vec,arr.ind=TRUE)#note that these are not the id's but rather the col numbers(which could be the same as id)

	#find the founder of all the observed carriers which is the intersection of all of their ancestors.  If there are multiple, pick the first one.  Then find the lineages going from each observed carrier to the founder to find the minimal pedigree containing all the observed values.  Then we need to chop off the top of the tree

	#finding the common ancestral founder
	temp.vec=array(TRUE,dim=c(number.people))
	for(i in 1:number.people){
		if(ped$genotype[i]==1){
			temp.vec=temp.vec&ancestor.descendent.array[i,]
		}
	}
	if(sum(temp.vec)>0){
		temp.founder=which(temp.vec)[1]
	}else{ #there is no ancestral founder for all the observed carriers...impossible under the model
		print("The observed pedigree does not follow the assumptions of the model.  There can't be a single founder for all the carriers.")
		return()
	}

	#Assuming the temp.founder is the true founder, find all lineages from the observed values to that founder and combine them to make the minimal pedigree containing all the observed carriers.
	temp.vec=ancestor.descendent.array[,temp.founder]
	minimal.affected.carrier.pedigree=array(FALSE,dim=c(number.people))
	minimal.carrier.pedigree=array(FALSE,dim=c(number.people))
	for(i in 1:number.people){
		if(ped$genotype[i]==1){
			temp.lineage=ancestor.descendent.array[i,]&temp.vec
			minimal.carrier.pedigree=minimal.carrier.pedigree|temp.lineage
			if(affected.vector[i]==1){
				minimal.affected.carrier.pedigree=minimal.affected.carrier.pedigree|temp.lineage
			}
		}
	}

	#start from the founder and check if he has 2 offspring that are carriers or if he is an observed carrier.  If not, cut him off, find his carrier offspring and check them.  Repeat until offspring either has 2 carriers or is observed.
  if(sum(affected.vector==1 & ped$genotype==1)>1){ #if there is only one carrier affected then this is pointless.
  	current.top=temp.founder
  	if(current.top>0){
  		while(current.top>0 & !observed.vector[current.top]){
  			#store current.top's direct descendents.
  			temp.vec={ped$momrow==current.top}|{ped$dadrow==current.top}
  			#set NA's in temp.vec to FALSE
  			temp.vec[is.na(temp.vec)]=FALSE
  			temp.int=sum(temp.vec&minimal.affected.carrier.pedigree)
  			if(temp.int>1){#done... no need to continue.  The top of the tree is ok
  				current.top=0
  			break
  			}else if(temp.int==1){
  				minimal.affected.carrier.pedigree[current.top]=FALSE
  				current.top=which(temp.vec&minimal.affected.carrier.pedigree)[1]
  			}else{
  				print("Something went wrong.  Impossible pedigree under assumptions.  Contact maintainer.  Affected.carrier.pedigree")
          return(0)
  			}
  		}
  	}
  }else{
    minimal.affected.carrier.pedigree= {affected.vector==1 & ped$genotype==1}
  }
	#Repeat for minimal.carrier.pedigree
	current.top=temp.founder
	if(current.top>0){
		while(current.top>0 & !observed.vector[current.top]){
			#store current.top's direct descendents.
			temp.vec={ped$momrow==current.top}|{ped$dadrow==current.top}
			#set NA's in temp.vec to FALSE
			temp.vec[is.na(temp.vec)]=FALSE
			temp.int=sum(temp.vec&minimal.carrier.pedigree)
			if(temp.int>1){#done... no need to continue.  The top of the tree is ok
				minimal.carrier.pedigree[current.top]=FALSE #remove the top because this one is unsure...
				current.top=0
			break
			}else if(temp.int==1){
				minimal.carrier.pedigree[current.top]=FALSE
				current.top=which(temp.vec&minimal.carrier.pedigree)[1]
			}else{
				print("Something went wrong.  Impossible pedigree under assumptions.  Contact maintainer")
				return()
			}
		}
	}

	#Here we modify the observed values so that all founders that do not contain all the observed carriers as descendents are non-carriers.
	temp.descendents=minimal.carrier.pedigree
	for(i in 1:number.people){
		if(pedigree.founders[i]& !observed.vector[i]){
			temp.descendents=ancestor.descendent.array[,i]
			if(!all(temp.descendents[minimal.carrier.pedigree])){
				ped$genotype[i]=0
				observed.vector[i]=TRUE
			}
		}
	}

	likelihood.ratio=0
	observed.separating.meioses=sum(minimal.affected.carrier.pedigree)-1

	#print(c("observed.separating.meioses",observed.separating.meioses))
	#R stores ints as 32-bits.  This gives a max value of about 2 billion.  We need a larger version so we store it as 2 32-bits integer.
	number.genotypes.vec=array(0,dim=c(2))
  possible.founders={ped$genotype!=0 & pedigree.founders}
  number.possible.founders=sum(possible.founders)
  possible.founder.cols=which(possible.founders,arr.ind=TRUE)
  # print(c("possible.founders: ",possible.founders))
  #print(c("number.possible.founders: ",number.possible.founders))
  # print(c("possible.founder.cols: ",possible.founder.cols))

	lr.results=.Fortran("likelihood_ratio_main",NumberPeople=as.integer(number.people),NumberProbandFounders=as.integer(number.proband.founders),NumberPossibleFounders=as.integer(number.possible.founders),ObservedSeparatingMeioses=as.integer(observed.separating.meioses),NumberOffspring=as.integer(number.offspring), PedGenotype=as.integer(ped$genotype), FounderCols=as.integer(founder.cols), NumberGenotypesVec=as.integer(number.genotypes.vec),  AllPhenotypeProbabilities=as.double(all.phenotype.probabilities), ProbandAncestors=as.logical(proband.ancestors), ObservedVector=as.logical(observed.vector), AncestorDescendentArray=as.logical(ancestor.descendent.array), MomRow=as.integer(ped$momrow), DadRow=as.integer(ped$dadrow), MinimalObservedPedigree=as.logical(minimal.carrier.pedigree), LikelihoodRatio=as.double(likelihood.ratio),PACKAGE="CoSeg")

	number.genotypes=lr.results$NumberGenotypesVec[1]*2147483647+lr.results$NumberGenotypesVec[2]
	#return(list(likelihood.ratio=likelihood.ratio,reordering=order.vec,separating.meioses=observed.separating.meioses))
		return(list(likelihood.ratio=lr.results$LikelihoodRatio,separating.meioses=observed.separating.meioses, number.genotypes.found=number.genotypes))

}
