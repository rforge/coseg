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


.add.genotype=function(ped){
	#this function converts geno to genotype where geno = 1 if carrier and 0 if unknown and genotype = 0 if not carrier, 1 if carrier, and 2 if unknown
	ped$genotype=ped$geno
	ped$genotype[ped$geno==0]=2
	return(ped)
}

.add.pedigree.degree=function(ped){
	#this function modifies the pedigree so that it has the relevant degree information
	#this function uses ped$female vec which is 1 if the individual is female and 0 otherwise
	#this function adds degree information to the pedigree
	if(length(ped$degree)>0){
		print("Degree is already present. Overwriting degree information")
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



.BRCA.penetrance.prob=function(genotype,cancer,age,gender,
	fBRCA=c(52.3,13.89,0.821),mBRCA=c(63.48,12.24,0.021), #John's Estimate
	#fBRCA=c(53,16.5,0.96),mBRCA=c(94.5,20,0.0025), #BRCA1 Mohammadi
	#fBRCA=c(53.9,16.5,0.96),mBRCA=c(94.5,20,0.0025), #BRCA1 Jonker
	#fBRCA=c(58.5,13.8,1),mBRCA=c(58.5,13.8,0.15), #BRCA2 Mohammadi
	fNorm=c(64.02,10.38,0.091),mNorm=c(67.29,9.73,0.0015)){ #John's Estimate
	#fNorm=c(72,20,0.15),mNorm=c(94.5,20,0.0025)){ #Mohammadi listed
	#fNorm=c(66.3,14.9,0.08),mNorm=c(94.5,20,0.0025)){ #Jonker model 1
	#fNorm=c(72,16.5,0.10),mNorm=c(94.5,20,0.0025)){ #Jonker model 2
	#mBRCA1 is set to mNorm

	#this function returns the probability that an individual has or doesn't have cancer given their genotype, age, and gender.  Note penetrance is given as (\mu,\sigma,r) for the normal distribution

	number.people=length(genotype)
	prob=genotype*0
	if(length(cancer)!=number.people | length(age)!=number.people | length(gender)!=number.people){
		print("Vectors are different lengths")
		return()
	}
	for(i in 1:number.people){
		if(cancer[i]==1){
			#for some reason the paper says it is given by the derivative if individual has cancer...
			if(genotype[i]==1){
				if(gender[i]==1){#male
					temp=mBRCA[3]*dnorm(age[i],mBRCA[1],mBRCA[2])/mBRCA[2]
				} else {
					temp=fBRCA[3]*dnorm(age[i],fBRCA[1],fBRCA[2])/fBRCA[2]
				}
			} else {
				if(gender[i]==1){#male
					temp=mNorm[3]*dnorm(age[i],mNorm[1],mNorm[2])/mNorm[2]
				} else {
					temp=fNorm[3]*dnorm(age[i],fNorm[1],fNorm[2])/mNorm[2]
				}
			}
			prob[i]=temp
		} else {
			if(genotype[i]==1){
				if(gender[i]==1){#male
					temp=mBRCA[3]*pnorm(age[i],mBRCA[1],mBRCA[2])
				} else {
					temp=fBRCA[3]*pnorm(age[i],fBRCA[1],fBRCA[2])
				}
			} else {
				if(gender[i]==1){#male
					temp=mNorm[3]*pnorm(age[i],mNorm[1],mNorm[2])
				} else {
					temp=fNorm[3]*pnorm(age[i],fNorm[1],fNorm[2])
				}
			}
			prob[i]=1-temp
		}
	}
	return(prob)
}





CalculateLikelihoodRatio=function(ped,affected.boolean){
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
		ped=.add.genotype(ped)
	}
	observed.vector=array(FALSE,dim=c(number.people))
	observed.vector[ped$genotype!=2]=TRUE
	#print(c("observed.vector", observed.vector))

	# if(length(ped$degree)==0){
		# ped=add.pedigree.degree(ped)
	# }
	#add degree information regardless...
	ped=.add.pedigree.degree(ped)

	#Here, we find all ancestors and descendents of each individual
	#Note that ancestor.descendent.array[i,j]=TRUE if j is an ancestor of i or i is a descendent of j
	#Also note that for convenience in the code, ancestor.descendent.array[i,i]=TRUE
	ancestor.descendent.array=array(FALSE,dim=c(number.people,number.people))
	pedigree.founder={ped$momrow==0}
	for(i in 1:number.people){
		ancestor.vec=array(FALSE,dim=c(number.people))
		future.vec=ancestor.vec
		ancestor.vec[i]=TRUE
		current.vec=ancestor.vec
		changes=TRUE
		while(changes){
			for(j in 1:number.people){
				if(current.vec[j] & !pedigree.founder[j]){
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

			#we also need to reorder the ancestor.descendent.array, observed.vector, affected.boolean, momrow, dadrow, and pedigree.founder
			print(order.vec)
			ancestor.descendent.array=ancestor.descendent.array[order.vec,order.vec]
			observed.vector=observed.vector[order.vec]
			pedigree.founder=pedigree.founder[order.vec]
			affected.boolean=affected.boolean[order.vec]
			ped=.add.parent.cols(ped)

			rm(old.ped)
			counter=1#here we set the counter back at 1 to recheck the ordering
		}
	}



	# print(ped)
	# print("observed.vector")
	# print(observed.vector)
	# print("pedigree.founder")
	# print(pedigree.founder)
	# print("ancestor.descendent.array")
	# print(ancestor.descendent.array)
	#now, we enumerate all possible phenotype probabilities (2 per individual)
	#note that 1 is not carrier and 2 is carrier
	all.phenotype.probabilities=array(0,dim=c(2,number.people))
	for(i in 1:2){
		all.phenotype.probabilities[i,]=.BRCA.penetrance.prob(genotype=rep({i-1},times=number.people),cancer=affected.boolean,age=ped$age,gender={ped$female+1})
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
	# print("proband.ancestors")
	# print(proband.ancestors)
	temp.vec={proband.ancestors & pedigree.founder}
	number.proband.founders=sum(temp.vec)
	founder.cols=which(temp.vec,arr.ind=TRUE)#note that these are not the id's but rather the col numbers(which could be the same as id)

	#Here we calculate the number of meioses that separate all of the observed individuals for use in numerator.genotype.probability later on
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
	minimal.observed.pedigree=array(FALSE,dim=c(number.people))
	for(i in 1:number.people){
		if(ped$genotype[i]==1){
			temp.lineage=ancestor.descendent.array[i,]&temp.vec
			minimal.observed.pedigree=minimal.observed.pedigree|temp.lineage
		}
	}

	#start from the founder and check if he has 2 offspring that are carriers or if he is an observed carrier.  If not, cut him off, find his carrier offspring and check them.  Repeat until offspring either has 2 carriers or is observed.
	current.top=temp.founder
	while(current.top>0 & !observed.vector[current.top]){
		#store current.top's direct descendents.
		temp.vec={ped$momrow==current.top}|{ped$dadrow==current.top}
		#set NA's in temp.vec to FALSE
		temp.vec[is.na(temp.vec)]=FALSE
		temp.int=sum(temp.vec&minimal.observed.pedigree)
		# print(temp.vec)
		# print(minimal.observed.pedigree)
		# print(temp.int)
		if(temp.int>1){#done... no need to continue.  The top of the tree is ok
			current.top=0
		}else if(temp.int==1){
			minimal.observed.pedigree[current.top]=FALSE
			current.top=which(temp.vec&minimal.observed.pedigree)[1]
		}else{
			print("Something went wrong.  Impossible pedigree under assumptions.  Contact maintainer")
			return()
		}
	}

	likelihood.ratio=0
	observed.separating.meioses=sum(minimal.observed.pedigree)-1
	print(c("observed.separating.meioses",observed.separating.meioses))
	#R stores ints as 32-bits.  This gives a max value of about 2 billion.  We need a larger version so we store it as 2 32-bits integer.
	number.genotypes.vec=array(0,dim=c(2))

	lr.results=.Fortran("likelihood_ratio_main",NumberPeople=as.integer(number.people),NumberProbandFounders=as.integer(number.proband.founders),ObservedSeparatingMeioses=as.integer(observed.separating.meioses),NumberOffspring=as.integer(number.offspring), PedGenotype=as.integer(ped$genotype), FounderCols=as.integer(founder.cols), NumberGenotypesVec=as.integer(number.genotypes.vec),  AllPhenotypeProbabilities=as.double(all.phenotype.probabilities), ProbandAncestors=as.logical(proband.ancestors), ObservedVector=as.logical(observed.vector), AncestorDescendentArray=as.logical(ancestor.descendent.array), LikelihoodRatio=as.double(likelihood.ratio),PACKAGE="CoSeg")

	number.genotypes=lr.results$NumberGenotypesVec[1]*2147483647+lr.results$NumberGenotypesVec[2]
	#return(list(likelihood.ratio=likelihood.ratio,reordering=order.vec,separating.meioses=observed.separating.meioses))
		return(list(likelihood.ratio=lr.results$LikelihoodRatio,separating.meioses=observed.separating.meioses, number.genotypes.found=number.genotypes))

}
