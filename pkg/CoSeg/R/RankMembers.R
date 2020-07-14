#This is R code for CoSeg made by John Michael O. Ranola ranolaj@uw.edu

.CalculateAncestorDescendentArray=function(ped){
  #Here, we find all ancestors and descendents of each individual
  #Note that ancestor.descendent.array[i,j]=TRUE if j is an ancestor of i or i is a descendent of j
  #Also note that for convenience in the code, ancestor.descendent.array[i,i]=TRUE
  number.people=length(ped$id)
  ped=.add.parent.cols(ped)

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
  return(ancestor.descendent.array)
}


.AnalyzePedigreeGenotypes=function(ped){
  #this function figures out which members of the pedigree truly have unknown genotypes according to the assumptions of the model (i.e. single founder).  It returns a vector of length number.people which is 0 for individuals that are known or implied non-carriers, 1 for individuals that are known or implied carriers, and 2 for untyped individuals whose genotype cannot be inferred.  We do this by finding out which individuals may be carriers and which individuals must be carriers.  This means that the rest of the individuals with unknown genotypes must be non-carriers.
  # print("Beginning AnalyzePedigreeGenotypes")
  number.people=length(ped$id)
  #first we make sure the pedigrees have the cols we need
  ped=.add.parent.cols(ped)
  if(length(ped$genotype)<=1){
    # ped=.add.genotype(ped)
    print("Error: Pedigree has no genotype information.")
    return(0)
  }

  if(sum(ped$genotype==2)==0){
    print("No members with unknown genotype to analyze.")
    return(ped$genotype)
  }


  ancestor.descendent.array=.CalculateAncestorDescendentArray(ped)
  pedigree.founders={ped$momrow==0}

  # Let's do an initial filtering for unusual cases.

  # Here we find all non-carriers by propagating all carriers genotype up the pedigree (through carriers and unknown genotype) to capture all possible founder and then propagating down from there to capture all possible carriers.  Those excluded are non carriers.
  #propagating up to find possible founders
  changes=TRUE
  temp=ped$genotype
  while(changes){
    changes=FALSE
    for(i in 1:number.people){
      if(ped$genotype[i]==1){
        if(ped$momrow[i]!=0 & ped$dadrow[i]!=0){
          if(ped$genotype[ped$momrow[i]]==2){
            ped$genotype[ped$momrow[i]]=1
            changes=TRUE
          }
          if(ped$genotype[ped$dadrow[i]]==2){
            ped$genotype[ped$dadrow[i]]=1
            changes=TRUE
          }
        }
      }
      # print(i)
      # print(ped$genotype)
    }
  }
  # propagating down
  # print("down")
  changes=TRUE
  while(changes){
    changes=FALSE
    for(i in 1:number.people){
      if(ped$genotype[i]==1){
        for(j in 1:number.people){ #Slow step... need to optimize
          if(ped$momrow[j]==i | ped$dadrow[j]==i){
            if(ped$genotype[j]==2){
              ped$genotype[j]=1
              changes=TRUE
            }
          }
        }
      }
      # print(i)
      # print(ped$genotype)
    }
  }
  temp[ped$genotype==2]=0
  ped$genotype=temp



  # print("next")
  changes=TRUE
  counter=0
  while(changes){
    changes=FALSE
    # print("test")
    # print(ped$momrow)
    # print(ped$dadrow)
    # print(ped$genotype)
    # Here we set genotype to 0 (non-carrier) if both parents are non-carriers.
    for(i in 1:number.people){
      # print(i)
      if(ped$genotype[i]==2){
        # print("test2")
        if(ped$momrow[i]!=0 & ped$dadrow[i]!=0){
          # print("parent genotypes")
          # print(ped$momrow[i])
          # print(ped$dadrow[i])
          # print(ped$genotype[ped$momrow[i]])
          # print(ped$genotype[ped$dadrow[i]])
          if(ped$genotype[ped$momrow[i]]==0 & ped$genotype[ped$dadrow[i]]==0){
            ped$genotype[i]=0
            # print("test4")
            # print(ped$genotype)
            changes=TRUE
          }
        }
      }
    }
    # print(ped$genotype)
    # print(ped$id)

    # Here if a carrier has one parent of known genotype, we set the other since we are assuming a single founder.
    for(i in 1:number.people){
      if(ped$genotype[i]==1){
        if(ped$momrow[i]>0 & ped$dadrow[i]>0){
          if(ped$genotype[ped$momrow[i]]==1 & ped$genotype[ped$dadrow[i]]==2){
            ped$genotype[ped$dadrow[i]]=0
            changes=TRUE
          }else if(ped$genotype[ped$momrow[i]]==0 & ped$genotype[ped$dadrow[i]]==2){
            ped$genotype[ped$dadrow[i]]=1
            changes=TRUE
          }else if(ped$genotype[ped$momrow[i]]==2 & ped$genotype[ped$dadrow[i]]==1){
            ped$genotype[ped$momrow[i]]=0
            changes=TRUE
          }else if(ped$genotype[ped$momrow[i]]==2 & ped$genotype[ped$dadrow[i]]==0){
            ped$genotype[ped$momrow[i]]=1
            changes=TRUE
          }
        }
      }
    }
    # print(counter)
    # print(ped$genotype)
    # print(ped$id)

  }
  # print(ped$genotype)
  # print("original")
  #find all possible founders containing all carriers
  if(sum(pedigree.founders & ped$genotype==1)>0){
    possible.founders=pedigree.founders & ped$genotype==1
  } else{
    possible.founders=ancestor.descendent.array[ped$proband==1,]
    for(i in 1:number.people){
      if(ped$genotype[i]==1){
        possible.founders=possible.founders & ancestor.descendent.array[i,]
      }
    }
    #this line makes sure that there aren't any non-carrier founders in the list.
    possible.founders=possible.founders & {ped$genotype!=0}
    possible.founders=possible.founders & pedigree.founders
  }
  # print(possible.founders)


  #Here we go through the unknown genotypes finding all possible carriers.  If an individual is a possible carrier(has a carrier parent), we set them to carrier and repeat.
  temp.genotype=ped$genotype
  temp.genotype[possible.founders]=1 #sets all the possible founders to carrier status
  temp.genotype[!possible.founders & pedigree.founders]=0 #sets founders that don't contain all carriers in the lineages to be non-carriers.
  changes=TRUE
  while(changes){
    changes=FALSE
    for(i in 1:number.people){
      if(temp.genotype[i]==2 & ped$momrow[i]>0){
        #check if parent is a carrier
        # print(ped$genotype)
        # print(temp.genotype[ped$momrow[i]])
        # print(temp.genotype[ped$dadrow[i]])
        # print(ped$momrow[i])
        # print(ped$dadrow[i])
        if(temp.genotype[ped$momrow[i]]==1 | temp.genotype[ped$dadrow[i]]==1){
          temp.genotype[i]=1
          changes=TRUE
        }
      }
    }
  }
  temp.positions=0*temp.genotype
  counter=0
  for(i in 1:number.people){
    if(ped$genotype[i]==2 & temp.genotype[i]==1){
      counter=counter+1
      temp.positions[counter]=i
    }
  }
  number.unknown.genotypes=sum(temp.positions>0)
  unknown.genotype.positions=temp.positions[1:number.unknown.genotypes]
  #temp.genotype has all possible carriers set to 1.


  #now we find all individuals who must be carriers according to the model (single founder) and known genotypes.
  #we can do that by going through all the known carriers up to each possible founder and taking the intersection of all these vectors.
  minimal.pedigree=array(1,dim=c(number.people))
  for(i in 1:number.people){
    if(possible.founders[i]==1){
      #find the minimal pedigree assuming this founder.
      temp.minimal.pedigree=array(0,dim=c(number.people))
      for(j in 1:number.people){
        if(ped$genotype[j]==1){
          #this finds the line connecting the carrier to the assumed founder
          temp.lineage=ancestor.descendent.array[,i] & ancestor.descendent.array[j,]
          temp.minimal.pedigree=temp.minimal.pedigree | temp.lineage
        }
      }
      minimal.pedigree=minimal.pedigree & temp.minimal.pedigree
    }
  }
  #minimal.pedigree has all implied carriers set to 1.
  #This is R code for CoSeg made by John Michael O. Ranola ranolaj@uw.edu

  implied.genotype=ped$genotype
  for(i in 1:number.people){
    if(ped$genotype[i]==2){
      if(minimal.pedigree[i]==1){ #must be carrier
        implied.genotype[i]=1
      }else if(temp.genotype[i]==0){ #must be non-carrier
        implied.genotype[i]=0
      }
    }
  }

  return(implied.genotype)
  # print("Ending AnalyzePedigreeGenotypes")

}
# test=AnalyzePedigreeGenotypes(temp2)
# temp2$genotype
# test

DebugPrint=function(ToPrint){
  # print(deparse(substitute(ToPrint)))
  # print(ToPrint)
}

# RenameID=function(ped){
#   number.people=length(ped$id)
#   old.id=ped$id
#   old.momid=ped$momid
#   old.dadid=ped$dadid
#   ped$id=1:number.people
#   for(i in 1:number.people){
#     ped$momid[old.momid==old.id[i]]=i
#     ped$dadid[old.dadid==old.id[i]]=i
#   }
#   return(ped)
# }

RankMembers=function(ped,affected.vector,gene="BRCA1", legend.location="topleft", legend.radius=0.1){
  #ped should have id, momid, dadid, age, y.born, female, geno and or genotype,
  #In this function we rank the members of the pedigree with unknown genotype according to how much the likelihood ratio changes if this person were to be genotyped.  Note that we take the average of them being a carrier and non-carrier.
  number.people=length(ped$id)
  old.ped=ped
  # ped=RenameID(ped)
  #first we check if the pedigrees have the cols we need
  ped=.add.parent.cols(ped)

  if(length(ped$genotype)==0){
    ped$genotype=ped$geno
    # print("Error: Pedigree has no genotype information.")
    # return(0)
  }

  ped$genotype=.AnalyzePedigreeGenotypes(ped)
  number.unknown.genotypes=sum(ped$genotype==2)
  unknown.genotype.positions=which(ped$genotype==2,arr.ind=TRUE)

  if(sum(ped$genotype==2)==0){
    print("Error: No members with unknown genotype to rank.")
    return(0)
  }

  #Note that ancestor.descendent.array[i,j]=TRUE if j is an ancestor of i or i is a descendent of j
  #Also note that for convenience in the code, ancestor.descendent.array[i,i]=TRUE
  #here we find the probability that each individual with an unknown genotype is a carrier.
  ancestor.descendent.array=.CalculateAncestorDescendentArray(ped)
  probability.carrier=array(0,dim=c(2,number.unknown.genotypes))
  tempprob=0
  for(i in 1:number.unknown.genotypes){
    tempint=unknown.genotype.positions[i]
    DebugPrint(tempint)
    DebugPrint(ped$id[tempint])
    for(j in 1:number.people){
      tempprob=0
      if(ped$genotype[j]==1){
        DebugPrint(j)
        DebugPrint(ped$id[j])
        if(ancestor.descendent.array[tempint,j]==1){
          tempprob=2^{1-sum(ancestor.descendent.array[tempint,] & ancestor.descendent.array[,j])}
          DebugPrint(ancestor.descendent.array[tempint,]&ancestor.descendent.array[,j])
        }else if(ancestor.descendent.array[j,tempint]==1){
          tempprob=2^{1-sum(ancestor.descendent.array[j,] & ancestor.descendent.array[,tempint])}
          DebugPrint(ancestor.descendent.array[j,]&ancestor.descendent.array[,tempint])
        }else{
          possible.founders=ancestor.descendent.array[j,] & ancestor.descendent.array[tempint,] & ped$momrow==0
          current.min=2*number.people
          current.count=0
          for(k in 1:number.people){
            if(possible.founders[k]){
              DebugPrint(k)
              DebugPrint(ped$id[k])
              current.count=sum(ancestor.descendent.array[j,]&ancestor.descendent.array[,k])+sum(ancestor.descendent.array[tempint,]&ancestor.descendent.array[,k])
              DebugPrint(ancestor.descendent.array[j,]&ancestor.descendent.array[,k])
              DebugPrint(ancestor.descendent.array[tempint,]&ancestor.descendent.array[,k])
              DebugPrint(current.count)
              current.min=min(current.min,current.count)
            }
          }
          DebugPrint(current.min)
          tempprob=2^{3-current.min}
        }
      }
      if(tempprob>probability.carrier[2,i]){
        probability.carrier[1,i]=1-tempprob #non-carrier
        probability.carrier[2,i]=tempprob #carrier
        DebugPrint(probability.carrier)
      }
    }
  }
  # print(c("tempprob",tempprob))


  #here we cycle through all the unknown genotypes making each one in turn a carrier and then a non-carrier and find the likelihood ratio.
  print(paste0("Iterating through unknowns. ", number.unknown.genotypes, " total"))
  original.lr=CalculateLikelihoodRatio(ped,affected.vector,gene=gene)$likelihood.ratio
  temp.results=array(0,dim=c(2,number.unknown.genotypes))#2 rows
  for(i in 1:number.unknown.genotypes){
    print(c("Current iteration number: ",i))
    for(j in 0:1){#non-carrier, carrier
      temp.ped=ped
      # print("temp.ped$id: ")
      # print(temp.ped$id)
      # print("temp.ped$genotype: before change")
      # print(temp.ped$genotype)
      temp.ped$genotype[unknown.genotype.positions[i]]=j
      # print(paste0(i, ",", j, ","))
      # print("temp.ped$genotype: ")
      # print(temp.ped$genotype)
      temp.ped$genotype=.AnalyzePedigreeGenotypes(temp.ped)
      # print("temp.ped$genotype After analyze: ")
      # print(temp.ped$genotype)
      temp=CalculateLikelihoodRatio(temp.ped,affected.vector,gene=gene)$likelihood.ratio
      # print(paste0(i, ",", j, ",",temp))
      temp.results[j+1,i]=temp
    }
  }

  #here we plot the results.
  # print("temp.results: ")
  # print(temp.results)
  DebugPrint(probability.carrier)
  average.lr.changes=colSums(abs(log10(temp.results*probability.carrier)-log10(original.lr)))/2
  DebugPrint(average.lr.changes)
  temp.changes=array(NA,dim=c(number.people))
  temp.changes[unknown.genotype.positions]=average.lr.changes
  ped2<-pedigree(id=ped$id,dadid=ped$dadid,momid=ped$momid,sex={ped$female+1},affected=cbind(Proband=ped$proband,Carrier={ped$genotype==1},Affected={affected.vector==1}))
	plot(ped2, id=paste0(ped$id, "\n", round(ped$age), "\n", round(temp.changes,digits=2)))
	#title(main=paste0("Pedigree with highlighted proband, genotype, and affected status"))#,sub="Label is age, age at death, or age of onset")
	title(main="Pedigree with highlighted proband, carriers, and affection status", sub=paste0("Label is ID, age, and average likelihood ratio change. Original LR: ",round(original.lr,digits=2)))
  pedigree.legend(ped2, location=legend.location, radius=legend.radius)

	return(list(unknown.genotypes=ped$id[unknown.genotype.positions],modified.lr=temp.results,original.lr=original.lr))

}
