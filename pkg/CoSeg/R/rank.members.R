
rank.members=function(ped,affected.vector,gene="BRCA1"){
  #ped should have id, momid, dadid, age, y.born, female, geno and or genotype,
  #In this function we rank the members of the pedigree with unknown genotype according to how much the likelihood ratio changes if this person were to be genotyped.  Note that we take the average of them being a carrier and non-carrier.
  number.people=length(ped$id)
  #first we check if the pedigrees have the cols we need
  ped=.add.parent.cols(ped)
  if(length(ped$genotype)==0){
    # ped=.add.genotype(ped)
    print("Error: Pedigree has no genotype information.")
    return(0)
  }

  ped$genotype=analyze.pedigree.genotypes(ped)
  number.unknown.genotypes=sum(ped$genotype==2)
  unknown.genotype.positions=which(ped$genotype==2,arr.ind=TRUE)

  if(sum(ped$genotype==2)==0){
    print("Error: No members with unknown genotype to rank.")
    return(0)
  }

  # #The first thing we need to do is find all possible carriers given the known carriers and non-carriers.
  # #We first find all possible founders whose descendents contain all the known carriers
  #
  # #Here, we find all ancestors and descendents of each individual
  # #Note that ancestor.descendent.array[i,j]=TRUE if j is an ancestor of i or i is a descendent of j
  # #Also note that for convenience in the code, ancestor.descendent.array[i,i]=TRUE
  # ancestor.descendent.array=array(FALSE,dim=c(number.people,number.people))
  # pedigree.founders={ped$momrow==0}
  # for(i in 1:number.people){
  #   ancestor.vec=array(FALSE,dim=c(number.people))
  #   future.vec=ancestor.vec
  #   ancestor.vec[i]=TRUE
  #   current.vec=ancestor.vec
  #   changes=TRUE
  #   while(changes){
  #     for(j in 1:number.people){
  #       if(current.vec[j] & !pedigree.founders[j]){
  #         future.vec[ped$dadrow[j]]=TRUE
  #         future.vec[ped$momrow[j]]=TRUE
  #       }
  #     }
  #     if(sum(future.vec)>0){
  #       ancestor.vec=ancestor.vec|future.vec
  #       current.vec=future.vec
  #       future.vec[]=FALSE
  #     } else {
  #       changes=FALSE
  #     }
  #   }
  #   ancestor.descendent.array[i,]=ancestor.vec[]
  # }
  #
  # #find all possible founders containing all carriers
  # possible.founders=ancestor.descendent.array[ped$proband==1,]
  # for(i in 1:number.people){
  #   if(ped$genotype[i]==1){
  #     possible.founders=possible.founders & ancestor.descendent.array[i,]
  #   }
  # }
  # #this line makes sure that there aren't any non-carrier founders in the list.
  # possible.founders=possible.founders & {ped$genotype!=0}
  #
  #
  # #Here we go throgh the unknown genotypes finding possible carriers.  If an individual is a possible carrier(has a carrier parent), we set them to carrier and repeat.
  # temp.genotype=ped$genotype
  # temp.genotype[possible.founders]=1 #sets all the possible founders to carrier status
  # changes=TRUE
  # while(changes){
  #   changes=FALSE
  #   for(i in 1:number.people){
  #     if(temp.genotype[i]==2){
  #       #check if parent is a carrier
  #       if(temp.genotype[ped$momrow[i]]==1 | temp.genotype[ped$dadrow[i]]==1){
  #         temp.genotype[i]=1
  #         changes=TRUE
  #       }
  #     }
  #   }
  # }
  # temp.positions=0*temp.genotype
  # counter=1
  # for(i in 1:number.people){
  #   if(ped$genotype[i]==2 & temp.genotype[i]==1){
  #     temp.positions[counter]=i
  #     counter=counter+1
  #   }
  # }
  # number.unknown.genotypes=sum(temp.positions>0)
  # unknown.genotype.positions=temp.positions[1:number.unknown.genotypes]

  # print(c("number.unknown.genotypes: ",number.unknown.genotypes))
  # print(c("unknown.genotype.positions: ", unknown.genotype.positions))

  #here we cycle through all the unknown genotypes making each one in turn a carrier and then a non-carrier and find the likelihood ratio.
  original.lr=calculate.likelihood.ratio(ped,affected.vector,gene=gene)$likelihood.ratio
  temp.results=array(0,dim=c(2,number.unknown.genotypes))#2 rows
  for(i in 1:number.unknown.genotypes){
    for(j in 0:1){#non-carrier, carrier
      temp.ped=ped
      temp.ped$genotype[unknown.genotype.positions[i]]=j
      temp.ped$genotype=analyze.pedigree.genotypes(temp.ped)
      temp.results[j+1,i]=calculate.likelihood.ratio(temp.ped,affected.vector,gene=gene)$likelihood.ratio
    }
  }

  #here we plot the results.
  average.lr.changes=colSums(abs(log10(temp.results)-log10(original.lr)))/2
  temp.changes=ped$id*0
  temp.changes[unknown.genotype.positions]=average.lr.changes
  ped2<-pedigree(id=ped$id,dadid=ped$dadid,momid=ped$momid,sex={ped$female+1},affected=cbind(ped$proband,ped$genotype==1,affected.vector==2))
	plot(ped2, id=paste0(ped$id, "\n", round(ped$age), "\n", round(temp.changes,digits=2)))
	#title(main=paste0("Pedigree with highlighted proband, genotype, and affected status"))#,sub="Label is age, age at death, or age of onset")
	title(main="Pedigree with highlighted proband, carriers, and affection status", sub=paste0("Label is ID, age, and average likelihood ratio change. Original LR: ",round(original.lr,digits=2)))

	return(list(unknown.genotypes=ped$id[unknown.genotype.positions],modified.lr=temp.results,original.lr=original.lr))

}
