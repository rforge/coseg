#This is R code for CoSeg made by John Michael O. Ranola ranolaj@uw.edu
#
# .CalculateAncestorDescendentArray=function(ped){
#   #Here, we find all ancestors and descendents of each individual
#   #Note that ancestor.descendent.array[i,j]=TRUE if j is an ancestor of i or i is a descendent of j
#   #Also note that for convenience in the code, ancestor.descendent.array[i,i]=TRUE
#   number.people=length(ped$id)
#   ped=.add.parent.cols(ped)
#
#   ancestor.descendent.array=array(FALSE,dim=c(number.people,number.people))
#   pedigree.founders={ped$momrow==0}
#   for(i in 1:number.people){
#     ancestor.vec=array(FALSE,dim=c(number.people))
#     future.vec=ancestor.vec
#     ancestor.vec[i]=TRUE
#     current.vec=ancestor.vec
#     changes=TRUE
#     while(changes){
#       for(j in 1:number.people){
#         if(current.vec[j] & !pedigree.founders[j]){
#           future.vec[ped$dadrow[j]]=TRUE
#           future.vec[ped$momrow[j]]=TRUE
#         }
#       }
#       if(sum(future.vec)>0){
#         ancestor.vec=ancestor.vec|future.vec
#         current.vec=future.vec
#         future.vec[]=FALSE
#       } else {
#         changes=FALSE
#       }
#     }
#     ancestor.descendent.array[i,]=ancestor.vec[]
#   }
#   return(ancestor.descendent.array)
# }
#
#
# .add.parent.cols=function(ped){
# 	#this function finds all parent row numbers and saves it to the pedigree
# 	number.people=length(ped$id)
# 	ped$momrow=0
# 	ped$dadrow=0
# 	for(i in 1:number.people){
# 		temp=which(ped$id==ped$momid[i],arr.ind=TRUE)
# 		if(length(temp)>0){
# 			ped$momrow[i]=temp
# 			ped$dadrow[i]=which(ped$id==ped$dadid[i],arr.ind=TRUE)
# 		}
# 	}
# 	return(ped)
# }
#
#
# plot.ped=function (ped, affected.vector = NULL, affected.vector2 = NULL, legend.location = "topleft", legend.radius = 0.1){
#   if (ped$proband[1] == -1) {
#     proband.vec = 0
#   }
#   else {
#     proband.vec = ped$proband
#   }
#   if (is.null(affected.vector)) {
#     ped2 <- pedigree(id = ped$id, dadid = ped$dadid, momid = ped$momid, sex = {ped$female+1}, affected = cbind(proband = proband.vec))
#   }
#   else {
#     ped2 <- pedigree(id = ped$id, dadid = ped$dadid, momid = ped$momid, sex = {ped$female+1}, affected = cbind(proband = proband.vec, affection = {affected.vector == 1}, keep = affected.vector2))
#   }
#   label = paste0(ped$id, "\n", ped$age, "\n", ped$geno)
#   plot(ped2, id = label)
#   title(main = "Pedigree with highlighted proband and affected status",
#       sub = "Label is ID, age, and genotype.  Age is current age, age at death, or age of onset.")
#   pedigree.legend(ped2, location = legend.location, radius = legend.radius)
# }


PrunePedigree=function(ped,affected.vector,pruning.level=1){
  #ped should have id, momid, dadid, age, y.born, female, geno and or genotype,
  #In this function we prune the pedigree of any individuals that are not affected, carriers, or first degree relatives of affecteds or carriers.
	original.ped=ped

	#Saving needed variables
	number.people=length(ped$id)
	ped=.add.parent.cols(ped)
  # ped=.add.pedigree.degree(ped)

  if(length(ped$genotype)==0){
		ped$genotype=ped$geno
    # stop("Error, pedigree has no genotype information.")
	}

  TempPed=ped
  AncestorDescendentArray=.CalculateAncestorDescendentArray(ped)
  #Note that ancestor.descendent.array[i,j]=TRUE if j is an ancestor of i or i is a descendent of j

  #Start by keeping all affecteds/carriers
	CarrierOrAffected={affected.vector==1 | ped$genotype==1}
	KeepVec=CarrierOrAffected

	if(pruning.level<=1){
	  #Now include first degree relatives of affecteds/carriers
	  OldKeepVec=KeepVec
	  for(i in 1:number.people){
	    if(OldKeepVec[i]){
				if(ped$momrow[i]!=0){
		      #saving i's parents
		      KeepVec[ped$momrow[i]]=TRUE
		      KeepVec[ped$dadrow[i]]=TRUE
					#saving i's siblings
					SameParents={ped$momrow==ped$momrow[i] & ped$dadrow==ped$dadrow[i]}
					KeepVec[SameParents]=TRUE
				}
	      #saving i's potential children
	      KeepVec[ped$momrow==i]=TRUE
	      KeepVec[ped$dadrow==i]=TRUE
	    }
	  }
	}else{
		#Here we are just keeping the affected or carriers but we include children of a spouse who is affected, if they have no affected/carrier kids, so that they remain connected to the pedigree.
		for(i in 1:number.people){
			if(KeepVec[i] & ped$momrow[i]==0){
				if(sum(KeepVec[ped$momrow==i])+sum(KeepVec[ped$dadrow==i])==0){
					#saving i's potential children
					KeepVec[ped$momrow==i]=TRUE
					KeepVec[ped$dadrow==i]=TRUE
				}
			}
		}
	}



  #Cycle through the pedigree including individuals who will connect the pedigree
  #Start by including all ancestors of the set.
  OldKeepVec=KeepVec
  for(i in 1:number.people){
    if(OldKeepVec[i]==1){
      KeepVec=KeepVec|AncestorDescendentArray[i,]
    }
  }

  #Chop off the top of the pedigree if it has only a single unaffected/non-carrier offspring
  #look at the founders first
  #save a vector of the spouses of all individuals
  Spouses=array(0,number.people)
  for(i in 1:number.people){
    if(ped$momrow[i]!=0){
      Spouses[ped$momrow[i]]=ped$dadrow[i]
      Spouses[ped$dadrow[i]]=ped$momrow[i]
    }
  }
  Tempmomrow=ped$momrow
	changes=TRUE
	counter=0
	while(changes & counter<=number.people){
		counter=counter+1
		changes=FALSE
	  OldKeepVec=KeepVec
	  for(i in 1:number.people){
			#start with the founder couples and cut them off if they don't have either two offspring that are keepers or one offspring that is a CarrierOrAffected
			if(!CarrierOrAffected[i]){
	    #if we are keeping this person, they and their spouse is a founder, and they are female...
		    if(OldKeepVec[i]==1 & Tempmomrow[i]==0 & Spouses[i]!=0 & ped$female[i]==1){
		      #if only one offspring is in KeepVec and it is neither a carrier nor affected then remove
		      if(Tempmomrow[Spouses[i]]==0 & sum(OldKeepVec[Tempmomrow==i])==1){
						if(sum(CarrierOrAffected[Tempmomrow==i])==0){
			        #remove the individual
							KeepVec[i]=0
							Tempmomrow[Tempmomrow==i]=0
							KeepVec[Spouses[i]]=0
							Tempmomrow[Tempmomrow==Spouses[i]]=0
							changes=TRUE
						}
		      }
		    }
			}
	  }
	}

  #add in the spouses of all remaining individuals
  # TempPed=.RemoveUnconnectedIndividuals(TempPed)
  # return(KeepVec)

	if(sum(KeepVec)<=1){
		KeepVec={ped$proband==1}
		KeepVec[ped$momrow[ped$proband==1]]=TRUE
		KeepVec[ped$dadrow[ped$proband==1]]=TRUE
	}

	sub.ped=subset(original.ped,KeepVec==1)
	for(i in 1:length(sub.ped$id)){
		if(!is.na(sub.ped$momid[i])){
			if(sum(sub.ped$id==sub.ped$momid[i])<1){
				sub.ped$momid[i]=NA
			}
		}
		if(!is.na(sub.ped$dadid[i])){
			if(sum(sub.ped$id==sub.ped$dadid[i])<1){
				sub.ped$dadid[i]=NA
			}
		}
	}
	return(sub.ped)
}

#
# familyid=708
# sub.ped=subset(peds,famid==familyid)
# sub.ped$geno=sub.ped$geno-1
# sub.ped$female=sub.ped$female-1
# KeepVec=PrunePedigree(sub.ped,affected.vector={sub.ped$brst.d==2|sub.ped$ovar.d==2})
# plot.ped(sub.ped, affected.vector={sub.ped$brst.d==2|sub.ped$ovar.d==2}, affected.vector2=KeepVec)
#
# KeepVec=PrunePedigree(ped2,affected.vector={ped2$brst.d==2|ped2$ovar.d==2})
# plot.ped(ped2, affected.vector={ped2$brst.d==2|ped2$ovar.d==2}, affected.vector2=KeepVec)
