
FormatWebToCoSeg=function(ped){
  #this function takes a pedigree in the format for analyze.myvariant.org which has ordered columns with names that may not be meaningful to one with meaningful names for use in CoSeg
	names(ped)[1]="famid"
	names(ped)[2]="id"
	names(ped)[3]="dadid"
	names(ped)[4]="momid"
  ped$female=ped[,5]*0
  ped$female[ped[,5]==1]=1
	names(ped)[6]="affection"
	names(ped)[7]="age"
	ped$genotype=ped[,9]
	for(i in 1:length(ped$genotype)){
		if(ped[i,8]==0 & ped[i,9]==0){
			ped$genotype[i]=2 #unknown genotype
		}else if(ped[i,8]==2 & ped[i,9]==2){
			ped$genotype[i]=0 #non-carrier
		}else if(ped[i,8]==1 & ped[i,9]==2){
			ped$genotype[i]=1 #carrier
		}else if(ped[i,8]==2 & ped[i,9]==1){
			ped$genotype[i]=1 #carrier
		}else{
			print("Genotypes do not fit model assumptions.  Individuals should be completely typed or untyped.  Also there should not be any individuals homozygous for the mutant allele.")
			return(0)
		}
	}
	names(ped)[10]="proband"
	return(ped)
}
