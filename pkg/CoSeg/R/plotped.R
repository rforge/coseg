plot.ped <-
function(ped){
	#function that easily plots the pedigrees
	#displays the pedigree...
	# if(length(.genotype.vector)>0 & length(affected.age)>0){
	# 	ped2<-pedigree(id=ped$id,dadid=ped$dadid,momid=ped$momid,sex={ped$female+1},affected=cbind(affected.bool,.genotype.vector,ped$proband))
	# } else if(length(affected.age)>0){
	# 	ped2<-pedigree(id=ped$id,dadid=ped$dadid,momid=ped$momid,sex={ped$female+1},affected=cbind(affected.bool,ped$proband))
	# } else if(length(.genotype.vector)>0){
	# 	ped2<-pedigree(id=ped$id,dadid=ped$dadid,momid=ped$momid,sex={ped$female+1},affected=cbind(.genotype.vector,ped$proband))
	# }else{
	# 	ped2<-pedigree(id=ped$id,dadid=ped$dadid,momid=ped$momid,sex={ped$female+1},affected=cbind(ped$proband))
	# }
  #
	# newid=round(ped$age)
	# if(length(affected.age)>0){
	# 	newid=round(ped$age*{affected.age==0}+affected.age)
	# }
  #There are 11 columns before the disease columns...
  number.affections={length(ped)-11}/2
  affected.status=rep(0,length(ped$id))
  for(i in 1:number.affections){
    affected.status=affected.status+ped[,10+i]
  }

  if(ped$proband[1]==-1){
    proband.vec=0
  }else{
    proband.vec=ped$proband
  }
  ped2<-pedigree(id=ped$id,dadid=ped$dadid,momid=ped$momid,sex={ped$female+1},affected=cbind(proband.vec,ped$geno,affected.status))
	plot(ped2, id=ped$id)
	#title(main=paste0("Pedigree with highlighted proband, genotype, and affected status"))#,sub="Label is age, age at death, or age of onset")
	title(main="Pedigree with highlighted proband, genotype, and affected status")
}
