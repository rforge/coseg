PlotPedigree <- function(ped, affected.vector=NULL, legend.location="topleft", legend.radius=0.1){
	#function that easily plots the pedigrees

  if(ped$proband[1]==-1){
    proband.vec=0
  }else{
    proband.vec=ped$proband
  }
  if(is.null(affected.vector)){
    ped2<-pedigree(id=ped$id,dadid=ped$dadid,momid=ped$momid,sex={ped$female+1},affected=cbind(proband=proband.vec))
  }else{
    ped2<-pedigree(id=ped$id,dadid=ped$dadid,momid=ped$momid,sex={ped$female+1},affected=cbind(proband=proband.vec,affection={affected.vector==1}))
  }

  label=paste0(ped$id, "\n", ped$age, "\n", ped$geno)
	plot(ped2, id=label)
	#title(main=paste0("Pedigree with highlighted proband, genotype, and affected status"))#,sub="Label is age, age at death, or age of onset")
	title(main="Pedigree with highlighted proband and affected status", sub="Label is ID, age, and genotype.  Age is current age, age at death, or age of onset.")
  pedigree.legend(ped2, location=legend.location, radius=legend.radius)
}
