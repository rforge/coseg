add.affected.to.tree <-
function(tree.f,frequencies.df=NULL,g=4){
  #frequencies.df is a data.frame with columns age(int 1-100), cancer.type(char),female(bool), carrier(bool), and frequencies(real).
  if(is.null(tree.f$famid)){
    tree.f$famid=1
  }
  if(is.null(frequencies.df)){
    print("No frequencies given.  Using BRCA1frequencies.df")
    frequencies.df=BRCA1frequencies.df
  }
  risk<-p.risk<-risk.f<-NULL
  size <- nrow(tree.f)
  progress <- 0
  for (i in 1:size){
    ### The next few lines are just to see how this is progressing to make sure it is progressing for larger samples.
    tprogress <- round(i/size, 2)
    if(tprogress > progress){
    	progress <- tprogress
    	print(progress)
    }

    #risk <- crisk1(tree.f$age[i], tree.f$female[i], tree.f$geno[i])  # BRCA1
    #risk <- crisk2(tree.f$age[i], tree.f$female[i], tree.f$geno[i])  # BRCA2
    #risk <- criskS(tree.f$age[i], tree.f$female[i], tree.f$geno[i])   # SEER breast ovarian(population risk, variant is benign)
    #risk <- criskLS(tree.f$age[i], tree.f$female[i], tree.f$geno[i])  # MLH1, MSH2
    #risk <- criskLSS(tree.f$age[i], tree.f$female[i], tree.f$geno[i])  #SEER colon, endometrial, and minor LS tumors
    risk <- .crisk(tree.f$age[i], tree.f$female[i], tree.f$geno[i], frequencies.df)  # general risk function requiring frequencies.df

  p.risk[i]<-list(risk)
  }
  risk.f <- as.data.frame(do.call("rbind", p.risk))
  #colnames(risk.f)<-c('brst.d','ovar.d',"brst.aoo", "ovar.aoo")  # use with crisk1, crisk2, and criskS
  #colnames(risk.f)<-c('crc.d','endo.d','minor.d','crc.aoo','endo.aoo','minor.aoo')  # use with criskLS and criskLSS
  cancer.names=unique(frequencies.df$cancer.type)
  colnames(risk.f)<-c(paste0(cancer.names,".d"),paste0(cancer.names,".aoo"))  # use with .crisk (general method)
  tree.f2<-cbind(tree.f, risk.f)
  #return(tree.f2)

  #### choosing proband from only affected, variant carriers based on carrier status, disease status, and generation of the pedigree
  f <- unique(tree.f2$famid[tree.f2$degree==g])
  first <- noproband <- first.temp <- c()

  cancer.names.d=paste0(cancer.names,".d")
  for (i in 1:length(f)) {
    #temp<-subset(subset(subset(subset(subset(subset(tree.f2,brst.d==1 | ovar.d == 1), geno ==1), degree <= g), degree > g-3) , dead == 0), famid==f[i])  ## for BRCA1, BRCA2 and SEER breast/ovar families
    #temp<-subset(subset(subset(subset(subset(subset(tree.f2,crc.d==1 | endo.d == 1), geno ==1), degree <= g), degree > g-3), dead == 0), famid==f[i])  ## for LS and SEER colon/endometrial/LSminor families

    ### for general risk method
    #here we set up temp.text to be cancer1.d==1|cancer2.d==1|cancer3.d...
    temp.text=paste0(cancer.names.d[1],"==1")
    if(length(cancer.names.d)>1){
      for(j in 2:length(cancer.names.d)){
        temp.text=paste0(temp.text,"|",cancer.names.d[j],"==1")
      }
    }
    temp<-subset(subset(subset(subset(subset(subset(tree.f2,eval(parse(text=temp.text))), geno ==1), degree <= g), degree > g-3), dead == 0), famid==f[i])

    if (nrow(temp) > 0){
        first.temp <- temp[sample.int(nrow(temp),1),]
        first<-rbind(first, first.temp)
		}
	if (nrow(temp) == 0){
		noproband <- c(noproband, i)
	}
  }

  proband<-ifelse(rownames(tree.f2) %in% rownames(first),1,0)
  tree.f2<-cbind(tree.f2, proband)
  for (i in 1:nrow(tree.f2)){

    if(tree.f2$famid[i] %in% noproband){tree.f2$proband[i] <- -1}
    }
  print(noproband)
return(tree.f2)
}
