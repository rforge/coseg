AddAffectedToTree <-
function(tree.f,frequencies.df=NULL,g=4,benign.boolean=FALSE){
  #frequencies.df is a data.frame with columns age(int 1-100), cancer.type(char),female(boolean), carrier(boolean), and frequencies(real).
  geno <- degree <- dead <- famid <- NULL #this line is here to appease R CMD Check
  if(is.null(tree.f$famid)){
    tree.f$famid=1
  }
  if(length(unique(tree.f$famid))>1){
    stop("Error: AddAffectedToTree only works on single trees")
  }
  if(is.null(frequencies.df)){
    print("No frequencies given.  Using BRCA1Frequencies.df")
    frequencies.df = BRCA1Frequencies.df
  }

  if(benign.boolean){
    print("Simulating a benign variant.")
  }

  risk<-p.risk<-risk.f<-tree.f2<-NULL
  size <- nrow(tree.f)
  progress <- 0

  no.proband.logical=TRUE
  no.possible.proband=FALSE
  counter=0
  while(no.proband.logical){
    for (i in 1:size){
      ### The next few lines are just to see how this is progressing to make sure it is progressing for larger samples.
      tprogress <- round(i/size, 2)
      if(tprogress > progress){
      	progress <- tprogress
      	print(progress)
      }

      #this is where affection status is calculated...
      #risk <- crisk1(tree.f$age[i], tree.f$female[i], tree.f$geno[i])  # BRCA1
      #risk <- crisk2(tree.f$age[i], tree.f$female[i], tree.f$geno[i])  # BRCA2
      #risk <- criskS(tree.f$age[i], tree.f$female[i], tree.f$geno[i])   # SEER breast ovarian(population risk, variant is benign)
      #risk <- criskLS(tree.f$age[i], tree.f$female[i], tree.f$geno[i])  # MLH1, MSH2
      #risk <- criskLSS(tree.f$age[i], tree.f$female[i], tree.f$geno[i])  #SEER colon, endometrial, and minor LS tumors

      #this allows for simulating a benign variant by setting all individuals to non-carriers.
      if(benign.boolean){
        risk <- .crisk(tree.f$age[i], tree.f$female[i], 0, frequencies.df)  # general risk function requiring frequencies.df
      }else{
        risk <- .crisk(tree.f$age[i], tree.f$female[i], tree.f$geno[i], frequencies.df)  # general risk function requiring frequencies.df
      }

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
      #possible.probands<-subset(subset(subset(subset(subset(subset(tree.f2,brst.d==1 | ovar.d == 1), geno ==1), degree <= g), degree > g-3) , dead == 0), famid==f[i])  ## for BRCA1, BRCA2 and SEER breast/ovar families
      #possible.probands<-subset(subset(subset(subset(subset(subset(tree.f2,crc.d==1 | endo.d == 1), geno ==1), degree <= g), degree > g-3), dead == 0), famid==f[i])  ## for LS and SEER colon/endometrial/LSminor families

      ### for general risk method
      #here we set up temp.text to be cancer1.d==1|cancer2.d==1|cancer3.d...
      temp.text=paste0(cancer.names.d[1],"==1")
      if(length(cancer.names.d)>1){
        for(j in 2:length(cancer.names.d)){
          temp.text=paste0(temp.text,"|",cancer.names.d[j],"==1")
        }
      }
      possible.probands<-subset(subset(subset(subset(subset(subset(tree.f2,eval(parse(text=temp.text))), geno ==1), degree <= g), degree > g-3), dead == 0), famid==f[i])

      if (nrow(possible.probands) > 0){
          first.temp <- possible.probands[sample.int(nrow(possible.probands),1),]
          first<-rbind(first, first.temp)
  		}

    	if (nrow(possible.probands) == 0){
    		noproband <- c(noproband, i)
        no.possible.proband=TRUE
    	}
    }

    proband<-ifelse(rownames(tree.f2) %in% rownames(first),1,0)
    tree.f2<-cbind(tree.f2, proband)
    if(no.possible.proband){
      tree.f2$proband <- -1
    }

    # for (i in 1:nrow(tree.f2)){
    #   if(tree.f2$famid[i] %in% noproband){
    #     tree.f2$proband[i] <- -1
    #   }
    # }
    # print(noproband)

    counter=counter+1
    if(max(tree.f2$proband)==1 | counter>100){
      no.proband.logical=FALSE
    }
    if(counter %% 20 == 0){
      print(c("AddAffectedToTree counter: ",counter))
    }
    # print(c("min(tree.f2$proband):",min(tree.f2$proband)))
    # print(c("tree.f2$proband:",tree.f2$proband))
  }

  if(counter>100){
    print("Error, proband not found for a pedigree after 100 tries.")
  }

  return(tree.f2)
}



AddAffectedToTrees <-
function(tree.f, frequencies.df=NULL,g=4,benign.boolean=FALSE){
  famid=benign=NULL #this line is here to appease R CMD Check
    temp.tree1 <- temp.tree2 <- trees <- NULL
    f <- unique(tree.f$famid)
    for (i in 1:length(f)) {
        temp.tree1<-subset(tree.f, famid==f[i])
        temp.tree2 <- AddAffectedToTree(tree.f=temp.tree1,frequencies.df=frequencies.df,g=g, benign.boolean=benign)
        trees <- rbind(trees,temp.tree2)
    }
    return(trees)
}
