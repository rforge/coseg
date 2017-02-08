MakeTree <-
function(g=4, gdown=2, seed.age=50, demographics.df=NULL){

  y.birth <- nat <- NULL

  if(is.null(demographics.df)){
    print("No demographics given.  Using USDemographics.df")
    demographics.df=USDemographics.df
  }

  y.birth[g]<-ceiling(2010-seed.age) 			#seed tested to be 'g'th generation
  fem.prob<-sample(0:1, 1, replace=TRUE)      #randomly assign gender for individual tested to up the pedgree to the founder

  no.possible.proband=TRUE
  counter=0
  while(no.possible.proband){
    #demographic info and birth year for all generations going back from the individual tested
    nat[g]<-list(.demog.nat(as.numeric(y.birth[g]), fem.prob, demographics.df))  #age marriage, death age, family size for individual tested from national data "deg.1.d
    #nat[g]<-list(.demog.nat.china(as.numeric(y.birth[g]), fem.prob))
    ## generate birth year for ancestors going up the pedigree
    eeg <- g-1
    for (i in eeg:1){
      y.birth[i]<-list(ceiling(as.numeric(y.birth[[i+1]])-nat[[i+1]]$age.mar))
      nat[i]<-list(.demog.nat(as.numeric(y.birth[[i]]), fem.prob, demographics.df))
      #nat[i]<-list(.demog.nat.china(as.numeric(y.birth[[i]]), fem.prob))
    }
    #print(paste("y.birth=",y.birth))
    #generation 1 ###########
    #generate this ancestor going back from "seed.age" "g" generations choose gender at random
    degree<-1
    momid<-NA
    dadid<-NA
    age<-dead<-geno<-off.s<-age.temp<-female<-y.born<-age.temp<-d.age.limit<-NULL
    female<-sample(0:1,1) #assign random gender to founder
    y.born <-y.birth[[1]] #initialize from y.birth assign birth years for founder individual
    age.temp<-ceiling(2010-y.born)  #age if still alive today
    d.age.limit<-round(rnorm(1, 81.1, 5),3)  # assign an upper limit to age # why is this not from demographics
     if (age.temp< d.age.limit){ #took the older age  current death age in 2010
        age<-age.temp
        dead<-0
      }else{
        age<-.demog.nat(y.born,female, demographics.df)$age.death
        #age<-.demog.nat.china(y.born,female)$age.death
        dead<-1
        }
    id <-1
    geno <- 1 # founder has variant by definition
    deg.1<-as.data.frame(cbind(degree, momid, dadid, id, age, female, y.born, dead, geno))

    deg <- deg.1
    tg <- g+gdown
  	for (i in 2:tg){
  		deg <- .grow.p(deg, demographics.df)
  		#print(deg)
  	}

    ############################################
    #here we only output individuals with age>0
    deg=subset(deg,age>0)
    ############################################

    counter=counter+1
    if(sum({!deg$dead}&deg$geno&{deg$age>25})>0 | counter>100){
      no.possible.proband=FALSE
    }
  }

if(counter>100){
  print("Error: no liver person with genotype afeter 100 iterations.")
}
return(deg)

}



MakeTrees <-
function(n = 1,g = 4, gdown = 2, demographics.df=NULL){

  tree.f <- age.prob <- NULL

  #create seen individuals 25 years or older ,
  if(is.null(demographics.df)){
      print("No demographics given.  Using USDemographics.df")
      demographics.df=USDemographics.df
  }

  while(length(age.prob)<n){
      age.temp <- rsnorm(1, mean = 51.49, sd = 10, xi =0.8)   #skewed normal distribution for age of seed individual from individuals tested for hereditary cancer at UW
      if (age.temp > 25){
          age.prob <- c(age.prob,age.temp)
      }   # make sure individual tested was at least 25 years old in 2010
  }
  ## This loop will create the pedigrees calling the MakeTree function
  for(i in 1:n){
      print(i)
      t.tree <- MakeTree(g, gdown, seed.age=age.prob[i], demographics.df)
      t.tree <- cbind(famid = i, t.tree)
      tree.f <- rbind(tree.f, t.tree)
  }
return(tree.f)
}
