make.affected.trees <-
function(n = 1,g = 4, gdown = 2,frequencies.df=NULL, demographics.df=NULL){
  #create seen individuals 25 years or older ,
  if(is.null(frequencies.df)){
    print("No frequencies given.  Using BRCA1frequencies.df")
    frequencies.df=BRCA1frequencies.df
  }
  if(is.null(demographics.df)){
    print("No demographics given.  Using USDemographics.df")
    demographics.df=USDemographics.df
  }
  tree.f <- age.prob <- NULL
  while(length(age.prob)<n){
  	age.temp <- rsnorm(1, mean = 51.49, sd = 10, xi =0.8)   #skewed normal distribution for age of seed individual from individuals tested for hereditary cancer at UW
      if (age.temp > 25){
      	age.prob <- c(age.prob,age.temp)
      }   # make sure individual tested was at least 25 years old in 2010
  }
  ## This loop will create the pedigrees calling the make.tree function
  for(i in 1:n){
    print(i)
    t.tree <- make.tree(g, gdown, seed.age=age.prob[i], demographics.df)
    t.tree <- cbind(famid = i, t.tree)
    tree.f <- rbind(tree.f, t.tree)
    #print(tree)
  }

  #return(tree.f)
  tree.f2=add.affected.to.tree(tree.f,frequencies.df,g=g)
return(tree.f2)
}
