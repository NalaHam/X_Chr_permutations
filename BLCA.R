#section off df for each cancer 

BLCA <- gene_x_chr_mut[,c(1,4,5)] #select BLCA 

names(BLCA)[2:3] <- c("female", "male") #change column names

BLCA$total <- rowSums(BLCA[,2:3]) #sum mutations for total mutations

BLCA["total"][BLCA["total"] == 0] <- NA #get totals that =0 to be na

BLCA <- BLCA[complete.cases(BLCA), ] #remove the genes where there is no mutations

BLCA$male_ratio <- BLCA$male/BLCA$total

BLCA$female_ratio <- BLCA$female/BLCA$total

mean(BLCA$male_ratio) #0.602
mean(BLCA$female_ratio) #0.398

write.csv(BLCA, file = "BLCA.csv")

##permutations

expected <- sum(BLCA$male)/sum(BLCA$total) 

BLCA$p_value <- NA

row.names(BLCA) <- c(1:645)

#Permutation test for X chr and BLCA

x_gene_list <- unique(BLCA$Gene)

x_gene_list <- x_gene_list[50:249]

df_for_perm <- subset(x_chr_mutations, x_chr_mutations$Project == "BLCA")

df_for_perm <- df_for_perm[, c(1,3,6)]

for (x in x_gene_list) {
  data <- subset(df_for_perm, df_for_perm$Gene == x)
  k <- length(data$gender)
  row.names(data) <- c(1:k)
  test.stat1 <- abs((sum(data[which(data$gender=='male'),2 ])/ sum(data$Total)) - 
                      expected) #expect that this should be 0 or close to it. if it is greater than or equal, lets look into it
  #values
  n <- length(data$gender) 
  
  P <- 500000
  
  variable <- data$Total
  
  #make empty df
  PermSamples <- matrix(0, nrow=n, ncol=P)
  
  for(i in 1:P){
    PermSamples[,i] <- sample(variable, size= n, replace=FALSE)
  }
  
  #make empty vector
  Perm.test.stat1 <- rep(0, P)
  
  # loop thru, and calculate the test-stats
  for (j in 1:P){
    # calculate the perm-test-stat1 and save it
    Perm.test.stat1[j] <- abs((sum(PermSamples[which(data$gender=='male'),j ])/ sum(data$Total)) - 
                                expected)
  }
  
  #calculate the p-value, for all P
  c <- mean(Perm.test.stat1 >= test.stat1)
  y <- which(grepl( x , BLCA$Gene))
  
  BLCA[ y , 7] <- c
}

BLCA_sig <- subset(BLCA, BLCA$p_value <= 0.05)
