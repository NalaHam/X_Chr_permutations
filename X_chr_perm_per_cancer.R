#gene_mutations is a df that includes the combined non synonymous mutations for each person 
#we want just the x chromosome from gene_mutations
library(dplyr)

x_chr_mutations <- subset(gene_mutations, gene_mutations$Chromosome == "chrX")

x_chr_mutations$Project <- substr(x_chr_mutations$Project, 6, 9) #keeps coln1 first 12 characters for each row

#for each 

#sum total mutations by gene and gender
gene_x_chr_mut <- data.frame(tapply(x_chr_mutations$Total, INDEX = list(x_chr_mutations$Gene, x_chr_mutations$gender, x_chr_mutations$Project),
                                    FUN = sum, default = 0))


gene_x_chr_mut$Gene <- row.names(gene_x_chr_mut) #make a row with the genes

gene_x_chr_mut <- gene_x_chr_mut[, c(51, 1:50)]

row.names(gene_x_chr_mut) <- c(1:817)

#section off df for each cancer 

ACC <- gene_x_chr_mut[,1:3] #select ACC 

names(ACC)[2:3] <- c("female", "male") #change column names

ACC$total <- rowSums(ACC[,2:3]) #sum mutations for total mutations

ACC["total"][ACC["total"] == 0] <- NA #get totals that =0 to be na

ACC <- ACC[complete.cases(ACC), ] #remove the genes where there is no mutations

ACC$male_ratio <- ACC$male/ACC$total

ACC$female_ratio <- ACC$female/ACC$total

mean(ACC$male_ratio) #0.144
mean(ACC$female_ratio) #0.856

write.csv(ACC, file = "ACC.csv")

##permutations

expected <- sum(ACC$male)/sum(ACC$total) 

ACC$p_value <- NA

row.names(ACC) <- c(1:249)

#Permutation test for X chr and ACC

x_gene_list <- unique(ACC$Gene)

x_gene_list <- x_gene_list[50:249]

df_for_perm <- subset(x_chr_mutations, x_chr_mutations$Project == "ACC")

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
  y <- which(grepl( x , ACC$Gene))
  
  ACC[ y , 7] <- c
}

ACC_sig <- subset(ACC, ACC$p_value <= 0.05)
























