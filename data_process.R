library(tidyr)
df<-separate(df,b,into = c("b","c"),sep = "_")

temp_id = data.frame(all_fam_info[which((all_fam_info$Fam_id==all_id2$fam_id[1]) & (all_fam_info$Role_id==all_id2$role_id[1])),])
for (i in 2:318) {
  tryCatch({
  temp = data.frame(all_fam_info[which((all_fam_info$Fam_id==all_id2$fam_id[i]) & (all_fam_info$Role_id==all_id2$role_id[i])),])
  pair_record[i,1] <- i
  pair_record[i,2] <- which((all_fam_info$Fam_id==all_id2$fam_id[i]) & (all_fam_info$Role_id==all_id2$role_id[i]))
  temp_id <- rbind(temp_id, temp)
  }, error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
}

#To ANNOVA analysis between Age and MethSites
#Create Frequency Table between Age and X3
b1 <- xtabs(~Age+X3, data = part_all_info2)
age_range <- max(part_all_info2$Age, na.rm = T)-min(part_all_info2$Age, na.rm = T)
labels <- c("< 20", "20 - 30","30 - 40", "40 - 50", "50 - 60", "60 - 70", ">= 70")
breaks <- c(1,20,30,40,50,60,70,100)
mytable_agecut<- cut(part_all_info2$Age, breaks = breaks, labels = labels, right = TRUE )
df <- as.data.frame(table(Age=mytable_agecut))
#names(df)[1] <- c('Age')
table(Age=mytable_agecut)
part_all_info2['AgeRange'] <- mytable_agecut
b1 <- xtabs(~AgeRange+X3, data = part_all_info2)
b1_prop <- prop.table(b1)
library(gmodels)
b1 <- CrossTable(part_all_info2$AgeRange,part_all_info2$X3)
chisq.test(b1$t)
dtemp <- c(1:252698)
b <- c(1:252698)
for (i in 1:100) {
b1 <- CrossTable(part_all_info2[,18],part_all_info2[,18+i-1])
temp <- chisq.test(b1$t)
b[i] <- temp$p.value
}
for (i in 1:100) {
  b1 <- CrossTable(part_all_info2[,18],part_all_info2[,18+i-1])
  temp <- fisher.test(b1$t)
  dtemp[i] <- temp$p.value
}

#run in server
dtemp <- c(1:252698)
b <- c(1:252698)
for (i in 1:252698) {
  b1 <- CrossTable(part_all_info2[,18],part_all_info2[,18+i-1])
  temp <- fisher.test(b1$t)
  dtemp[i] <- temp$p.value
}
b <- dtemp
for (j in 2:252698) {
  
for (i in 2:252698) {
  b1 <- CrossTable(part_all_info2[,18+j-1],part_all_info2[,18+i-1])
  temp <- fisher.test(b1$t)
  dtemp[i] <- temp$p.value
  
}
  b <- cbind(b, dtemp)
}
write.table(b,file = "meth_sites_fisher_matrix.txt", sep="\t")
