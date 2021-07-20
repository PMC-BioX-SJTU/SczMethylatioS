#! /usr/bin/env Rscript
all_106_batman <- read.table("all_106fam_batman.txt", header = TRUE)#read datafile
probe.number <- nrow(all_106_batman)
fam.number <- ncol(all_106_batman)
dist_all = array(data = NA, dim = c(probe.number-1, 27))
dist_all[is.na(dist_all)] <- 0
colnames(dist_all) <- c("111","110","11-1","101","100","10-1","1-11","1-10","1-1-1","011","010","01-1","001","000","00-1","0-11","0-10","0-1-1","-111","-110","-11-1","-101","-100","-10-1","-1-11","-1-10","-1-1-1")
dict_comb <- c("111" = 1,
"110" = 2,
"11-1" = 3,
"101" = 4,
"100" = 5,
"10-1" = 6,
"1-11" = 7,
"1-10" = 8,
"1-1-1" = 9,
"011" = 10,
"010" = 11,
"01-1" = 12,
"001" = 13,
"000" = 14,
"00-1" = 15,
"0-11" = 16,
"0-10" = 17,
"0-1-1" = 18,
"-111" = 19,
"-110" = 20,
"-11-1" = 21,
"-101" = 22,
"-100" = 23,
"-10-1" = 24,
"-1-11" = 25,
"-1-10" = 26,
"-1-1-1" =27)
#Define a matrix include all probe and all combination as 111 110 11-1 101 100 10-1 1-11 1-10 1-1-1 011 010 01-1 001 000 00-1 0-11 0-10 0-1-1 -111 -110 -11-1 -101 -100 -10-1 -1-11 -1-10 -1-1-1
#start to calculate each probe's distribution in all families
for (i in 2:10){
for (j in 1:5){
  temp_fa <- all_106_batman[i,(3*(j-1)+1)]
temp_ma <- all_106_batman[i,(3*(j-1)+2)]
temp_of <- all_106_batman[i,(3*(j-1)+3)]
combt <- paste (as.character(temp_fa),as.character(temp_ma),as.character(temp_of),sep="")
D <- dict_comb[as.character(combt)]
if (D == 1) {dist_all[i-1,1] = dist_all[i-1,1] + 1}
else if (D == 2) {dist_all[i-1,2] = dist_all[i-1,2] + 1}
else if (D == 3) {dist_all[i-1,3] = dist_all[i-1,3] + 1}
else if (D == 4) {dist_all[i-1,4] = dist_all[i-1,4] + 1}
else if (D == 5) {dist_all[i-1,5] = dist_all[i-1,5] + 1}
else if (D == 6) {dist_all[i-1,6] = dist_all[i-1,6] + 1}
else if (D == 7) {dist_all[i-1,7] = dist_all[i-1,7] + 1}
else if (D == 8) {dist_all[i-1,8] = dist_all[i-1,8] + 1}
else if (D == 9) {dist_all[i-1,9] = dist_all[i-1,9] + 1}
else if (D == 10) {dist_all[i-1,10] = dist_all[i-1,10] + 1}
else if (D == 11) {dist_all[i-1,11] = dist_all[i-1,11] + 1}
else if (D == 12) {dist_all[i-1,12] = dist_all[i-1,12] + 1}
else if (D == 13) {dist_all[i-1,13] = dist_all[i-1,13] + 1}
else if (D == 14) {dist_all[i-1,14] = dist_all[i-1,14] + 1}
else if (D == 15) {dist_all[i-1,15] = dist_all[i-1,15] + 1}
else if (D == 16) {dist_all[i-1,16] = dist_all[i-1,16] + 1}
else if (D == 17) {dist_all[i-1,17] = dist_all[i-1,17] + 1}
else if (D == 18) {dist_all[i-1,18] = dist_all[i-1,18] + 1}
else if (D == 19) {dist_all[i-1,19] = dist_all[i-1,19] + 1}
else if (D == 20) {dist_all[i-1,20] = dist_all[i-1,20] + 1}
else if (D == 21) {dist_all[i-1,21] = dist_all[i-1,21] + 1}
else if (D == 22) {dist_all[i-1,22] = dist_all[i-1,22] + 1}
else if (D == 23) {dist_all[i-1,23] = dist_all[i-1,23] + 1}
else if (D == 24) {dist_all[i-1,24] = dist_all[i-1,24] + 1}
else if (D == 25) {dist_all[i-1,25] = dist_all[i-1,25] + 1}
else if (D == 26) {dist_all[i-1,26] = dist_all[i-1,26] + 1}
else if (D == 27) {dist_all[i-1,27] = dist_all[i-1,27] + 1}
}
}
write.table(dist_all, "distribution_allcomb.txt", sep ="\t")
print("Done!")
