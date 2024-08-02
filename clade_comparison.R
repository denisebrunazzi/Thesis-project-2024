'Comparison of gene presence/absence between Wibowo and the other genes in their clusters'

table_pr_ab = read.table("gene_presence_absence.Rtab", header=T, row.names=1) #from ROARY 
# Dimensions: 66955 rows (the genes) and 170 columns (the genomes)

#We are interested in considering only the two clades containing the 3 Wibowo genomes: Clade E (C3) and Clade A (C6)
clades = read.delim("clade_division_env.txt", header = T, row.names = 1, sep=" ")  
clade_E = row.names(subset(clades, x=="C3")) # 14 genomes, 14 is Wibowo_3
clade_A = row.names(subset(clades, x=="C6")) # 34 genomes, 33 is Wibowo_16, 34 is Wibowo_28

table_E = table_pr_ab[clade_E]  # Dimensions: 66955 rows (the genes) and 14 columns (the genomes of clade E)
table_A = table_pr_ab[clade_A] # Dimensions: 66955 rows (the genes) and 34 columns (the genomes of clade A)


'CLADE E'

'Present in all genomes of Clade E, but absent in Wibowo'
present_clade_E <- apply(table_E[, 1:13], 1, function(row) all(row == 1))  # margin = 1 -> the function is applied across rows
absent_W_E <- table_E[, 14] == 0
pc_aw_E_all <- rownames(table_E)[present_clade_E & absent_W_E] # 466 genes
write.table(pc_aw_E_all, file = "PC_AW_E.txt", row.names = FALSE, col.names = FALSE)
PC_AW_E_nogroup <- pc_aw_E_all[!grepl("group", pc_aw_E_all)] # 225
write.table(PC_AW_E_nogroup, file = "PC_AW_E_nogroup.txt", row.names = FALSE, col.names = FALSE)

'Absent in all genomes of Clade E, but present in Wibowo'
absent_clade_E <- apply(table_E[, 1:13], 1, function(row) all(row == 0)) 
present_W_E <- table_E[, 14] == 1
ac_pw_E_all <- rownames(table_E)[absent_clade_E & present_W_E] # 1307 genes
write.table(ac_pw_E_all, file = "AC_PW_E.txt", row.names = FALSE, col.names = FALSE)
AC_PW_E_nogroup <- ac_pw_E_all[!grepl("group", ac_pw_E_all)] # 67
write.table(AC_PW_E_nogroup, file = "AC_PW_E_nogroup.txt", row.names = FALSE, col.names = FALSE)

'Present in all genomes of Clade E, also Wibowo'
present_E <- apply(table_E, 1, function(row) all(row == 1))
all_E <- rownames(table_E)[present_E]
write.table(all_E, file = "all_E.txt", row.names = FALSE, col.names = FALSE)
all_E_nogroup <- all_E[!grepl("group", all_E)] 
write.table(all_E_nogroup, file = "all_E_nogroup.txt", row.names = FALSE, col.names = FALSE)


'CLADE A'

'Present in all genomes of Clade A, but absent in Wibowo'
present_clade_A <- apply(table_A[, 1:32], 1, function(row) all(row == 1)) 
absent_W_A <- apply(table_A[, 33:34], 1, function(row) all(row == 0)) 
pc_aw_A_all <- rownames(table_A)[present_clade_A & absent_W_A] # 297 genes
write.table(pc_aw_A_all, file = "PC_AW_A.txt", row.names = FALSE, col.names = FALSE)
PC_AW_A_nogroup <- pc_aw_A_all[!grepl("group", pc_aw_A_all)] # 136 genes
write.table(PC_AW_A_nogroup, file = "PC_AW_A_nogroup.txt", row.names = FALSE, col.names = FALSE)

'Absent in all genomes of Clade A, but present in Wibowo'
absent_clade_A <- apply(table_A[, 1:32], 1, function(row) all(row == 0)) 
present_W_A <- apply(table_A[, 33:34], 1, function(row) all(row == 1))
ac_pw_A_all <- rownames(table_A)[absent_clade_A & present_W_A] # 42 genes
write.table(ac_pw_A_all, file = "AC_PW_A.txt", row.names = FALSE, col.names = FALSE)
AC_PW_A_nogroup <- ac_pw_A_all[!grepl("group", ac_pw_A_all)] # 1 egene
write.table(AC_PW_A_nogroup, file = "AC_PW_A_nogroup.txt", row.names = FALSE, col.names = FALSE)

'Present in all genomes of Clade A, also Wibowo'
present_A <- apply(table_A, 1, function(row) all(row == 1))
all_A <- rownames(table_A)[present_A]
write.table(all_A, file = "all_A.txt", row.names = FALSE, col.names = FALSE)
all_A_nogroup <- all_A[!grepl("group", all_A)] 
write.table(all_A_nogroup, file = "all_A_nogroup.txt", row.names = FALSE, col.names = FALSE)
