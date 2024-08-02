'PHYLOGENETIC TREE VIEW'

'ANI distance and clade division' 
# Read files
library(ape) 
length = read.delim('ANIb_alignment_lengths.tab', header = T, row.names=1) 
identity = read.delim('ANIb_percentage_identity.tab', header=T, row.names=1) 

# Convert identity table in distance table and build the tree on the clusters
distance = matrix(ncol=ncol(identity), nrow=nrow(identity))
for(j in 1:170) {
  for (i in 1:170){
    distance[j,i] <- 1-as.numeric(identity[j,i]) } } 
colnames(distance) = colnames(identity) 
rownames(distance) = colnames(identity) 
tree_cluster = hclust(as.dist(distance), method='ward.D2') 
par(bg="white", mar=c(3, 4, 5, 2))    # A numerical vector of the form c(bottom, left, top, right)
plot(tree_cluster, lwd=2, hang=0.03, main="Phylogenetic tree", cex=0.6)


'Cluster identification'
group1 <- cutree(tree_cluster, k=6) # 5 clades + OS 
group1
for (i in 1:ncol(identity)) 
  {group1[i]=paste("C", group1[i], sep="")} 
group1 = as.factor(group1)
write.table(group1, file = "clade_division_env.txt", append=F, sep=" ", dec=".", row.names=T, col.names=T) 
save <- as.phylo(tree_cluster)
write.tree(save, file = "tree_cluster_1st_try.tre") 

# Identification of clades considering only comparisons with sequence alignments > 500’000bp'
identity_3 = matrix(ncol=ncol(identity), nrow=nrow(identity), data=NA) 
colnames(identity_3) = colnames(identity) 
rownames(identity_3) = rownames(identity) 
for (i in 1:ncol(identity))  {
  for (j in 1:nrow(identity)){
    if (length[j,i]>500000) {identity_3[j,i]=identity[j,i]} } } 
env = read.delim("clade_division_env.txt", header = T, row.names = 1, sep=" ")  
rownames(env) == rownames(identity_3) 
clade=env$x
clade[env$x=="C1"]<-"B"
clade[env$x=="C2"]<-"C"
clade[env$x=="C3"]<-"E"
clade[env$x=="C4"]<-"D"
clade[env$x=="C6"]<-"A"
clade[env$x=="C5"]<-"OS"

# Genetic distances in terms of ANI
A=(1-identity_3[clade=="A",clade=="A"]) [lower.tri(identity_3[clade=="A",clade=="A"],diag = F)] 
B=(1-identity_3[clade=="B",clade=="B"]) [lower.tri(identity_3[clade=="B",clade=="B"],diag = F)] 
C=(1-identity_3[clade=="C",clade=="C"]) [lower.tri(identity_3[clade=="C",clade=="C"],diag = F)] 
D=(1-identity_3[clade=="D",clade=="D"]) [lower.tri(identity_3[clade=="D",clade=="D"],diag = F)] 
E=(1-identity_3[clade=="E",clade=="E"]) [lower.tri(identity_3[clade=="E",clade=="E"],diag = F)] 
AB=as.numeric(1-identity_3[clade=="A",clade=="B"]) 
AC=as.numeric(1-identity_3[clade=="A",clade=="C"]) 
AD=as.numeric(1-identity_3[clade=="A",clade=="D"]) 
AE=as.numeric(1-identity_3[clade=="A",clade=="E"]) 
BC=as.numeric(1-identity_3[clade=="B",clade=="C"]) 
BD=as.numeric(1-identity_3[clade=="B",clade=="D"]) 
BE=as.numeric(1-identity_3[clade=="B",clade=="E"]) 
CD=as.numeric(1-identity_3[clade=="C",clade=="D"]) 
CE=as.numeric(1-identity_3[clade=="C",clade=="E"]) 
DE=as.numeric(1-identity_3[clade=="D",clade=="E"]) 
AOS=as.numeric(1-identity_3[clade=="A",clade=="OS"]) 
BOS=as.numeric(1-identity_3[clade=="B",clade=="OS"]) 
COS=as.numeric(1-identity_3[clade=="C",clade=="OS"]) 
DOS=as.numeric(1-identity_3[clade=="D",clade=="OS"]) 
EOS=as.numeric(1-identity_3[clade=="E",clade=="OS"]) 

'PLOT - ANI DISTANCE'
par(bg="white", mar=c(8,8,8,8), font.axis=2, cex.axis=1, font.lab=2) # A numerical vector of the form c(bottom, left, top, right)
boxplot(A*100, B*100, C*100, D*100, E*100, AB*100, AC*100, AD*100, AE*100, BC*100, BD*100, BE*100, CD*100, CE*100,DE*100, AOS*100, BOS*100, COS*100, DOS*100, EOS*100,
        col = c("green3", "dodgerblue2", "firebrick3", "gold2", "maroon3","forestgreen", "forestgreen", "forestgreen","forestgreen","forestgreen","forestgreen","forestgreen","forestgreen" ,"forestgreen","forestgreen","palevioletred3", "palevioletred3", "palevioletred3", "palevioletred3", "palevioletred3"), 
        main="Whole genome distance", cex.main=2,
        border = c("green3", "dodgerblue2", "firebrick3", "gold2", "maroon3", "forestgreen", "forestgreen", "forestgreen","forestgreen","forestgreen","forestgreen","forestgreen","forestgreen" ,"forestgreen","forestgreen", "palevioletred3", "palevioletred3", "palevioletred3", "palevioletred3", "palevioletred3"), 
        names = c("Clade A", "Clade B", "Clade C", "Clade D", "Clade E","A vs B", "A vs C", "A vs D", "A vs E", "B vs C", "B vs D", "B vs E", "C vs D", "C vs E", "D vs E", "A vs OS", "B vs OS", "C vs OS", "D vs OS", "E vs OS"), 
        medcol = "black", staplecol = "black", whiskcol = "black", medlwd = 3, staplelwd = 2, whisklwd = 2, whisklty=3, outpch = 21, outcex=0.4, las=2, ylim=c(-3,27), yaxt="n", boxlwd=1, axes=F) 

abline(h=6, lty=3, lwd=1.2, col=c("grey"), xpd=F, )  # dotted line representing the 6% ANI distance threshold
axis(2, at=c(-2,0, 5, 10, 15, 20, 25, 27), lab=c("", 0, 5, 10, 15, 20, 25, " "), las=T, pos=0, lwd=2.5, lwd.ticks=0, font.axis=3) 
axis(1, at=c(0, 1, 2, 3, 4, 5,6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 21.5), 
     lab=c(" ", "Clade A", "Clade B", "Clade C", "Clade D", "Clade E","A vs B", "A vs C", "A vs D", "A vs E", "B vs C", "B vs D", "B vs E", "C vs D", "C vs E", "D vs E","A vs OS", "B vs OS", "C vs OS", "D vs OS", "E vs OS", " ", " "), 
     lwd=2.5, lwd.ticks=0, font.axis=2, cex.axis=1, las=2, srt=45, pos=-2) 

mtext("ANI distance (%)", side=2, line=2, font=2, cex=1.5) 
axis(1, at=c(0.7,5.2), lab=FALSE, lwd=1.5, lwd.ticks=0, pos=27.5) 
axis(1, at=c(5.7,15.2), lab=FALSE, lwd=1.5, lwd.ticks=0, pos=27.5) 
axis(1, at=c(15.7,20.2), lab=FALSE, lwd=1.5, lwd.ticks=0, pos=27.5) 
mtext("Inter Clade", side=3, line=-0.5, font=1, cex=1.5) 
text(3,27.5, "Intra Clade", pos=3, font=1, cex=1.5, xpd=T) 
text(18,27.5, "Inter species", pos=3, font=1, cex=1.5, xpd=T) 


'Jaccard dissimilarity index and clade division'
tab = read.table("gene_presence_absence.Rtab", header=T, row.names=1) #from ROARY 
colnames(tab)==rownames(env)
group = as.character(clade) 
remove = grep("OS",group) 
tab=as.matrix(tab)
A=tab[,group=="A"] 
B=tab[,group=="B"] 
C=tab[,group=="C"] 
D=tab[,group=="D"] 
E=tab[,group=="E"] 
tab1=tab[group!="OS",group!="OS"]

library(vegan)
jac_1 = vegdist(t(tab),"jaccard")  # not ordered
jac = jac_no[rownames(env), rownames(env)] # ordered
colnames(jac) == rownames(env) 

# Genetic distances in terms of Jaccard distances
A1 = (jac[clade=="A",clade=="A"]) [lower.tri(jac[clade=="A",clade=="A"], diag=F)] 
B1 = (jac[clade=="B",clade=="B"]) [lower.tri(jac[clade=="B",clade=="B"], diag=F)] 
C1 = (jac[clade=="C",clade=="C"]) [lower.tri(jac[clade=="C",clade=="C"], diag=F)] 
D1 = (jac[clade=="D",clade=="D"]) [lower.tri(jac[clade=="D",clade=="D"], diag=F)] 
E1 = (jac[clade=="E",clade=="E"]) [lower.tri(jac[clade=="E",clade=="E"], diag=F)] 
AB1 = as.numeric(unlist(jac[clade=="A",clade=="B"]))
AC1 = as.numeric(unlist(jac[clade=="A",clade=="C"]))
AD1 = as.numeric(unlist(jac[clade=="A",clade=="D"]))
AE1 = as.numeric(unlist(jac[clade=="A",clade=="E"]))
BC1 = as.numeric(unlist(jac[clade=="B",clade=="C"]))
BD1 = as.numeric(unlist(jac[clade=="B",clade=="D"]))
BE1 = as.numeric(unlist(jac[clade=="B",clade=="E"])) 
CD1 = as.numeric(unlist(jac[clade=="C",clade=="D"]))
CE1 = as.numeric(unlist(jac[clade=="C",clade=="E"]))
DE1 = as.numeric(unlist(jac[clade=="D",clade=="E"]))


'PLOT - JACCARD DISTANCE'
par(bg="white", mar=c(8,10,8,10), font.axis=2, cex.axis=1, font.lab=2) # A numerical vector of the form c(bottom, left, top, right)
boxplot(A1*10, B1*10, C1*10, D1*10, E1*10, AB1*10, AC1*10, AD1*10, AE1*10, BC1*10, BD1*10, BE1*10, CD1*10, CE1*10, DE1*10,
        col = c("green3", "dodgerblue2", "firebrick3", "gold2", "maroon3","forestgreen", "forestgreen", "forestgreen","forestgreen","forestgreen","forestgreen","forestgreen","forestgreen" ,"forestgreen","forestgreen","palevioletred3"), 
        main="Gene content distance", cex.main=2,
        border = c("green3", "dodgerblue2", "firebrick3", "gold2", "maroon3", "forestgreen", "forestgreen", "forestgreen","forestgreen","forestgreen","forestgreen","forestgreen","forestgreen" ,"forestgreen","forestgreen"), 
        names = c("Clade A", "Clade B", "Clade C", "Clade D", "Clade E","A vs B", "A vs C", "A vs D", "A vs E", "B vs C", "B vs D", "B vs E", "C vs D", "C vs E", "D vs E"), 
        medcol = "black", staplecol = "black", whiskcol = "black", medlwd=3, staplelwd=2, whisklwd=2, whisklty=3, outpch=21, outcex=0.4, las=2, ylim=c(0,10), yaxt="n", boxlwd=1, axes=F) 

#Dotted lines - xpd' maintains the lines within the region of the plot
abline(h=2, lty=3, lwd=1.2, col=c("grey"), xpd=F)  
abline(h=4, lty=3, lwd=1.2, col=c("grey"), xpd=F)
abline(h=6, lty=3, lwd=1.2, col=c("grey"), xpd=F)
# Left axis
axis(2, at=c(0, 2, 4, 6, 8, 10), lab=c(0, 0.2, 0.4, 0.6, 0.8, 1), las=T, pos=0, lwd=2.5, lwd.ticks=0, font.axis=3) 
# Bottom axis
axis(1, at=c(0, 1, 2, 3, 4, 5,6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), 
     lab=c(" ", "Clade A", "Clade B", "Clade C", "Clade D", "Clade E","A vs B", "A vs C", "A vs D", "A vs E", "B vs C", "B vs D", "B vs E", "C vs D", "C vs E", "D vs E", ""), 
     lwd=2.5, lwd.ticks=0, font.axis=2, cex.axis=1, las=2, srt=45, srt=45, pos=0) 
#Texts
mtext("Jaccard distance", side=2, line=2.5, font=2, cex=2) 
axis(1, at=c(0.6,5.2), lab=FALSE, lwd=1.5, lwd.ticks=0, pos=10.2) 
axis(1, at=c(5.7,15), lab=FALSE, lwd=1.5, lwd.ticks=0, pos=10.2)  
mtext("                                         Inter Clades", side=3, line=0, font=1, cex=1.5)  # aggiungi spazio prima per centrarlo alla linea che ha sotto
text(3,10.2, "Intra Clade", pos=3, font=1, cex=1.5, xpd=T) 
