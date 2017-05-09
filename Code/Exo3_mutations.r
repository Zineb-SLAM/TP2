
library(MASS)
library(ggplot2)
library(GGally)
library(ggfortify)
library(ggtern)
library(mclust)
library(factoextra)

mut = read.csv("data/mutations2.csv", header=T, row.names=1)
mut = as.dist(mut, diag=T, upper=T)
#mut
nrow(mut)
summary(mut)
#Calculer une représentation euclidienne des données en d = 2 variables par AFTD, et 
#l’afficher; 
#-->takes a set of dissimilarities and returns a set of points such that the distances 
#between the points are approximately equal to the dissimilarities



#On calculera tout d’abord une représentation des données mutations dans un espace de dimension
#d = 5. On utilisera par la suite la fonction kmeans sur ces données.

mut.aftd = cmdscale(mut, k=5, eig = TRUE)
mut.aftd 

shepard = Shepard(mut, mut.aftd$points)
#png("/Users/zineb/Desktop/Studies/GI05/SY09/SY09_TPs/TP2/Figures/Mutations2_1/shepard2.png")
plot(shepard, asp =1)
abline(0,1)

mutation.kmeans = kmeans(mut.aftd$points, centers=3)
#adjustedRandIndex(mutation.kmeans$cluster, mut)
mat = data.frame(mutation.kmeans$cluster)
mat[,1]
#shepard = Shepard(mut.aftd$points, mutation.kmeans$cluster)
rownames(mut.aftd$points)

ggplot(mut.aftd$points, aes(mut.aftd$points[,1], mut.aftd$points[,2], color= mat[,1], 
                            label= rownames(mut.aftd$points))) + geom_point() 
+                     geom_text(aes(label=rownames(mut.aftd$points)),hjust=0, vjust=0) 
+ scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442"))



N = 1000
var.within = matrix(nrow=N, ncol = 1)
for (i in 1: N)
{
    mutation.kmeans = kmeans(mut.aftd$points, centers=3)
    var.within[i,1] = mutation.kmeans$tot.withinss  
}
min(var.within[,1])

unique(var.within[,1])


