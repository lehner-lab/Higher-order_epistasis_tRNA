### 
# (A) Circular fitness landscape (Extended Data Figure 1c)
# (B) Linear fitness landscape (Figure 1e)

# JDE, September 2017

# Load data
setwd("007-Landscape_topography/")
load(file = "../001-Data/001-GenotypesFitness.RData")


library(ggplot2)
library(stringdist)
library(parallel)
library(aqfig)
library(geoR)
library(geoR)
library(scales)


# Initiate vars
colors = data.frame(green = "#33A02C", yellow="gold", orange = "#FF7F00", red = "#E31A1C",  purple= "#992ca0", blue = "#1F78B4", brown="#b45b1f", grey="#6f6f6f")
P$id[P$class_var == "WT"] <- ""
vars_P_dna = c("G1A", "T2C", "T2G", "G6A", "G6T", "C27T", "A43C", "C46T", "C66T", "C66A", "A69G", "A70G", "A70T", "C71T")
vars_P = gsub("T", "U", vars_P_dna)


## Generate DATA ##

seq = as.character(P$seq)
id = as.character(P$id)

# Creates a matrix of all possible hamming distances between mutants
m = outer(seq,seq,function(x1,x2){stringdist(x1,x2,method="hamming")})

# For each index (column 1, id_bckg), add at column two the ids of mutants at hamming distance one (the ones to who its connected in the network, id)
edges = setNames(data.frame(do.call("rbind",mclapply(as.list(1:ncol(m)),function(x){
  cbind(id[x],id[m[,x] == 1])
},mc.cores=8))), c("id_bckg", "id"))


# Remove repeated one step mutations and keep the ones where the bckg has equal or smaller HD to the one step genotype
Phylo1 = merge(edges, P[,c("id", "num_vars")])
Phylo1 = merge(Phylo1, P[,c("id", "num_vars")], by.x = "id_bckg", by.y = "id", suffixes = c("","_bckg"))
idx = Phylo1$num_vars_bckg <= Phylo1$num_vars
edges = Phylo1[idx,]


# Get range of colors (Number of breaks per color proportional to the magnitud)
r = range(P$fitness)
rmax = ceiling(r[2]*100)/100
rmin = floor(r[1]*100)/100
n = 10 # min number of breaks
color = c(colorRampPalette(c(as.character(colors[,"red"]),"grey90"))(n*((0 - rmin)/rmax))[1:(n*((0 - rmin)/rmax))-1], colorRampPalette(c("grey90",as.character(colors[,"blue"])))(n))


# Determine positions of dots
# Split P to get one data for each HD to s.cer
sP = split(P,P$num_vars)


# Order by position? If not, alphabetical order.
order_by_pos = F

# Determine position in plot (Circular from center to exterior)
pos = do.call("rbind",lapply(sP,function(x){
  if (!order_by_pos){
    x = x[order(x$id),]
  }
  n = x$num_vars[1]
  s = seq(0,2*pi,length.out = nrow(x)+1)
  pos = setNames(cbind(as.data.frame(t(sapply(1:(length(s)-1),function(y){
    c(n*cos(s[y]), n*sin(s[y]))
  }))), x$id), c("posY","posX","id"))
}))
pos = merge(pos, P[,c("id", "fitness")])

# Merge positions of nodes with each extreme of edges
edges2 = merge(edges,pos,by.x="id",by.y="id") # pos.x are id
edges2 = merge(edges2,pos,by.x="id_bckg",by.y="id") # pos.y are id_bckg


# Find WTs in pos DS
# Species IDs
scer = ""
tb = "U2C,C27U,A43C,C46U,A70G"
cg = "U2G,A70U"
kn = "G1A,G6A,C66U,C71U"
ka = "C46U"
nd = "G6U,C66A"
sc = "G6A,C66U,A69G"
sp_ids = c(scer, tb, cg, kn, ka, nd, sc)

idx_sp = which(pos$id %in% sp_ids)



 ### (A) Extended Data Figure 1c: Circular fitness landscape colouring all edges when mutation U2C is adquired

nbreaks = (abs(rmin) + rmax)/length(color)
idx_eg_U2C = !grepl("U2C", edges2$id_bckg) & grepl("U2C", edges2$id) 

# Plot figure
plot(1,type="n",xlim=c(-max(P$num_vars),max(P$num_vars)),ylim=c(-max(P$num_vars),max(P$num_vars)),axes=F,xlab="",ylab="")
with(edges2[!idx_eg_U2C,],segments(posX.x,posY.x,posX.y,posY.y,col="#00000008"))
with(edges2[idx_eg_U2C,],segments(posX.x,posY.x,posX.y,posY.y,col=alpha("#ffef7f66", 0.2)))
points(pos$posX,pos$posY,col=color[cut(pos$fitness,breaks=seq(rmin,rmax, length.out = length(color)+1),include.lowest = T)],pch=16,cex=0.6)
points(pos$posX[idx_sp],pos$posY[idx_sp],bg=color[cut(pos$fitness[idx_sp],breaks=seq(rmin,rmax, length.out = length(color)+1),include.lowest = T)],pch=21,cex=1, lwd=1.5, col="black")
vertical.image.legend(c(rmin,rmax),col=color)



### (B) Figure 1e: Fitness landscape with fitness on the y axis and HD from S.cer in the x axis showing the distribution of genotypes fitness at each distance.

# Using R
# plot(1,type="n",xlim=c(min(P$num_vars),max(P$num_vars)),ylim=c(min(P$fitness),max(P$fitness)),axes=T,xlab="Hamming distance from S. cerevisiae",ylab="Fitness")
# with(edges2,segments(num_vars,fitness.x,num_vars_bckg,fitness.y,col=alpha("grey50", 0.05)))

# Using ggplot
data_sp = P[P$id %in% sp_ids,]
data_sp$fitness.y = data_sp$fitness
f1e <- ggplot(edges2, aes(x = num_vars_bckg, y=fitness.y)) + ylab(expression(paste("Fitness relative to ", italic("S. cerevisiae")))) + 
  xlab(expression(paste("Hamming distance from ", italic("S. cerevisiae")))) +
  geom_violin(aes(group=num_vars_bckg), scale="width", fill="#FDFCDF", size=0.2, alpha=0.8) +
  geom_segment(aes(xend = num_vars, yend = fitness.x), lineend = "round", size=0.1, alpha=0.08) + 
  geom_point(size=1, aes(color=fitness.y)) + theme_classic() +
  scale_color_gradientn("Fitness", colours = c(as.character(colors[,"red"]),"grey90",as.character(colors[,"blue"])), 
                        values = rescale(c(rmin,0,rmax)),
                        guide = "colorbar", limits=c(rmin,rmax)) +
  scale_fill_gradientn("Fitness", colours = c(as.character(colors[,"red"]),"grey90",as.character(colors[,"blue"])), 
                       values = rescale(c(rmin,0,rmax)),
                       guide = "colorbar", limits=c(rmin,rmax)) +
  guides(color = guide_colorbar(barwidth = 0.8, barheight = 15)) +
  scale_x_continuous(breaks = c(0:10)) +
  geom_point(data=data_sp, aes(x=num_vars, y=fitness, fill=fitness.y), shape=21, color=as.character(colors[,"red"]), size=3) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.line.y = element_line(size=0.2))
f1e


