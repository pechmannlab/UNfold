# UNfold (2024)
# (c) S Pechmann
# Figure 1

library(ggplot2)
library(igraph)
library(ggraph)
library(reshape2)

setwd("~/UNfold/")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pymol structures and 2d secondary structure diagram 



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
svg("figures/figure1/B_alignment.svg", height=4, width=7)

# load and parse data
aln <- read.table("data/processed/c3718_aln.txt", sep='\t') 
names <- aln[,1]
aln <- aln[,-1]
nr <- nrow(aln)
nc <- ncol(aln)

ss <- as.data.frame(read.table("data/processed/consensus_aln_ss.txt", header=T))
ss$aln <- ss$aln + 1 # adjust python indexing

ssdf <- data.frame()
ssclass <- c()
for (i in unique(ss$element)){
  if (i != "none"){
    current_idx = ss[ss$element==i,1]
    current_a = min(current_idx)
    current_b = max(current_idx)
    current_class = as.vector(ss[ss$aln==current_a,2])
    ssdf <- rbind(ssdf, c(current_a, current_b))
    ssclass <- c(ssclass, current_class)
    }
}
colnames(ssdf) <- c("start", "end")
ssdf$class <- ssclass


# plot sequence alignment
par(mar=c(4,4,2,1))
plot(-5, -5, xlim=c(0, nc), ylim=c(0, nr+5), xlab="", ylab="", main="", pch=15, col="white", axes=F)
axis(1, line=0, cex.axis=1)

for (i in 1:nc){
  if (length(which(ss$aln == i)) > 0 ){
    idx = which(ss$aln == i)
    if (ss$ss[idx] == "loop"){ colo2 = "#777777"}
    else if (ss$ss[idx] == "helix"){ colo2  = "#ff0000"}
    else if (ss$ss[idx] == "strand"){ colo2  = "#b735b7"}
    else { colo2 ="black"}
    }
  else { colo2  = "#DEDEDE"}
  
  for (j in 1:nr){
    if (j <= 11){offset <- 0}
    else {offset <- 0.5}
    
    current_aa = aln[j,i]
    if (current_aa == '-'){ colo <- "white"}		
    if (current_aa != '-'){ colo <- colo2}
    polygon(x=c(i, i+1, i+1, i), y=c(j, j, j+1, j+1), col=colo, border=F)
  }
}

legend(440, 4.7, legend=c("strand", "helix", "loop", "not aligned"), pch=15, col=c("#b735b7", "#ff0000", "#777777", "#DEDEDE"), bty='n', xpd=T, cex=0.75)
mtext("Alignment position", side=1, line=2, cex=1.2)
mtext("Proteins", side=2, line=1.5, cex=1.2)

# add interaction network
iArrows <- igraph:::igraph.Arrows
network <- read.table("data/processed/network_repclusterp.txt")
network <- as.matrix(network) + 1
for (i in 1:nrow(network)){
  iArrows(network[i,1], nr+2, network[i,2], nr+2, h.lwd=1, sh.lwd=1, sh.col="#999999", curve=0.02, width=1, size=0)
}

# add secondary structure diagram
lines(1:nc,rep(nr+2, nc), lwd=2 )
plotSS <- function(mat){
  for (i in 1:nrow(mat)){
    start <- mat[i,1]
    end <- mat[i,2] +1
    if (mat[i,3] == "strand"){ colo="#b735b7" }
    else if (mat[i,3] == "helix") {colo="#ff0000"}
    polygon( x=c(start:end, end:start), y=c(rep(nr+2-0.2, (end-start+1)), rep(nr+2+0.2, (end-start+1))), border=F, col=colo)
  }
}
plotSS(ssdf)

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

svg("figures/figure1/C_contacts.svg", height=4, width=4)

Tf <- as.data.frame(read.table("data/processed/Tf_contacts2.txt", header=T))
Tf2 <- Tf[,c(2, 4,5)]
Tf2.m <- melt(Tf2, id="Tf")


ggplot(Tf2.m, aes(x=Tf, y=value)) + 
  geom_point(color="#222222", shape=19) + 
  geom_smooth(method=lm, aes(color=variable, fill=variable), se=T) + 
  labs(x="Folding temperature Tf", y="Number of contacts") + 
  scale_fill_manual(values=c("#222222", "#29ABE2")) +   
  scale_color_manual(values=c("#222222", "#555555")) + 
  theme_classic() +
  theme(
    text = element_text(size=20), 
    legend.position = c(0.8, 0.3), 
    legend.title = element_blank()
  )

dev.off()


  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

svg("figures/figure1/D_ss.svg", height=4, width=4)

contss <- as.data.frame(read.table("data/processed/contacts_aligned_ss.txt", header=T))
names <- contss[,1]  
contss <- contss[,-1]   #remove protein names

ss <- as.data.frame(read.table("data/processed/consensus_aln_ss.txt", header=T))
ssnorm <- rep(NA, length(colnames(contss)))
for (i in 1:length(colnames(contss))){
  current_ss <- colnames(contss[i])
  current_count <- sum(ss$element == current_ss)
  ssnorm[i] <- current_count
}
contssn <- contss/ssnorm          # normalize by number AAs
contssn <- contssn[,c(1, 7, 2, 3, 8, 4, 9, 5, 10, 6, 11)]       # bring in order
contss.m <- melt(contssn[,c(1:11)])

sstype <- rep(NA, nrow(contss.m))   #add ss type for coloring
for (i in 1:nrow(contss.m)){
  current_ss <- as.vector(contss.m$variable[i])
  current_type <- as.vector(ss[ss$element == current_ss, 2])
  current_type <- current_type[1]
  sstype[i] <- current_type
}
contss.m$type <- sstype


ggplot(contss.m, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=type)) + 
  geom_jitter(size=0.5, alpha=0.6, width=0.2, color="#333333") + 
  labs(y="Rel. number of contacts", x="") + 
  scale_fill_manual(values=c("#FF0000", "#b735b7")) +   
  theme_classic() + 
  coord_flip() + 
  theme(
  text = element_text(size=20), 
  axis.line.y = element_blank(), 
  axis.ticks.y = element_blank(),
  legend.position = c(0.83, 0.83), 
  legend.title = element_blank()
  )

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

svg("figures/figure1/E_hex.svg", height=3, width=4)

clustdist <- as.data.frame(read.table("data/processed/clusters_avgdists.txt", header=T))

ggplot(clustdist, aes(x=dist_aa, y=dist_res)) + 
  geom_hex() + 
  labs(x="Distance AAs", y="Distance contacts") + 
  theme_classic() + 
  theme(
    text = element_text(size=20),
    legend.position = 'none' 
  ) 

dev.off()




svg("figures/figure1/E_hist.svg", height=1.2, width=4)

clustdist2 <- clustdist[clustdist$pdb=="YOR089C",]    # only one set of clusters

ggplot(clustdist2, aes(x=N)) + 
  scale_x_continuous(labels=c(3, 6, 9, 12), breaks=c(3, 6, 9, 12)) + 
  scale_y_continuous(labels=c(0, 10, 20), breaks=c(0, 10, 20)) + 
  geom_histogram(color="white", fill="#597995", size=0.8, breaks=c(1:12)+0.5) + 
  labs(x="Cluster size", y="Count") + 
  theme_classic() +
  theme(
    text = element_text(size=20)
  )

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL F ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pymol structures
