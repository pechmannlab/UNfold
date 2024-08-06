# UNfold (2024)
# (c) S Pechmann
# SI

library(ggplot2)
library(igraph)
library(ggraph)
library(reshape2)
library(cowplot)
library(ggridges)
library(plyr)

setwd("~/M2/UNfold/")


# Additional supplementary figures not included in the R files for the main figures


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL S1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1. plot RMSD matrix

structaln <- as.data.frame(read.table("data/SI/structaln.txt", header=T))

structaln$Orf1 <- factor(structaln$Orf1, levels=unique(structaln$Orf1))
structaln$Orf2 <- factor(structaln$Orf2, levels=rev(unique(structaln$Orf2)))


svg(file = "figures/SI/rmsd.svg", height = 4, width = 4)

ggplot(structaln, aes(x = Orf1, y = Orf2)) +
  geom_tile(aes(fill=RMS),color="white", size=0.1) + 
  labs(x=NULL, y=NULL) + 
  scale_fill_viridis_c(option="G") + 
  #scale_fill_continuous(type = "viridis") + 
  theme_classic() + 
  theme(
    text = element_text(size=16),
    axis.text.x = element_text(angle=90), 
    axis.line=element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.9, 0.8)
    ) 

dev.off()




svg(file = "figures/SI/rmsd_hist.svg", height = 4, width = 3)

ggplot(structaln, aes(x = RMS)) +
  geom_histogram(color="white", fill="#386e9e") + 
  labs(x="RMSD", y="Density") + 
  theme_classic() + 
  theme(
    text = element_text(size=16)
  ) 

dev.off()





svg(file = "figures/SI/seqid_hist.svg", height = 4, width = 3)

ggplot(structaln, aes(x = SeqID)) +
  geom_histogram(color="white", fill="#5858DD") + 
  labs(x="Sequence identity", y="Density") + 
  theme_classic() + 
  theme(
    text = element_text(size=16)
  ) 

dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL S2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2. pywham illustration + resuls

ex_rough <- as.data.frame(read.table("data/SI/pywham_YBR164C_rough.out", col.names=c("T", "HC")))
ex_refine <- as.data.frame(read.table("data/SI/pywham_YBR164C_refine.out", col.names=c("T", "HC")))
ex_fine <- as.data.frame(read.table("data/SI/pywham_YBR164C_fine.out", col.names=c("T", "HC")))


p1 <- ggplot(ex_rough, aes(x=T, y=HC)) + 
  geom_line() +
  geom_point(col=2, size=2) + 
  scale_x_continuous(limits=c(120, 180)) +
  theme_classic() +
  theme(
    text = element_text(size=16)
  )

p2 <- ggplot(ex_refine, aes(x=T, y=HC)) + 
  geom_line() +  
  geom_point(col=2, size=2) + 
  scale_x_continuous(limits=c(120, 180)) +
  theme_classic() +
  theme(
    text = element_text(size=16)
  )

p3 <- ggplot(ex_fine, aes(x=T, y=HC)) + 
  geom_line() +
  geom_point(col=2, size=2) + 
  scale_x_continuous(limits=c(120, 180)) +
  theme_classic() +
  theme(
    text = element_text(size=16)
  )




svg(file = "figures/SI/A_Tf.svg", height = 5, width = 4)

plot_grid(p1, p2, p3, labels ="", ncol = 1, align = 'v')

dev.off()


#####

# load all files and arrange in ascending order

hc_YBR164C <- as.data.frame(read.table("data/pywham/pywham_YBR164C.out", col.names=c("T", "HC"))) 
hc_YFL005W <- as.data.frame(read.table("data/pywham/pywham_YFL005W.out", col.names=c("T", "HC")))   
hc_YKR014C <- as.data.frame(read.table("data/pywham/pywham_YKR014C.out", col.names=c("T", "HC")))   
hc_YNL093W <- as.data.frame(read.table("data/pywham/pywham_YNL093W.out", col.names=c("T", "HC"))) 
hc_YBR264C <- as.data.frame(read.table("data/pywham/pywham_YBR264C.out", col.names=c("T", "HC")))   
hc_YFL038C <- as.data.frame(read.table("data/pywham/pywham_YFL038C.out", col.names=c("T", "HC")))   
hc_YLR293C <- as.data.frame(read.table("data/pywham/pywham_YLR293C.out", col.names=c("T", "HC")))   
hc_YOR089C <- as.data.frame(read.table("data/pywham/pywham_YOR089C.out", col.names=c("T", "HC"))) 
hc_YCR027C <- as.data.frame(read.table("data/pywham/pywham_YCR027C.out", col.names=c("T", "HC"))) 
hc_YHR022C <- as.data.frame(read.table("data/pywham/pywham_YHR022C.out", col.names=c("T", "HC")))   
hc_YMR138W <- as.data.frame(read.table("data/pywham/pywham_YMR138W.out", col.names=c("T", "HC")))   
hc_YOR094W <- as.data.frame(read.table("data/pywham/pywham_YOR094W.out", col.names=c("T", "HC"))) 
hc_YDL192W <- as.data.frame(read.table("data/pywham/pywham_YDL192W.out", col.names=c("T", "HC")))   
hc_YKL154W <- as.data.frame(read.table("data/pywham/pywham_YKL154W.out", col.names=c("T", "HC")))   
hc_YNL090W <- as.data.frame(read.table("data/pywham/pywham_YNL090W.out", col.names=c("T", "HC")))   
hc_YOR185C <- as.data.frame(read.table("data/pywham/pywham_YOR185C.out", col.names=c("T", "HC"))) 

hc_YDL192W$y <- 1
hc_YBR264C$y <- 2
hc_YNL093W$y <- 3
hc_YMR138W$y <- 4
hc_YCR027C$y <- 5
hc_YOR089C$y <- 6
hc_YHR022C$y <- 7
hc_YLR293C$y <- 8  
hc_YKL154W$y <- 9
hc_YOR185C$y <- 10
hc_YBR164C$y <- 11
hc_YOR094W$y <- 12 
hc_YKR014C$y <- 13
hc_YFL005W$y <- 14
hc_YNL090W$y <- 15
hc_YFL038C$y <- 16

data <- rbind(hc_YBR164C, hc_YFL005W, hc_YKR014C, hc_YNL093W, hc_YBR264C, hc_YFL038C, hc_YLR293C, hc_YOR089C, 
              hc_YCR027C, hc_YHR022C, hc_YMR138W, hc_YOR094W, hc_YDL192W, hc_YKL154W, hc_YNL090W, hc_YOR185C)





svg(file = "figures/SI/Tf_perprotein.svg", height = 4, width = 4)

ggplot(data, aes(x=T, y=y, height = 0.01*HC, group = y)) + 
  geom_ridgeline(aes(fill=y)) +
  #scale_fill_discrete() + 
  labs(x="Temperature", y="Heat capacity") + 
  theme_classic() +
  theme(
    text = element_text(size=20), 
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )

dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL S3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d <- as.data.frame(read.table("data/processed/assign_bpinit.txt", header=T))

countdf <- count(d, c("protein", "class")  )
countdf$class <- factor(countdf$class, levels=c('a2/b2', 'other', 'a1', 'a3', 'a4/b5', 'a5', order ))


svg(file = "figures/SI/unfolding_counts.svg", height = 4, width = 10)

ggplot(countdf, aes(x=protein, y=freq)) + 
  geom_col(aes(fill=class), position=position_dodge2() ) +
  scale_fill_manual(values=c("orange", "#777777", "blue", "red", "darkgreen", "purple" )) + 
  labs(x="", y="Count") +
  theme_classic() + 
  theme(
    text = element_text(size=20),
    axis.text.x = element_text(angle=90),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )

dev.off()





#3. unfolding curves

b2_YDL192W <- as.matrix(read.table("data/SI/ufmat_b2_YDL192W.txt"))
b2_YHR022C <- as.matrix(read.table("data/SI/ufmat_b2_YHR022C.txt"))

b2df <- data.frame(index=c(), folded=c(), pdb=c(), traj=c())
for (i in 1:nrow(b2_YDL192W)){
  tmpdf <- data.frame(index=c(1:ncol(b2_YDL192W)), folded=b2_YDL192W[i,], pdb=rep("YDL192W", ncol(b2_YDL192W)), traj=rep(paste("YDL192W",i, sep="_"), ncol(b2_YDL192W)))
  b2df <- rbind(b2df, tmpdf)
}
for (i in 1:nrow(b2_YHR022C)){
  tmpdf <- data.frame(index=c(1:ncol(b2_YHR022C)), folded=b2_YHR022C[i,], pdb=rep("YHR022C", ncol(b2_YHR022C)), traj=rep(paste("YHR022C",i, sep="_"), ncol(b2_YHR022C)))
  b2df <- rbind(b2df, tmpdf)
}

colfunc_a <- colorRampPalette(c("red", "orange"))
col_a <- colfunc_a(nrow(b2_YDL192W))

colfunc_b <- colorRampPalette(c("#29ABE2", "blue"))
col_b <- colfunc_b(nrow(b2_YHR022C))
col_all <- c(col_a, col_b)


svg(file = "figures/SI/unfolding_b2.svg", height = 4, width = 5)

ggplot(b2df, aes(x=index)) + 
  geom_line(aes(y=folded, color=factor(traj)), size=0.3, show.legend=F) + 
  labs(x="% unfolded protein", y=paste("% unfolded ", "\u03B2", "2", sep=""), title="YDL192W vs. YHR022C" )+
  scale_x_continuous(breaks=c(0, 25, 50, 75,  100), labels=c(0, 25, 50, 75, 100), limits=c(0, 102)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 50, 100)) + 
  scale_color_manual(values=col_all) + 
  geom_line(data=data.frame(pos=c(1:99), avg=as.vector(colMeans(b2_YDL192W))), inherit.aes=F, aes(x=pos, y=avg), col="red", size=1.3) + 
  geom_line(data=data.frame(pos=c(1:99), avg=as.vector(colMeans(b2_YHR022C))), inherit.aes=F, aes(x=pos, y=avg), col="blue", size=1.3) + 
  theme_classic() +
  theme(
    text = element_text(size=18)
  )

dev.off()


#####


svg(file = "figures/SI/unfolding_b2_box.svg", height = 4, width = 3)

boxdf <- data.frame(auc=c(rowSums(b2_YDL192W), rowSums(b2_YHR022C)), orf=c( rep("YDL192W", nrow(b2_YDL192W)), rep("YHR022C", nrow(b2_YHR022C)) )    )
ggplot(boxdf, aes(x=orf, y=auc)) + 
  geom_boxplot(aes(fill=orf), show.legend=F) + 
  labs(x="", y=paste("AUC (\u03B2", "2 unfolding)", sep=""), title="YDL192W vs. YHR022C" ) +
  scale_fill_manual(values = c("red", "blue" )) + 
  theme_classic() +
  theme(
    text = element_text(size=18), 
    axis.text.x = element_text(angle=30),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  )

dev.off()


### -----------------------------------------------------------------
# always reload data, as var names reused!

b2_a2b2 <- read.table("data/SI/ufmat_a2b2_beta2.txt")
b3_a2b2 <- read.table("data/SI/ufmat_a2b2_beta3.txt")
b5_a2b2 <- read.table("data/SI/ufmat_a2b2_beta5.txt")
b6_a2b2 <- read.table("data/SI/ufmat_a2b2_beta6.txt")

df <- data.frame(b2=rowSums(b2_a2b2), b3=rowSums(b3_a2b2), b5=rowSums(b5_a2b2), b6=rowSums(b6_a2b2))
colnames(df) <- c(paste("\u03B2", "2", sep=""), paste("\u03B2", "3", sep=""), paste("\u03B2", "5", sep=""), paste("\u03B2", "6", sep=""))
df.m <- melt(df)

svg(file = "figures/SI/unfolding_a2b2_beta.svg", height = 4, width = 4)

ggplot(df.m, aes(x=variable, y=value)) + 
  geom_violin(aes(fill=variable), show.legend=F) + 
  geom_jitter(size=0.8, alpha=0.5, width=0.2, color="#BBBBBB") + 
  labs(x="", y="AUC (unfolding)", title="a2/b2") +
  scale_fill_manual(values = c("#E527F9", "#E527F9", "#81128D", "#81128D" )) + 
  theme_classic() +
  theme(
    text = element_text(size=18), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  )

dev.off()

boxplot(rowSums(b2_a2b2), rowSums(b3_a2b2), rowSums(b5_a2b2), rowSums(b6_a2b2) )
wilcox.test( c(rowSums(b2_a2b2), rowSums(b3_a2b2)), c(rowSums(b5_a2b2), rowSums(b6_a2b2)) )$p.value


##
b2_a1 <- read.table("data/SI/ufmat_a1_beta2.txt")
b3_a1 <- read.table("data/SI/ufmat_a1_beta3.txt")
b5_a1 <- read.table("data/SI/ufmat_a1_beta5.txt")
b6_a1 <- read.table("data/SI/ufmat_a1_beta6.txt")
boxplot(rowSums(b2_a1), rowSums(b3_a1), rowSums(b5_a1), rowSums(b6_a1) )
wilcox.test( c(rowSums(b2_a1), rowSums(b3_a1)), c(rowSums(b5_a1), rowSums(b6_a1)) )$p.value

##
b2_a3 <- read.table("data/SI/ufmat_a3_beta2.txt")
b3_a3 <- read.table("data/SI/ufmat_a3_beta3.txt")
b5_a3 <- read.table("data/SI/ufmat_a3_beta5.txt")
b6_a3 <- read.table("data/SI/ufmat_a3_beta6.txt")
boxplot(rowSums(b2_a3), rowSums(b3_a3), rowSums(b5_a3), rowSums(b6_a3) )
wilcox.test( c(rowSums(b2_a3), rowSums(b3_a3)), c(rowSums(b5_a3), rowSums(b6_a3)) )$p.value

##
b2_a4b5 <- read.table("data/SI/ufmat_a4b5_beta2.txt")
b3_a4b5 <- read.table("data/SI/ufmat_a4b5_beta3.txt")
b5_a4b5 <- read.table("data/SI/ufmat_a4b5_beta5.txt")
b6_a4b5 <- read.table("data/SI/ufmat_a4b5_beta6.txt")
boxplot(rowSums(b2_a4b5), rowSums(b3_a4b5), rowSums(b5_a4b5), rowSums(b6_a4b5) )
wilcox.test( c(rowSums(b2_a4b5), rowSums(b3_a4b5)), c(rowSums(b5_a4b5), rowSums(b6_a4b5)) )$p.value

##
b2_a5 <- read.table("data/SI/ufmat_a5_beta2.txt")
b3_a5 <- read.table("data/SI/ufmat_a5_beta3.txt")
b5_a5 <- read.table("data/SI/ufmat_a5_beta5.txt")
b6_a5 <- read.table("data/SI/ufmat_a5_beta6.txt")
boxplot(rowSums(b2_a5), rowSums(b3_a5), rowSums(b5_a5), rowSums(b6_a5) )
wilcox.test( c(rowSums(b2_a5), rowSums(b3_a5)), c(rowSums(b5_a5), rowSums(b6_a5)) )$p.value


df <- data.frame(b2=rowSums(b2_a5), b3=rowSums(b3_a5), b5=rowSums(b5_a5), b6=rowSums(b6_a5))
colnames(df) <- c(paste("\u03B2", "2", sep=""), paste("\u03B2", "3", sep=""), paste("\u03B2", "5", sep=""), paste("\u03B2", "6", sep=""))
df.m <- melt(df)

svg(file = "figures/SI/unfolding_a5_beta.svg", height = 4, width = 4)

ggplot(df.m, aes(x=variable, y=value)) + 
  geom_violin(aes(fill=variable), show.legend=F) + 
  geom_jitter(size=0.8, alpha=0.5, width=0.2, color="#BBBBBB") + 
  labs(x="", y="AUC (unfolding)", title="a5") +
  scale_fill_manual(values = c("#E527F9", "#E527F9", "#81128D", "#81128D" )) + 
  theme_classic() +
  theme(
    text = element_text(size=18), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  )

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL S4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4. refolding


unfold_b2 <- as.matrix(read.table("data/SI/ufmat_b2_YDL192W.txt"))
refold_b2 <- t(as.matrix(read.table("data/SI/refoldfmat_b2_YDL192W.txt")))

unfold_b6 <- as.matrix(read.table("data/SI/ufmat_b6_YDL192W.txt"))
refold_b6 <- t(as.matrix(read.table("data/SI/refoldfmat_b6_YDL192W.txt")))

UNdf_b2 <- data.frame(index=c(), folded=c(),traj=c())
for (i in 1:nrow(unfold_b2)){
  tmpdf <- data.frame(index=c(1:ncol(unfold_b2)), folded=unfold_b2[i,], traj=rep(paste("unfold",i, sep="_"), ncol(unfold_b2)))
  UNdf_b2 <- rbind(UNdf_b2, tmpdf)
}
UNdf_b6 <- data.frame(index=c(), folded=c(),traj=c())
for (i in 1:nrow(unfold_b6)){
  tmpdf <- data.frame(index=c(1:ncol(unfold_b6)), folded=unfold_b6[i,], traj=rep(paste("unfold",i, sep="_"), ncol(unfold_b6)))
  UNdf_b6 <- rbind(UNdf_b6, tmpdf)
}


REdf_b2 <- data.frame(index=c(), folded=c(),traj=c())
for (i in 1:nrow(refold_b2)){
  tmpdf <- data.frame(index=c(1:ncol(refold_b2)), folded=refold_b2[i,], traj=rep(paste("refold",i, sep="_"), ncol(refold_b2)))
  REdf_b2 <- rbind(REdf_b2, tmpdf)
}
REdf_b6 <- data.frame(index=c(), folded=c(),traj=c())
for (i in 1:nrow(refold_b6)){
  tmpdf <- data.frame(index=c(1:ncol(refold_b6)), folded=refold_b6[i,], traj=rep(paste("refold",i, sep="_"), ncol(refold_b6)))
  REdf_b6 <- rbind(REdf_b6, tmpdf)
}


colfunc_a <- colorRampPalette(c("red", "orange"))
col_a <- colfunc_a(nrow(unfold_b2))

colfunc_b <- colorRampPalette(c("#29ABE2", "blue"))
col_b <- colfunc_b(nrow(refold_b2))

col_all <- c(col_a, col_b)




svg(file = "figures/SI/refolding.svg", height = 6, width = 6)

p.un_b2 <- ggplot(UNdf_b2, aes(x=index)) + 
  geom_line(aes(y=folded, color=factor(traj)), size=0.3, show.legend=F) + 
  labs(x="% unfolded protein", y=paste("% unfolded \u03B2", "2", sep=""))+
  scale_x_continuous(breaks=c(0, 25, 50, 75,  100), labels=c(0, 25, 50, 75, 100), limits=c(0, 102)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 50, 100)) + 
  scale_color_manual(values=col_a) + 
  geom_line(data=data.frame(pos=c(1:99), avg=as.vector(colMeans(unfold_b2))), inherit.aes=F, aes(x=pos, y=avg), col="red", size=1.3) + 
  #geom_line(data=data.frame(pos=c(1:99), avg=as.vector(colMeans(refold_b2))), inherit.aes=F, aes(x=pos, y=avg), col="blue", size=1.3) + 
  theme_classic() +
  theme(
    text = element_text(size=18)
  )

p.re_b2 <- ggplot(REdf_b2, aes(x=index)) + 
  geom_line(aes(y=folded, color=factor(traj)), size=0.3, show.legend=F) + 
  labs(x="% folded protein", y=paste("% folded \u03B2", "2", sep=""))+
  scale_x_continuous(breaks=c(0, 25, 50, 75,  100), labels=c(0, 25, 50, 75, 100), limits=c(0, 102)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 50, 100)) + 
  scale_color_manual(values=col_b) + 
  #geom_line(data=data.frame(pos=c(1:99), avg=as.vector(colMeans(unfold_b2))), inherit.aes=F, aes(x=pos, y=avg), col="red", size=1.3) + 
  geom_line(data=data.frame(pos=c(1:99), avg=as.vector(colMeans(refold_b2))), inherit.aes=F, aes(x=pos, y=avg), col="blue", size=1.3) + 
  theme_classic() +
  theme(
    text = element_text(size=18)
  )

p.un_b6 <- ggplot(UNdf_b6, aes(x=index)) + 
  geom_line(aes(y=folded, color=factor(traj)), size=0.3, show.legend=F) + 
  labs(x="% unfolded protein", y=paste("% unfolded \u03B2", "6", sep=""))+
  scale_x_continuous(breaks=c(0, 25, 50, 75,  100), labels=c(0, 25, 50, 75, 100), limits=c(0, 102)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 50, 100)) + 
  scale_color_manual(values=col_a) + 
  geom_line(data=data.frame(pos=c(1:99), avg=as.vector(colMeans(unfold_b6))), inherit.aes=F, aes(x=pos, y=avg), col="red", size=1.3) + 
  #geom_line(data=data.frame(pos=c(1:99), avg=as.vector(colMeans(refold_b2))), inherit.aes=F, aes(x=pos, y=avg), col="blue", size=1.3) + 
  theme_classic() +
  theme(
    text = element_text(size=18)
  )

p.re_b6 <- ggplot(REdf_b6, aes(x=index)) + 
  geom_line(aes(y=folded, color=factor(traj)), size=0.3, show.legend=F) + 
  labs(x="% folded protein", y=paste("% folded \u03B2", "6", sep=""))+
  scale_x_continuous(breaks=c(0, 25, 50, 75,  100), labels=c(0, 25, 50, 75, 100), limits=c(0, 102)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 50, 100)) + 
  scale_color_manual(values=col_b) + 
  #geom_line(data=data.frame(pos=c(1:99), avg=as.vector(colMeans(unfold_b2))), inherit.aes=F, aes(x=pos, y=avg), col="red", size=1.3) + 
  geom_line(data=data.frame(pos=c(1:99), avg=as.vector(colMeans(refold_b6))), inherit.aes=F, aes(x=pos, y=avg), col="blue", size=1.3) + 
  theme_classic() +
  theme(
    text = element_text(size=18)
  )


plot_grid(p.un_b2, p.re_b2, p.un_b6, p.re_b6, labels ="", ncol = 2, align = 'h')

dev.off()






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL S5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#5. N-terminal degron depletion


degron <- as.data.frame(read.table("data/processed/degron_probs.txt", header=T))


svg(file = "figures/SI/degron_Nter.svg", height = 4, width = 4)

R_Nd70 <- round(cor(degron$Nunfold, degron$d70_Np), 3)
ggplot(degron, aes(x=Nunfold, y=d70_Np)) + 
  geom_smooth(method=lm, color="orange", fill="orange", se=F, alpha=0.3) + 
  geom_point(size=2) + 
  geom_text(x=0.7, y=15, label=paste("R =", R_Nd70, sep=' '), size=7, stat="unique") + 
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c(0, 25, 50, 75, 100)) + 
  labs(x="% N-unfolding", y="N-ter degron score") +
  theme_classic() + 
  theme(
    text = element_text(size=20)
  )

dev.off()




###
# degron depletion by class

degron$class <- ifelse(degron$Cunfold > 0.5, "C", "N")


svg(file = "figures/SI/degron_control.svg", height = 6, width = 3)

p.Cter <- ggplot(degron, aes(x=class)) + 
  geom_boxplot(aes(y=d70_C, fill=class)) + 
  labs(x="", y="Degron score", title="C-ter") + 
  scale_fill_manual(values=c("#b735b7", "#ffd17f")) + 
  theme_classic() + 
  theme(
    text = element_text(size=16), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(), 
    legend.position = 'none'
  )

p.middle <- ggplot(degron, aes(x=class)) + 
  geom_boxplot(aes(y=d70_middle, fill=class)) + 
  labs(x="", y="Degron score", title="middle") + 
  scale_fill_manual(values=c("#b735b7", "#ffd17f")) + 
  theme_classic() + 
  theme(
    text = element_text(size=16), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(), 
    legend.position = 'none'
  )

p.Nter <- ggplot(degron, aes(x=class)) + 
  geom_boxplot(aes(y=d70_Np, fill=class)) + 
  labs(x="", y="Degron score", title="N-ter") + 
  scale_fill_manual(values=c("#b735b7", "#ffd17f")) + 
  theme_classic() + 
  theme(
    text = element_text(size=16), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(), 
    legend.position = 'none'
  )


plot_grid(p.Cter,  p.middle, p.Nter, labels ="", ncol = 1, align = 'v')

dev.off()


#wilcox.test(degron$d70_C[degron$class=="C"], degron$d70_C[degron$class=="N"])
#wilcox.test(degron$d70_middle[degron$class=="C"], degron$d70_middle[degron$class=="N"])
#wilcox.test(degron$d70_Np[degron$class=="C"], degron$d70_Np[degron$class=="N"])


svg(file = "figures/SI/degron_control2.svg", height = 4, width = 6)

plot_grid(p.Cter,  p.middle, p.Nter, labels ="", ncol = 3, align = 'h')

dev.off()
