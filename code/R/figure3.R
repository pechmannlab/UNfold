# UNfold (2024)
# (c) S Pechmann
# Figure 3

library(ggplot2)
library(dplyr)
library(igraph)
library(networkD3)
library(htmlwidgets)
library(cowplot)
library(reshape2)
library(ape)

setwd("~/UNfold/")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pymol structure


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

allDF <- as.data.frame(read.table("code/tmp.all", header=T))               #"../data/processed/allDF.txt"
bp <- allDF$YNL093W_07
clusters <- as.data.frame(read.table("data/processed/cluster_members.txt", header=T))

d <- as.matrix(read.table("code/YNL093W_07_distances.txt"))
qf <- scan("code/YNL093W_07_qf.txt")


sdat <- matrix(NA, nrow=length(qf), ncol=7)
for (i in 1:(length(qf)-20) ){
  sdat[i,1] <- mean(qf[i:(i+20)])
  sdat[i,2] <- mean(d[i:(i+20), 316])
  sdat[i,3] <- mean(d[i:(i+20), 476])
  sdat[i,4] <- mean(d[i:(i+20), 477])
  sdat[i,5] <- mean(d[i:(i+20), 327])
  sdat[i,6] <- mean(d[i:(i+20), 486])
  sdat[i,7] <- mean(d[i:(i+20), 487])
}
colnames(sdat) <- c('qfs', 'c41.1', 'c41.2', 'c41.3', 'c64.1', 'c64.2', 'c64.3')
sdat <- as.data.frame(sdat)
sdat$pos <- c(1:length(qf))
sdat$qf <- qf


svg("figures/figure3/B_qf_dist.svg", height=4, width=4)

p.qf <- ggplot(sdat, aes(x=pos)) + 
  geom_vline(xintercept=3053, color="red", size=0.7, linetype="dashed", alpha=0.5) + 
  geom_vline(xintercept=3352, color="darkblue", size=0.7, linetype="dashed", alpha=0.5) + 
  geom_line(aes(y=qf), alpha=0.2) + 
  geom_line(aes(y=qfs)) + 
  labs(x="", y="Qf") + 
  #scale_x_continuous(limits=c(1000, 6000)) + 
  theme_classic() + 
  theme(
    text = element_text(size=20),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )

p.d <- ggplot(sdat, aes(x=pos)) + 
  geom_vline(xintercept=3053, color="red", size=0.7, linetype="dashed", alpha=0.5) + 
  geom_vline(xintercept=3352, color="darkblue", size=0.7, linetype="dashed", alpha=0.5) + 
  geom_line(aes(y=c41.1), alpha=0.5, col="red") + 
  geom_line(aes(y=c41.2), alpha=0.5, col="red") + 
  geom_line(aes(y=c41.3), alpha=0.5, col="red") + 
  geom_line(aes(y=c64.1), alpha=0.5, col="darkblue") + 
  geom_line(aes(y=c64.2), alpha=0.5, col="darkblue") + 
  geom_line(aes(y=c64.3), alpha=0.5, col="darkblue") + 

  labs(x="Simulation time", y="Distance") + 
  theme_classic()+
  theme(
    text = element_text(size=20)
  )

plot_grid(p.qf, p.d, labels ="", ncol = 1, align = 'v')

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source('~/sw/R/ggsankey-main/R/utils-pipe.R')
source('~/sw/R/ggsankey-main/R/sankey.R')
#library(ggsankey)


svg("figures/figure3/C_sankey.svg", height=4, width=7)

d <- as.data.frame(read.table("data/processed/assign_bpinit.txt", header=T))
d <- d[,c(2,5)]

order <- rev(c("YDL192W", "YOR094W", "YNL093W", "YMR138W", "YBR164C", "YNL090W", "YHR022C", "YLR293C", "YOR089C", "YOR185C", "YFL005W", "YCR027C", "YKR014C", "YKL154W", "YBR264C", "YFL038C"))


df <- d %>% make_long(class, protein)

df$node <- factor(df$node, levels=c('a2/b2', 'other', 'a1', 'a3', 'a4/b5', 'a5', order ))
df$next_node <- factor(df$next_node, levels=c('a2/b2', 'other', 'a1', 'a3', 'a4/b5', 'a5', order))

ggplot(df, aes(x = x,
               next_x = next_x,
               node = node,
               next_node = next_node,
               fill = factor(node), 
               label=node)) +
  geom_sankey(flow.alpha=0.5, node.color="white", node.fill="grey50", show.legend=F) + 
  #geom_sankey_label() + 
  geom_sankey_text(hjust=0) + 
  scale_fill_manual(values=c("blue", "orange", "red", "darkgreen", "purple", "#777777", rep("grey50", 16))) + 
  labs(x="", y="") + 
  theme_classic() +
  theme(
    text = element_text(size=20),
    #axis.line.x = element_blank(), 
    #axis.ticks.x = element_blank(),
    axis.line = element_blank(), 
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = 'none'
  )

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bp_YNL093W_07 <- allDF$YNL093W_07
bp_YDL192W_51 <- allDF$YDL192W_51
bp_YCR027C_08 <- allDF$YCR027C_08

qf_YNL093W_07 <- scan("code/YNL093W_07_qf.txt")
qf_YDL192W_51 <- scan("code/YDL192W_51_qf.txt")
qf_YCR027C_08 <- scan("code/YCR027C_08_qf.txt")

qfDF <- data.frame(YNL093W_07=qf_YNL093W_07, YDL192W_51=qf_YDL192W_51, YCR027C_08=qf_YCR027C_08)
qfDF$pos <- c(1:length(qf_YNL093W_07))


svg("figures/figure3/D_qfex.svg", height=3, width=12)

p1 <- ggplot(qfDF, aes(x=pos)) + 
  geom_vline(xintercept=bp_YNL093W_07, color="orange", size=0.2, alpha=0.5) + 
  geom_line(aes(y=YNL093W_07), col="#222222") + 
  scale_x_continuous(limits=c(2000, 4000)) + 
  labs(x="Simulation time", y="Qf", title="YNL093W") +
  theme_classic() + 
  theme(
    text = element_text(size=16)
  )

p2 <- ggplot(qfDF, aes(x=pos)) + 
  geom_vline(xintercept=bp_YDL192W_51, color="purple", size=0.2, alpha=0.5) + 
  geom_line(aes(y=YDL192W_51), col="#222222") + 
  scale_x_continuous(limits=c(3000, 5500)) + 
  labs(x="Simulation time", y="Qf", title="YDL192W") +
  theme_classic() + 
  theme(
    text = element_text(size=16)
  )

p3 <- ggplot(qfDF, aes(x=pos)) + 
  geom_vline(xintercept=bp_YCR027C_08, color="red", size=0.2, alpha=0.5) + 
  geom_line(aes(y=YCR027C_08), col="#222222") + 
  scale_x_continuous(limits=c(2000, 4000)) + 
  labs(x="Simulation time", y="Qf", title="YCR027C") +
  theme_classic() + 
  theme(
    text = element_text(size=16)
  )


plot_grid(p1, p2, p3, labels ="", ncol = 3, align = 'h')

dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ssunf <- as.data.frame(read.table("data/processed/ssunfold_byprotein.txt", header=T))

a1 <- ssunf[ssunf$ss=="a1",]
a2 <- ssunf[ssunf$ss=="a2",]
a3 <- ssunf[ssunf$ss=="a3",]
a4 <- ssunf[ssunf$ss=="a4",]
a5 <- ssunf[ssunf$ss=="a5",]
b1 <- ssunf[ssunf$ss=="b1",]
b2 <- ssunf[ssunf$ss=="b2",]
b3 <- ssunf[ssunf$ss=="b3",]
b4 <- ssunf[ssunf$ss=="b4",]
b5 <- ssunf[ssunf$ss=="b5",]
b6 <- ssunf[ssunf$ss=="b6",]
nn <- ssunf[ssunf$ss=="none",]

colfunc_a <- colorRampPalette(c("red", "darkred"))
col_a <- colfunc_a(16)

colfunc_b <- colorRampPalette(c("#e527f9", "#81128d"))
col_b <- colfunc_b(16)


p.b1 <- ggplot(b1, aes(x=pos)) + 
  geom_line(aes(y=val, color=pdb), size=0.5, show.legend=F) + 
  labs(x="", y="", title="b1")+
  scale_x_continuous(breaks=c(0, 50, 100), labels=c(0, 50, 100), limits=c(0, 102)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 50, 100)) + 
  scale_color_manual(values=col_b) + 
  theme_classic() + 
  theme(
    text = element_text(size=16),
    axis.text.x = element_blank()
  )

p.a1 <- ggplot(a1, aes(x=pos)) + 
  geom_line(aes(y=val, color=pdb), size=0.5, show.legend=F) + 
  labs(x="", y="", title="a1")+
  scale_x_continuous(breaks=c(0, 50, 100), labels=c(0, 50, 100), limits=c(0, 102)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 50, 100)) + 
  scale_color_manual(values=col_a) + 
  theme_classic() +
  theme(
    text = element_text(size=16),
    axis.text = element_blank()
  )

p.b2 <- ggplot(b2, aes(x=pos)) + 
  geom_line(aes(y=val, color=pdb), size=0.5, show.legend=F) + 
  labs(x="", y="", title="b2")+
  scale_x_continuous(breaks=c(0, 50, 100), labels=c(0, 50, 100), limits=c(0, 102)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 50, 100)) + 
  scale_color_manual(values=col_b) + 
  theme_classic() +
  theme(
    text = element_text(size=16),
    axis.text = element_blank()
  )

p.b3 <- ggplot(b3, aes(x=pos)) + 
  geom_line(aes(y=val, color=pdb), size=0.5, show.legend=F) + 
  labs(x="", y="", title="b3")+
  scale_x_continuous(breaks=c(0, 50, 100), labels=c(0, 50, 100), limits=c(0, 102)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 50, 100)) + 
  scale_color_manual(values=col_b) + 
  theme_classic() +
  theme(
    text = element_text(size=16),
    axis.text = element_blank()
  )

p.a2 <- ggplot(a2, aes(x=pos)) + 
  geom_line(aes(y=val, color=pdb), size=0.5, show.legend=F) + 
  labs(x="", y="", title="a2")+
  scale_x_continuous(breaks=c(0, 50, 100), labels=c(0, 50, 100), limits=c(0, 102)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 50, 100)) + 
  scale_color_manual(values=col_a) + 
  theme_classic() +
  theme(
    text = element_text(size=16),
    axis.text.x = element_blank()
  )

p.b4 <- ggplot(b4, aes(x=pos)) + 
  geom_line(aes(y=val, color=pdb), size=0.5, show.legend=F) + 
  labs(x="", y="", title="b4")+
  scale_x_continuous(breaks=c(0, 50, 100), labels=c(0, 50, 100), limits=c(0, 102)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 50, 100)) + 
  scale_color_manual(values=col_b) + 
  theme_classic() +
  theme(
    text = element_text(size=16),
    axis.text = element_blank()
  )


p.a3 <- ggplot(a3, aes(x=pos)) + 
  geom_line(aes(y=val, color=pdb), size=0.5, show.legend=F) + 
  labs(x="", y="", title="a3")+
  scale_x_continuous(breaks=c(0, 50, 100), labels=c(0, 50, 100), limits=c(0, 102)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 50, 100)) + 
  scale_color_manual(values=col_a) + 
  theme_classic() +
  theme(
    text = element_text(size=16),
    axis.text = element_blank()
  )


p.b5 <- ggplot(b5, aes(x=pos)) + 
  geom_line(aes(y=val, color=pdb), size=0.5, show.legend=F) + 
  labs(x="", y="", title="b5")+
  scale_x_continuous(breaks=c(0, 50, 100), labels=c(0, 50, 100), limits=c(0, 102)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 50, 100)) + 
  scale_color_manual(values=col_b) + 
  theme_classic() +
  theme(
    text = element_text(size=16),
    axis.text = element_blank()
  )


p.a4 <- ggplot(a4, aes(x=pos)) + 
  geom_line(aes(y=val, color=pdb), size=0.5, show.legend=F) + 
  labs(x="", y="", title="a4") + 
  scale_x_continuous(breaks=c(0, 50, 100), labels=c(0, 50, 100), limits=c(0, 102)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 50, 100)) + 
  scale_color_manual(values=col_a) + 
  theme_classic() + 
  theme(
    text = element_text(size=16)
  )

p.b6 <- ggplot(b6, aes(x=pos)) + 
  geom_line(aes(y=val, color=pdb), size=0.5, show.legend=F) + 
  labs(x="", y="", title="b6") +
  scale_x_continuous(breaks=c(0, 50, 100), labels=c(0, 50, 100), limits=c(0, 102)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 50, 100)) + 
  scale_color_manual(values=col_b) + 
  theme_classic() + 
  theme(
    text = element_text(size=16),
    axis.text.y = element_blank()
  )

p.a5 <- ggplot(a5, aes(x=pos)) + 
  geom_line(aes(y=val, color=pdb), size=0.5, show.legend=F) + 
  labs(x="", y="", title="a5")+
  scale_x_continuous(breaks=c(0, 50, 100), labels=c(0, 50, 100), limits=c(0, 102)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 50, 100)) + 
  scale_color_manual(values=col_a) + 
  theme_classic() + 
  theme(
    text = element_text(size=16),
    axis.text.y = element_blank()
  )


p.nn <- ggplot(nn, aes(x=pos)) + 
  geom_line(aes(y=val, color=pdb), size=0.5, show.legend=F) + 
  scale_x_continuous(breaks=c(0, 50, 100), labels=c(0, 50, 100), limits=c(0, 102)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 50, 100)) + 
  scale_color_manual(values=rep("grey50", 16)) + 
  labs(x="", y="", title="other")+
  theme_classic() + 
  theme(
    text = element_text(size=16),
    axis.text.y = element_blank()
  )

#xlab: % unfolded protein
#ylab: % unfolded ss

svg("figures/figure3/E_ssunf.svg", height=6, width=8)
plot_grid(p.b1, p.a1, p.b2, p.b3, p.a2, p.b4, p.a3, p.b5, p.a4, p.b6, p.a5, p.nn, nrow=3, ncol = 4, align = 'hv')
dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL F ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

auc <- as.data.frame(read.table("data/processed/ssunfold_bystart_auc.txt", header=T))

b2 <- auc[auc$ss=="b2",]
b3 <- auc[auc$ss=="b3",]
b5 <- auc[auc$ss=="b5",]
b6 <- auc[auc$ss=="b6",]
bdf <- rbind(b2, b3, b5, b6)
bdf$start_ss <- factor(bdf$start_ss, levels=c("a1", "a2/b2", "a3", "a4/b5", "a5", "other"))



p.b <- ggplot(bdf, aes(x=start_ss)) + 
  geom_col(aes(y=auc, fill=ss), position=position_dodge2() ) + 
  geom_errorbar(aes(ymin=auc-auc_sd, ymax=auc+auc_sd), position=position_dodge2(width = 0.2, padding = 0.8), col="grey50", alpha=0.5) + 
  labs(y="Unfolding AUC", x="Initial unfolding") + 
  scale_y_continuous(limits=c(0, 100)) + 
  scale_fill_manual(values = c("#E527F9", "#E527F9", "#81128D", "#81128D" )) + 
  coord_flip() + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    axis.line.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    legend.title = element_blank(),
    legend.position = c(0.89, 0.83)
  )



a1 <- auc[auc$ss=="a1",]
a3 <- auc[auc$ss=="a3",]
a4 <- auc[auc$ss=="a4",]
a5 <- auc[auc$ss=="a5",]
adf <- rbind(a1, a3, a4, a5)
adf$start_ss <- factor(adf$start_ss, levels=c("a1", "a2/b2", "a3", "a4/b5", "a5", "other"))
adf$col <- factor(c(1:nrow(adf)))
alpha_sel <- c(0.8,0.8,1,0.8,0.8,0.8,  0.8,0.8,0.8,1,0.8,0.8,  0.8,0.8,0.8,0.8,0.8,1,  0.8,0.8,0.8,0.8,1,0.8)

p.a <- ggplot(adf, aes(x=start_ss)) + 
  geom_col(aes(y=auc, fill=ss, alpha=alpha_sel), position=position_dodge2() ) + 
  geom_errorbar(aes(ymin=auc-auc_sd, ymax=auc+auc_sd), position=position_dodge2(width = 0.2, padding = 0.8), col="grey50", alpha=0.5) + 
  labs(y="Unfolding AUC", x="Initial unfolding") + 
  scale_y_continuous(limits=c(0, 100)) + 
  scale_fill_manual(values = c("#FF0000", "#E00000", "#B80000","#900000")) + 
  coord_flip() + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    axis.line.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    legend.title = element_blank(),
    legend.position = c(0.92, 0.7)
  )

svg("figures/figure3/F_ssauc.svg", height=4, width=6)

plot_grid(  p.b, p.a, ncol = 2, align = 'h')

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL G ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bpinit <- as.data.frame(read.table("data/processed/assign_bpinit.txt", header=T))

pdblist <- unique(bpinit$protein)
df <- data.frame(pdb=c(), a2=c(), a5=c())
for (i in 1:length(pdblist)){
  current_pdb <- as.vector(pdblist[i])
  current_data <- bpinit[bpinit$protein==current_pdb,]
  current_class <- current_data$class
  current_a2 <- sum(current_class == "a2/b2")/length(current_class)
  current_a5 <- sum(current_class == "a5")/length(current_class)
  
  df <- rbind(df, data.frame(pdb=current_pdb, a2=current_a2, a5=current_a5))
}

preclass <- c()
for (i in 1:nrow(df)){
  if (df$a2[i] > 0.5){preclass <- c(preclass, "a2")}
  else if (df$a5[i] > 0.5){preclass <- c(preclass, "a5")}
  else {preclass <- c(preclass, "other")}
}
df$class <- preclass

selclust <- as.data.frame(read.table("code/tmp.selclust", header=T))



sel_46 <- selclust[selclust$cluster==46, ]
c46 <- cbind(df, sel_46[,c(2,3,4)])



svg("figures/figure3/G_a2b2_dist.svg", height=4, width=4)

ggplot(c46, aes(x=a2, y=mean)) + 
  geom_smooth(method=lm, alpha=0.1, size=1, color="#222222", fill="orange", se=T) + 
  geom_point(aes(color=class, shape=class), size=3) + 
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std),  col="grey50", alpha=0.5) +     # position=position_dodge2(width = 0.2, padding = 0.8),
  geom_text(aes(label=pdb), nudge_x=0.055, nudge_y=0.07, size=4) + 
  scale_x_continuous(breaks=c(0, 0.3, 0.6, 0.9), labels=c(0, 30, 60, 90), limits=c(0, 0.95)) +
  labs(x="% class a2/b2", y="Mean distance") + 
  scale_color_manual(values=c("orange", "purple", "#777777")) + 
  theme_classic() + 
  theme(
    text = element_text(size=16), 
    legend.title = element_blank(), 
    legend.position = c(0.15, 0.88)
    
  )

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL H ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

contss <- as.data.frame(read.table("data/processed/contacts_aligned_ss.txt", header=T))
aln <- read.table("data/processed/c3718_aln.txt", sep='\t') 


# !! df from panel g !! 
# !! alignment as a different order !!
Cter <- c()
for (i in 1:nrow(df)){
  current_pdb <- as.vector(df$pdb)[i]
  current_idx <- which(aln[,1]==current_pdb)
  current_cter <- aln[current_idx, 392:524]       # 1 column names + 390 aligned
  current_len <- sum(current_cter!='-')
  print(c(current_pdb, current_idx, current_len))
  Cter <- c(Cter, current_len)
}
df$Cter <- Cter
df$conta5 <- contss$a5




svg("figures/figure3/H_a5_cont.svg", height=4, width=5)

p1 <- ggplot(df, aes(x=a5, y=conta5)) + 
  geom_smooth(method=lm, alpha=0.1, size=1, color="#222222", fill="purple", se=T) + 
  geom_point(aes(color=class, shape=class), size=3) + 
  #geom_text(aes(label=pdb), nudge_x=0.055, nudge_y=0.07, size=4) + 
  scale_x_continuous(breaks=c(0, 0.3, 0.6, 0.9), labels=c(0, 30, 60, 90), limits=c(0, 0.9)) +
  labs(x="% class a5", y="Contacts (a5)") + 
  scale_color_manual(values=c("orange", "purple", "#777777")) + 
  theme_classic() + 
  theme(
    text = element_text(size=18), 
    legend.title = element_blank(), 
    legend.position = 'none'
  )


p2 <- ggplot(df, aes(x=a5, y=Cter)) + 
  geom_smooth(method=lm, alpha=0.1, size=1, color="#222222", fill="purple", se=T) + 
  geom_point(aes(color=class, shape=class), size=3) + 
  #geom_text(aes(label=pdb), nudge_x=0.055, nudge_y=0.07, size=4) + 
  scale_x_continuous(breaks=c(0, 0.3, 0.6, 0.9), labels=c(0, 30, 60, 90), limits=c(0, 0.9)) +
  labs(x="% class a5", y="C-ter (aa)") + 
  scale_color_manual(values=c("orange", "purple", "#777777")) + 
  theme_classic() + 
  theme(
    text = element_text(size=18), 
    legend.title = element_blank(), 
    legend.position = 'none'
  )

plot_grid(  p1, p2, ncol = 1, align = 'v')


dev.off()

