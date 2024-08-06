# UNfold (2024)
# (c) S Pechmann
# Figure 4

library(ggplot2)
library(ggridges)
library(viridis)
library(cowplot)


setwd("~/UNfold/")

# COMMENT: some figure panels use data from earlier panels, so reload all data 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pymol folded/unfolded aggregation surface example



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

surf <- as.data.frame(read.table("data/processed/surfaceAP_trajectory.txt", header=T))
qf <- as.data.frame(read.table("data/processed/qf_trajectory.txt", header=T))

# plotting first two trajs: "YOR089C_43" "YOR089C_54"
data_qf <- qf[,c(1,2)]
data_qf$position <- c(1:nrow(data_qf))

agg1_mean <- rep(NA, nrow(data_qf))
agg1_sd <- rep(NA, nrow(data_qf))
agg2_mean <- rep(NA, nrow(data_qf))
agg2_sd <- rep(NA, nrow(data_qf))

w <- 100
for (i in 1:(nrow(data_qf) - w)  ){
  agg1_mean[i] <- mean(surf[i:(i+w),1])
  agg1_sd[i] <- sd(surf[i:(i+w),1])
  agg2_mean[i] <- mean(surf[i:(i+w),2])
  agg2_sd[i] <- sd(surf[i:(i+w),2])
}
mmax <- max(c(max(surf[1]), max(surf[,2])))
data_qf$agg1_mean <- (agg1_mean / mmax) 
data_qf$agg1_sd <- (agg1_sd / mmax) 
data_qf$agg2_mean <- (agg2_mean / mmax) 
data_qf$agg2_sd <- (agg2_sd / mmax) 




p1 <- ggplot(data_qf, aes(x=position)) +
  geom_line(aes(y=YOR089C_43), color="grey20", size=0.8, alpha=0.5) +
  geom_ribbon(aes(ymax=agg1_mean+agg1_sd, ymin=agg1_mean-agg1_sd), fill="orange", alpha=0.2) +
  geom_line(aes(y=agg1_mean), color="orange", size=0.8, alpha=0.5) +
  labs(x="Simulation time", y="Qf") + 
  scale_y_continuous(sec.axis = sec_axis( ~ . * mmax, name="Aggregation score")) +
  theme_classic() + 
  theme(
     text = element_text(size=16)
  )

  
p2 <- ggplot(data_qf, aes(x=position)) +
    geom_line(aes(y=YOR089C_54), color="grey70", size=0.8, alpha=0.5) +
    geom_ribbon(aes(ymax=agg2_mean+agg2_sd, ymin=agg2_mean-agg2_sd), fill="red", alpha=0.2) +
    geom_line(aes(y=agg2_mean), color="red", size=0.8, alpha=0.5) +
    labs(x="Simulation time", y="Qf") + 
    scale_y_continuous(sec.axis = sec_axis( ~ . * mmax, name="Aggregation score")) +
    theme_classic() + 
    theme(
    text = element_text(size=16)
  )
  
  

  
  
svg(file = "figures/figure4/B_Qfagg.svg", height = 4, width = 5)

plot_grid(p1, p2, labels ="", ncol = 1, align = 'v')
  
dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

surf <- as.data.frame(read.table("data/processed/surfaceAP_trajectory.txt", header=T))
qf <- as.data.frame(read.table("data/processed/qf_trajectory.txt", header=T))
rg <- as.data.frame(read.table("data/processed/rg_trajectory.txt", header=T))

genes <- read.table("data/processed/gene_names.txt", header=T)


traj_names <- names(surf)
orf_names <- substr(traj_names, 1, 7)
gene_names <- rep(NA, length(orf_names))
for (i in 1:length(gene_names)){gene_names[i] <- as.character(genes$name[genes$orf==orf_names[i]])}

result_folded <- data.frame(orf=c(), ap=c())
result_unfolded <- data.frame(orf=c(), ap=c())
for (i in 1:length(traj_names)){                   #
    current_qf<- qf[,i]
    sel_folded <- current_qf > 0.7
    sel_unfolded <- current_qf < 0.3 
    current_data <- surf[,i]
    
    data_folded <- current_data[sel_folded]  
    #data_folded <- sample(data_folded, floor(length(data_folded)/10) )
                          
    data_unfolded <- current_data[sel_unfolded]   
    #data_unfolded <- sample(data_unfolded, floor(length(data_unfolded)/10) )
    
    current_df_folded <- data.frame(orf=rep(factor(gene_names[i]), length(data_folded)), ap=data_folded )
    current_df_unfolded <- data.frame(orf=rep(factor(gene_names[i]), length(data_unfolded)), ap=data_unfolded )
    result_folded <- rbind(result_folded, current_df_folded)
    result_unfolded <- rbind(result_unfolded, current_df_unfolded)
  } 



#list_orf <- unique(orf_names)
list_orf <- unique(gene_names)
sortdf <- data.frame(orf=c(), mean=c())
for (i in 1:length(list_orf)){
  current_orf <- list_orf[i]
 
  current_data <- result_unfolded[result_unfolded$orf==current_orf,]
  current_mean <- mean(current_data$ap)
  current_df <- data.frame(orf=current_orf, mean=current_mean)
  #current_df <- data.frame(orf=current_gn, mean=current_mean)
  sortdf <- rbind(sortdf, current_df)
}
sortdf <- sortdf[sort(sortdf$mean, index.return=T, decreasing=T)$ix, ]

result_folded$orf <- factor(result_folded$orf, levels=sortdf$orf)
result_unfolded$orf <- factor(result_unfolded$orf, levels=sortdf$orf)



p.folded <- ggplot(result_folded, aes(x=ap, y=orf, fill= ..x..) ) +         
  geom_density_ridges_gradient( scale=3, rel_min_height=0.01 ) +  
  scale_fill_viridis(option="G", begin=0.4, end=1) + 
  scale_x_continuous(limits=c(0, 15000)) + 
  labs(x="Aggregation score", y="", subtitle="folded") +
  theme_classic() + 
  theme(
    text = element_text(size=16),
    axis.line = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.text.y = element_text(size=7),
    legend.position = 'none'
  )
  


p.unfolded <- ggplot(result_unfolded, aes(x=ap, y=orf, fill= ..x..) ) +         
  geom_density_ridges_gradient( scale=3, rel_min_height=0.01 ) +  
  scale_fill_viridis(option="F", begin=0., end=1) + 
  scale_x_continuous(limits=c(0, 15000)) + 
  labs(x="Aggregation score", y="", subtitle="unfolded") +
  theme_classic() + 
  theme(
    text = element_text(size=16),
    axis.line = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.text.y = element_text(size=7),
    legend.position = 'none'
  )





svg(file = "figures/figure4/C_aggsurf.svg", height=4, width=4)

plot_grid(p.folded, p.unfolded, labels ="", ncol = 1, align = 'v')

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

imdf <- as.data.frame(read.table("data/processed/unfolding_intermediates.txt", header=T))
imdf$pdb <- substr(imdf$traj,1, 7)

genes <- read.table("data/processed/gene_names.txt", header=T)
gn <- rep(NA, length(imdf$pdb))
for (i in 1:length(gn)){gn[i] <- as.character(genes$name[genes$orf==orf_names[i]])}
imdf$name <- gn

pdblist <- unique(imdf$name)
m <- rep(NA, length(pdblist))
for (i in 1:length(pdblist)){
  current_pdb <- as.vector(pdblist)[i]
  current_dat <- imdf[imdf$name==current_pdb,]
  current_mean <- mean(current_dat$mindist)
  m[i] <- current_mean
}
order_pdb <- pdblist[sort(m, index.return=T)$ix]
imdf$name <- factor(imdf$name, levels=rev(order_pdb))


svg(file = "figures/figure4/D_unfoldint.svg", height=4, width=5.5)

p1 <- ggplot(imdf, aes(x=end-start, y=name) ) +         
  geom_boxplot(aes(fill=name)) +
  #scale_fill_viridis_d(option="C", begin=0.3, end=1) + 
  scale_fill_viridis_d(option="D", begin=0.4, end=0.95, direction=-1) + 
  scale_x_continuous(limits=c(0, 500), breaks=c(0, 250, 500)) + 
  labs(x="Unfolding [sim.time]", y="") +
  theme_classic() + 
  theme(
    text = element_text(size=14),
    axis.line.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.text.y = element_text(size=12),
    legend.position = 'none'
  )

p2 <- ggplot(imdf, aes(x=mindist, y=name) ) +         
  geom_boxplot(aes(fill=name)) +
  scale_fill_viridis_d(option="D", begin=0.4, end=0.95, direction=-1) + 
  scale_x_continuous(limits=c(0, 500), breaks=c(0, 250, 500)) + 
  labs(x="Intermediate [sim.time]", y="") +
  theme_classic() + 
  theme(
    text = element_text(size=14),
    axis.line.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.text.y = element_text(size=12),
    legend.position = 'none'
  )


plot_grid(p1, p2, labels ="", ncol = 2, align = 'h')

dev.off()








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ufint <- as.data.frame(read.table("data/processed/unfolding_intermediates.txt", header=T))
ufint_sel <- ufint[ufint$traj=='YDL192W_07',]

qf <- scan("data/processed/qf_YDL192W_07.txt")
qf20 <- rep(NA, length(qf))
for (i in 1:(length(qf)-20)){qf20[i] <- mean(qf[i:(i+20)])}
qfdf <- data.frame(pos=c(1:length(qf)), qf=qf, qf20=qf20)

polydf <- data.frame(x=c(0, 10000, 10000, 0, 0, 10000, 10000, 0), y=c(0.1, 0.1, 0.3, 0.3, 0.8, 0.8, 1, 1), class=c(rep("unfolded", 4), rep("folded", 4)))


svg(file = "figures/figure4/E_refold.svg", height=4, width=6)

ggplot(qfdf, aes(x=pos)) +
  geom_polygon(data=polydf, inherit.aes=F, aes(x=x, y=y, fill=class), show.legend= FALSE, alpha=0.1) + 
  geom_vline(xintercept=ufint_sel$start, color="#222222", size=0.5, linetype="dashed") + 
  geom_vline(xintercept=ufint_sel$end, color="#222222", size=0.5, linetype="dashed") + 
  geom_vline(xintercept=ufint_sel$refold, color="#222222", size=0.5, linetype="dashed") + 
  geom_line(aes(y=qf), col="#777777", alpha=0.7) + 
  geom_line(aes(y=qf20), col="purple") + 
  labs(x="Simulation time", y="Qf") + 
  scale_y_continuous(limits=c(0.1, 1)) + 
  scale_fill_manual(values=c("#c02a64", "#4675a7")) +
  theme_classic() +
  theme(
    text = element_text(size=20)
  )

dev.off()

data_sel <- imdf[imdf$pdb=="YDL192W",]
c_total <- nrow(data_sel)
c_rf <- nrow(data_sel[data_sel$refold > 0, ])
pct <- c_rf / c_total
# 6/65 is refolding events of YDL192W

rfdf <- data.frame(pdb=order_pdb, rfpct=c(rep(0.001, 15), pct), pos=c(1:16))


svg(file = "figures/figure4/E_refoldpct.svg", height=3, width=6)

ggplot(rfdf, aes(x=pos, y=rfpct)) + 
  geom_col() + 
  scale_x_continuous(breaks=rfdf$pos, labels=rfdf$pdb) + 
  scale_y_continuous(breaks=c(0, 0.05, 0.1), labels=c(0, 5, 10), limits=c(0, 0.1)) + 
  labs(x="", y="% refolding") + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), 
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank()
  )

dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL F ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


qf <- as.data.frame(read.table("data/processed/qf_trajectory.txt", header=T))
rg <- as.data.frame(read.table("data/processed/rg_trajectory.txt", header=T))


imdf <- as.data.frame(read.table("data/processed/unfolding_intermediates.txt", header=T))
imdf$pdb <- substr(imdf$traj,1, 7)






traj_names <- names(rg)
orf_names <- substr(traj_names, 1, 7)
result_folded <- data.frame(orf=c(), rg=c())
result_unfolded <- data.frame(orf=c(), rg=c())

for (i in 1:length(traj_names)){                   #
  current_qf<- qf[,i]
  sel_folded <- current_qf > 0.7 
  sel_unfolded <- current_qf < 0.3
  current_data <- rg[,i]
  
  data_folded <- current_data[sel_folded]  
  data_folded <- sample(data_folded, floor(length(data_folded)/10) )
  
  data_unfolded <- current_data[sel_unfolded]   
  data_unfolded <- sample(data_unfolded, floor(length(data_unfolded)/10) )
  
  current_df_folded <- data.frame(orf=rep(factor(orf_names[i]), length(data_folded)), rg=data_folded )
  current_df_unfolded <- data.frame(orf=rep(factor(orf_names[i]), length(data_unfolded)), rg=data_unfolded )
  result_folded <- rbind(result_folded, current_df_folded)
  result_unfolded <- rbind(result_unfolded, current_df_unfolded)
} 


list_orf <- unique(orf_names)
sortdf <- data.frame(orf=c(), mean=c())
for (i in 1:length(list_orf)){
  current_orf <- list_orf[i]
  current_data <- result_unfolded[result_unfolded$orf==current_orf,]
  current_mean <- mean(current_data$ap)
  current_df <- data.frame(orf=current_orf, mean=current_mean)
  sortdf <- rbind(sortdf, current_df)
}
sortdf <- sortdf[sort(sortdf$mean, index.return=T, decreasing=T)$ix, ]
result_folded$orf <- factor(result_folded$orf, levels=sortdf$orf)



svg(file = "figures/figure4/F_Rg.svg", height=2, width=3.5)

ggplot(result_folded, aes(x=rg, y=orf, fill= ..x..) ) +         
  geom_density_ridges_gradient( scale=3, rel_min_height=0.01 ) +  
  scale_fill_viridis(option="H", begin=0.1, end=0.7) + 
  scale_x_continuous(limits=c(1, 4)) + 
  labs(x="Radius of gyration (nm)", y="Proteins") +
  theme_classic() + 
  theme(
    text = element_text(size=16),
    axis.line.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.text.y = element_blank(),
    legend.position = 'none'
  )

dev.off()



####


ss <- as.data.frame(read.table("data/processed/consensus_aln_ss.txt", header=T))
ss$aln <- ss$aln + 1 # adjust python indexing

aln <- read.table("data/processed/c3718_aln.txt", sep='\t') 
names <- aln[,1]
aln <- aln[,-1]

Nter <- 37
Cter <- 392

bpinit <- as.data.frame(read.table("data/processed/assign_bpinit.txt", header=T))



N_data <- rowSums(aln[,1:37] != '-')
C_data <- rowSums(aln[,392:523] != '-')
termini <- N_data + C_data

rgm <- rep(NA, length(names))
initS <- rep(NA, length(names))
a5 <- rep(NA, length(names))
a2 <- rep(NA, length(names))
for (i in 1:length(names)){
  current_pdb <- names[i]
  rgm[i]<- mean(result_folded[result_folded$orf==current_pdb,2])
  current_init <- bpinit[bpinit$protein==current_pdb,'class']
  current_counts <- summary(current_init)
  current_prob <- current_counts / sum(current_counts)
  initS[i] <- - sum( current_prob * log(current_prob), na.rm=T)
  a5[i] <- sum(current_counts[3:5])/sum(current_counts)
  a2[i] <- sum(current_counts[1:2])/sum(current_counts)
}


rgdf <- data.frame(pdb=names, rg=rgm, initS=initS, C=C_data, N=N_data, ter=termini, a5=a5, a2=a2)



svg(file = "figures/figure4/F_correlations.svg", height=2, width=4)

p1 <- ggplot(rgdf, aes(x=rg, y=ter)) + 
  geom_point( size=2 ) + 
  geom_smooth(method=lm, alpha=0.1, size=1, color="#222222", fill="green", se=T) + 
  scale_x_continuous(limits=c(1.6, 2.5)) + 
  labs(y="Length termini", x="Radius of gyration") + 
  #geom_text(aes(label=pdb)) + 
  theme_classic() + 
  theme(
    text = element_text(size=12)
  )

# termini, not insertions

p2<- ggplot(rgdf, aes(x=a5, y=ter)) + 
  geom_point( ) + 
  geom_smooth(method=lm, alpha=0.1, size=1, color="#222222", fill="purple", se=T) + 
  scale_x_continuous(breaks=c(0.25, 0.5, 0.75, 1), labels=c(25, 50, 75, 100)) + 
  labs(x="% C-unfolding", y="Length termini") +
  #geom_text(aes(label=pdb)) + 
  theme_classic() + 
  theme(
    text = element_text(size=12)
  )

plot_grid(p1, p2, labels ="", ncol = 2, align = 'h')


dev.off()




# cor(rgdf$rg, rgdf$ter)
# 0.8329932
# fit <- lm(rgdf$ter ~ rgdf$rg)
# fit
# Coefficients:
#  (Intercept)      rgdf$rg  
# -43.47        39.82  


# cor(rgdf$a5, rgdf$ter)
# -0.5727268
# fit <- lm(rgdf$ter ~ rgdf$a5)
# fit
# Coefficients:
#  (Intercept)      rgdf$a5  
# 52.79       -27.28  


