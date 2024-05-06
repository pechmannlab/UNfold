# UNfold (2024)
# (c) S Pechmann
# Figure 5

library(ggplot2)
library(ggridges)
library(viridis)
library(cowplot)
library(reshape2)

setwd("~/UNfold/")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

svg(file = "figures/figure5/A_Cdegron.svg", height=6, width=4)


degron <- as.data.frame(read.table("data/processed/degron_probs.txt", header=T))


R_Nd70 <- round(cor(degron$Cunfold, degron$d70_N), 3)
pN1 <- ggplot(degron, aes(x=Cunfold, y=d70_N)) + 
  geom_smooth(method=lm, color="orange", fill="orange", se=F, alpha=0.3) + 
  geom_point(size=2) + 
  geom_text(x=0.8, y=6, label=paste("R =", R_Nd70, sep=' '), size=7, stat="unique") + 
  scale_x_continuous(breaks=c(0.25, 0.5, 0.75, 1), labels=c(25, 50, 75, 100)) + 
  scale_y_continuous(limits=c(-0.5, 6.5)) + 
  labs(x="% C-unfolding", y="N-ter degron score") +
  theme_classic() + 
  theme(
    text = element_text(size=20)
  )

R_Nmean <- round(cor(degron$Cunfold, degron$mean_N), 3)
pN2 <- ggplot(degron, aes(x=Cunfold, y=mean_N)) + 
  geom_smooth(method=lm, color="orange", fill="orange", se=F, alpha=0.3) + 
  geom_point(size=2) + 
  geom_text(x=0.8, y=0.6, label=paste("R =", R_Nmean, sep=' '), size=7, stat="unique") + 
  scale_x_continuous(breaks=c(0.25, 0.5, 0.75, 1), labels=c(25, 50, 75, 100)) + 
  labs(x="% C-unfolding", y="C-ter avrg. degron prob.") +
  theme_classic() + 
  theme(
    text = element_text(size=20)
  )

R_Cd70 <- round(cor(degron$Cunfold, degron$d70_C), 3)
pC1 <- ggplot(degron, aes(x=Cunfold, y=d70_C)) + 
  geom_smooth(method=lm, color="purple", fill="purple", se=F, alpha=0.3) + 
  geom_point(size=2) + 
  geom_text(x=0.8, y=6, label=paste("R =", R_Cd70, sep=' '), size=7, stat="unique") + 
  scale_x_continuous(breaks=c(0.25, 0.5, 0.75, 1), labels=c(25, 50, 75, 100)) + 
  scale_y_continuous(limits=c(-0.5, 6.5)) + 
  labs(x="% C-unfolding", y="C-ter degron score") +
  theme_classic() + 
  theme(
    text = element_text(size=20)
  )

R_Cmean <- round(cor(degron$Cunfold, degron$mean_C), 3)
pC2 <- ggplot(degron, aes(x=Cunfold, y=mean_C)) + 
  geom_smooth(method=lm, color="purple", fill="purple", se=F, alpha=0.3) + 
  geom_point(size=2) + 
  geom_text(x=0.8, y=0.6, label=paste("R =", R_Cmean, sep=' '), size=7, stat="unique") + 
  scale_x_continuous(breaks=c(0.25, 0.5, 0.75, 1), labels=c(25, 50, 75, 100)) + 
  labs(x="% C-unfolding", y="C-ter avrg. degron prob.") +
  theme_classic() + 
  theme(
    text = element_text(size=20)
  )

#plot_grid(pN1, pC1, pN2, pC2, labels ="", ncol = 2, align = 'h')
plot_grid(pC1,  pC2, labels ="", ncol = 1, align = 'v')


dev.off()








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


svg(file = "figures/figure5/B_Ndegron.svg", height=6, width=3)


degron <- as.data.frame(read.table("data/processed/degron_probs.txt", header=T))
degron$substr <- factor(degron$substr, levels=c(1, 0))


p1 <- ggplot(degron, aes(y=d70_Np)) + 
  geom_boxplot(aes(x=substr, fill=substr )) + 
  labs(x="", y="Degron score") + 
  scale_x_discrete(breaks=c(1, 0), labels=c("Hsp70+", "Hsp70-")) + 
  scale_fill_manual(values=c("#29ABE2", "#FF0000")) + 
  theme_classic() +
  theme(
    text = element_text(size=16), 
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    legend.position = 'none'
  )

p2 <- ggplot(degron, aes(y=mean_Np)) + 
  geom_boxplot(aes(x=substr, fill=substr )) + 
  labs(x="", y="Avrg. degron prob.") + 
  scale_x_discrete(breaks=c(1, 0), labels=c("Hsp70+", "Hsp70-")) + 
  scale_fill_manual(values=c("#29ABE2", "#FF0000")) + 
  theme_classic() +
  theme(
    text = element_text(size=16), 
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    legend.position = 'none'
  )

p3 <- ggplot(degron, aes(y=N)) + 
  geom_boxplot(aes(x=substr, fill=substr )) + 
  labs(x="", y="Length N-ter") + 
  scale_x_discrete(breaks=c(1, 0), labels=c("Hsp70+", "Hsp70-")) + 
  scale_fill_manual(values=c("#29ABE2", "#FF0000")) + 
  theme_classic() +
  theme(
    text = element_text(size=16), 
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    legend.position = 'none'
  )



plot_grid(p1, p2, p3, labels ="", ncol = 1, align = 'v')

dev.off()




boxplot(degron$d70_Np[degron$substr==1], degron$d70_Np[degron$substr==0])
wilcox.test(degron$d70_Np[degron$substr==1], degron$d70_Np[degron$substr==0])

boxplot(degron$N[degron$substr==1], degron$N[degron$substr==0])
wilcox.test(degron$N[degron$substr==1], degron$N[degron$substr==0])





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# !! always reload vars, as re-using variable names for different data !! 

svg(file = "figures/figure5/C_1_unfolding.svg", height=3, width=6)

pred.fullfit <- as.data.frame(read.table("data/processed/unfold_fullfit.txt", header=T))
cor.ff <- round(cor(pred.fullfit$true, pred.fullfit$pred),3)

plot.ff <- ggplot(pred.fullfit, aes(x=true, y=pred)) + 
  labs(x="Unfolding (true)", y="Unfolding (fitted)") + 
  geom_point( size=2) + 
  geom_smooth(method=lm, alpha=0.2, size=1, color="#222222", fill="#777777", se=T) + 
  geom_text(x=0.4, y=0.8, label=paste("R =", cor.ff, sep=' '), size=5, stat="unique") + 
  theme_classic() +
  theme(
    text = element_text(size=16)
  )


pred.xval <- as.data.frame(read.table("data/processed/unfold_xval.txt", header=T))
res.0 <- pred.xval[1,]        # first data point not part of fitting !! 
pred.xval <- pred.xval[-1,]
cor.xval <- round(cor(pred.xval$true, pred.xval$pred),3)

plot.xval <- ggplot(pred.xval, aes(x=true, y=pred)) + 
  geom_point( size=2) + 
  labs(x="Unfolding (true)", y="Unfolding (predicted)") + 
  geom_point(x=res.0$true, y=res.0$pred, size=2, col="red", shape=17) + 
  geom_smooth(method=lm, alpha=0.2, size=1, color="#222222", fill="#0a5dd3", se=T) + 
  geom_text(x=0.4, y=0.535, label=paste("R =", cor.xval, sep=' '), size=5, stat="unique") + 
  theme_classic() +
  theme(
    text = element_text(size=16)
  )


plot_grid(plot.ff, plot.xval, labels ="", ncol = 2, align = 'h')

dev.off()


###


svg(file = "figures/figure5/C_2_unfold_fts.svg", height=3, width=6)

fa.fullfit <- scan("data/processed/unfold_fullfit_fa.txt")
fa.xval <- read.table("data/processed/unfold_xval_fa.txt")

fa.fullfit.sd <- rep(0, 5)
fa.xval.sd <- rep(NA, 5)
for (i in 1:5){fa.xval.sd[i] <- sd(fa.xval[,i])}

fa <- data.frame(ff=fa.fullfit, ff_sd=fa.fullfit.sd, xv=colMeans(fa.xval), xv_sd=fa.xval.sd)
fa$feature  <- c('Abundance', 'Aggregation', 'Hsp70', 'Tf', 'No. contacts')

fa.means <- fa[,c('feature', 'ff', 'xv')]
fa.means.m <- melt(fa.means, id='feature')
fa.sd <- fa[,c('feature', 'ff_sd', 'xv_sd')]
fa.sd.m <- melt(fa.sd, id='feature')
fa.means.m$sd <- fa.sd.m$value



ggplot(fa.means.m, aes(x=feature, y=value)) + 
  geom_col(aes(fill=variable), position="dodge2", alpha=0.7) + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), position=position_dodge2(width = 0.2, padding = 0.8), col="#222222", alpha=0.7) + 
  labs(x="", y="Feature importance") +
  scale_fill_manual(values=c("#777777", "#0a5dd3")) + 
  theme_classic() +
  theme(
    text = element_text(size=20), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle=30, hjust=0.5, vjust=0), 
    legend.position = 'none'
  )

dev.off()

### substr

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


svg(file = "figures/figure5/D_substr.svg", height=3, width=6)

substr.fullfit <- as.data.frame(read.table("data/processed/substr_fullfit.txt", header=T))
substr.xval <- as.data.frame(read.table("data/processed/substr_xval.txt", header=T))

tree.fullfit <- data.frame(class=c("+T", "+F", "-T", "-F"), 
           fit = rep("fullfit",4),
           pred =c(
             sum( substr.fullfit$pred[substr.fullfit$true==1] == 1 ),
             sum( substr.fullfit$pred[substr.fullfit$true==1] == 0 ),
             sum( substr.fullfit$pred[substr.fullfit$true==0] == 0 ),
             sum( substr.fullfit$pred[substr.fullfit$true==0] == 1 )
           ))
tree.fullfit$class <- factor(tree.fullfit$class, levels=c("+T", "+F", "-T", "-F"))

tree.xval <- data.frame(class=c("+T", "+F", "-T", "-F"), 
                           fit = rep("xval",4),
                           pred =c(
                             sum( substr.xval$pred[substr.xval$true==1] == 1 ),
                             sum( substr.xval$pred[substr.xval$true==1] == 0 ),
                             sum( substr.xval$pred[substr.xval$true==0] == 0 ),
                             sum( substr.xval$pred[substr.xval$true==0] == 1 )
                             
                           ))
tree.xval$class <- factor(tree.xval$class, levels=c("+T", "+F", "-T", "-F"))


t1 <- ggplot(tree.fullfit , aes(x=class, y=pred+0.05)) + 
  geom_col(aes(fill=class), position="dodge2") + 
  scale_fill_manual(values=c("#29ABE2", "grey50", "#FF0000",  "grey70")) + 
  labs(x="", y="Hsp70 int. (fitted)") +
  theme_classic() +
  theme(
    text = element_text(size=15),
    axis.text.x = element_text(angle=90, hjust=0.5, vjust=0.9),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  )

t2 <- ggplot(tree.xval , aes(x=class, y=pred+0.05)) + 
  geom_col(aes(fill=class), position="dodge2") + 
  scale_fill_manual(values=c("#29ABE2", "grey50", "#FF0000",  "grey70")) + 
  scale_y_continuous(breaks=c(0, 3, 6, 9), labels=c(0, 3, 6, 9)) + 
  labs(x="", y="Hsp70 int. (pred.)") +
  theme_classic()  +
  theme(
    text = element_text(size=15),
    axis.text.x = element_text(angle=90, hjust=0.5, vjust=0.9),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  )


fa.fullfit <- scan("data/processed/substr_fullfit_fa.txt")
fa.xval <- read.table("data/processed/substr_xval_fa.txt")
fa.fullfit.sd <- rep(0, 4)
fa.xval.sd <- rep(NA, 4)
for (i in 1:4){fa.xval.sd[i] <- sd(fa.xval[,i])}

fa <- data.frame(ff=fa.fullfit, ff_sd=fa.fullfit.sd, xv=colMeans(fa.xval), xv_sd=fa.xval.sd)
fa$feature  <- c('Abundance', 'Aggregation', 'Unfolding', 'No. contacts')

fa.means <- fa[,c('feature', 'ff', 'xv')]
fa.means.m <- melt(fa.means, id='feature')
fa.sd <- fa[,c('feature', 'ff_sd', 'xv_sd')]
fa.sd.m <- melt(fa.sd, id='feature')
fa.means.m$sd <- fa.sd.m$value



p.fts <- ggplot(fa.means.m, aes(x=feature, y=value)) + 
  geom_col(aes(fill=variable), position="dodge2", alpha=0.7) + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), position=position_dodge2(width = 0.2, padding = 0.8), col="#222222", alpha=0.7) + 
  labs(x="", y="Feature importance") +
  scale_fill_manual(values=c("#777777", "#0a5dd3")) + 
  theme_classic() +
  theme(
    text = element_text(size=16), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle=30, hjust=0.5, vjust=0.7), 
    legend.position = 'none', #c(0.8, 0.9), 
    legend.title = element_blank()
  )


plot_grid(t1,  t2, p.fts, labels ="", rel_widths = c(1, 1, 2), ncol = 3, align = 'h')

dev.off()


