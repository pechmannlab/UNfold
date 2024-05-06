# UNfold (2024)
# (c) S Pechmann
# Figure 2

library(ggplot2)
library(reshape2)


setwd("~/UNfold/")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#pymol structure

frust <- as.data.frame(read.table("data/processed/frustration_msaprofile.txt", header=T))
frust$ymin <- ifelse(frust$mean - frust$sd > 0, frust$mean - frust$sd, 0)
frust$ymax <- frust$mean + frust$sd

ss <- as.data.frame(read.table("data/processed/consensus_aln_ss.txt", header=T))
ss$aln <- ss$aln + 1 # adjust python indexing


# get start and end of secondary structure elements 
last_element <- "placeholder"
counter_loop <- 1
ssdf <- data.frame()
for (i in 1:(nrow(ss)-1)){
  current_element <- as.vector(ss$element)[i]
  next_element <- as.vector(ss$element)[i+1]
  current_ss <- as.vector(ss$ss)[i]
 
  if (current_element != last_element){
    res_start <- ss$aln[i]
    last_element <- current_element
    if (current_ss == "loop") {
      res_ss <- paste("loop", counter_loop, sep="")
      counter_loop <- counter_loop + 1
    } else {
        res_ss <- current_element
      }
    }
  if (current_element != next_element){
    res_end <- ss$aln[i]
    ssdf <- rbind(ssdf, data.frame(start=res_start, end=res_end, name=res_ss ) )
  }
  
  if (i == nrow(ss)-1){
    res_end <- ss$aln[i+1]
    ssdf <- rbind(ssdf, data.frame(start=res_start, end=res_end, name=res_ss ) )
  }
}

# parse coordinates for polygon plotting
polyx <- c()
polyy <- c()
polyn <- c()
for (i in 1:nrow(ssdf)){
  current_start <- ssdf$start[i]
  current_end <- ssdf$end[i]
  current_name <- as.vector(ssdf$name)[i]
  current_x <- c(current_start, current_end, current_end, current_start  )
  current_y <-  rep(c(0, 1), each=2)
  polyx <- c(polyx, current_x)
  polyy <- c(polyy, current_y)
  polyn <- c(polyn, current_name)
}

# parse data frame for polygon plotting and set levels/order
sspoly <- data.frame(id=rep(polyn, each=4), x=polyx, y=polyy )
colo = c()
level =c()
for (i in 1:(nrow(sspoly)-1) ){
  current_ss <- as.vector(sspoly$id)[i]
  next_ss    <- as.vector(sspoly$id)[i+1]

  if (current_ss != next_ss) {
    level <- c(level, current_ss)
    if (length(grep("loop", current_ss))==1){
      colo <- c(colo, "#777777")
    } else if (length(grep("b", current_ss))==1){
      colo <- c(colo, "#b735b7")
    } else if (length(grep("a", current_ss))==1){
      colo <- c(colo, "red")
    } else {colo <- c(colo, "white")}

  } else if (i == (nrow(sspoly)-1) ) {
    level <- c(level, next_ss)
    if (length(grep("loop", current_ss))==1){
      colo <- c(colo, "#777777")
    } else if (length(grep("b", current_ss))==1){
      colo <- c(colo, "#b735b7")
    } else if (length(grep("a", current_ss))==1){
      colo <- c(colo, "red")
    } else {colo <- c(colo, "white")}
  }
}
sspoly$id <- factor(sspoly$id, levels=level)




svg(file = "figures/figure2/B_frustration_aln.svg", height = 3, width = 6)

ggplot(frust, aes(x=pos, y=mean, ymin=ymin, ymax=ymax)) + 
  geom_polygon(data=sspoly, inherit.aes=F, aes(x=x, y=y, fill=id), show.legend= FALSE, alpha=0.1) + 
  scale_fill_manual(values=colo) +
  geom_line() + 
  geom_ribbon(alpha=0.5, fill="#006600") + 
  labs(x="Alignment position", y="Frustration") + 
  theme_classic() +
  theme(
    text = element_text(size=20)
  )

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

svg("figures/figure2/C_frustration_ss.svg", height=4, width=4)

frustss <- as.data.frame(read.table("data/processed/frustration_aligned_ss.txt", header=T))
frustss.m <- melt(frustss)

sstype <- rep(NA, nrow(frustss.m))   #add ss type for coloring
for (i in 1:nrow(frustss.m)){
  current_ss <- as.vector(frustss.m$variable[i])
  current_type <- as.vector(ss[ss$element == current_ss, 2])
  current_type <- current_type[1]
  sstype[i] <- current_type
}
frustss.m$type <- sstype


ggplot(frustss.m, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=type)) + 
  geom_jitter(size=0.5, alpha=0.6, width=0.2, color="#333333") + 
  labs(y="Average frustration", x="") + 
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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

svg("figures/figure2/D_frustration_aln.svg", height=4, width=2.4)


f2 <- as.data.frame(read.table("data/processed/frustration_profiles.txt", header=T))
sel_aln <- colSums(is.na(f2)==F) == 16
sel_not <- colSums(is.na(f2)==T)  > 0
f <- data.frame(pdb=row.names(f2), not=rowMeans(f2[,sel_not], na.rm=T), aln=rowMeans(f2[,sel_aln]))


f.m <- melt(f)
f.m$variable <- factor(f.m$variable, levels=c('aln', 'not'))

ggplot(f.m, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable), alpha=0.7) + 
  geom_jitter(size=0.8, alpha=0.6, width=0.2, color="red") + 
  scale_fill_manual(values=c("#006600", "#777777")) +   
  labs(x="", y="Avrg. frustration") + 
  theme_classic() + 
  theme(
    text = element_text(size=24), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'none' 
  )

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#pymol structure


rmsf <- as.data.frame(read.table("data/processed/rmsf_msaprofile.txt", header=T))
rmsf$ymin <- ifelse(rmsf$mean - rmsf$sd > 0, rmsf$mean - rmsf$sd, 0)
rmsf$ymax <- rmsf$mean + rmsf$sd

# !! requires sspoly from above panel B !! 
sspoly_rmsf <- sspoly
sspoly_rmsf$y <- ifelse(sspoly_rmsf$y==1, 6, 0)

svg(file = "figures/figure2/F_rmsf_aln.svg", height = 3, width = 6)

ggplot(rmsf, aes(x=pos, y=mean, ymin=ymin, ymax=ymax)) + 
  scale_y_continuous(limits=c(0, 6.4)) + 
  geom_polygon(data=sspoly_rmsf, inherit.aes=F, aes(x=x, y=y, fill=id), show.legend= FALSE, alpha=0.1) + 
  scale_fill_manual(values=colo) +
  geom_line() + 
  geom_ribbon(alpha=0.5, fill="#da9a24") + 
  labs(x="Alignment position", y="Rmsf") + 
  theme_classic() +
  theme(
    text = element_text(size=20)
  )

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

svg("figures/figure2/G_rmsf_ss.svg", height=4, width=4)

rmsfss <- as.data.frame(read.table("data/processed/rmsf_aligned_ss.txt", header=T))
rmsfss.m <- melt(rmsfss)

sstype <- rep(NA, nrow(rmsfss.m))   #add ss type for coloring
for (i in 1:nrow(rmsfss.m)){
  current_ss <- as.vector(rmsfss.m$variable[i])
  current_type <- as.vector(ss[ss$element == current_ss, 2])
  current_type <- current_type[1]
  sstype[i] <- current_type
}
rmsfss.m$type <- sstype


ggplot(rmsfss.m, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=type)) + 
  geom_jitter(size=0.5, alpha=0.6, width=0.2, color="#333333") + 
  labs(y="Average rmsf", x="") + 
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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL F ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

svg("figures/figure2/H_rmsf_aln.svg", height=4, width=2.4)


f2 <- as.data.frame(read.table("data/processed/rmsf_profiles.txt", header=T))
sel_aln <- colSums(is.na(f2)==F) == 16
sel_not <- colSums(is.na(f2)==T)  > 0
f <- data.frame(pdb=row.names(f2), not=rowMeans(f2[,sel_not], na.rm=T), aln=rowMeans(f2[,sel_aln]))

f.m <- melt(f)
f.m$variable <- factor(f.m$variable, levels=c('aln', 'not'))

ggplot(f.m, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable), alpha=0.7) + 
  geom_jitter(size=0.8, alpha=0.6, width=0.2, color="red") + 
  scale_fill_manual(values=c("#da9a24", "#777777")) +   
  labs(x="", y="Avrg. rmsf") + 
  theme_classic() + 
  theme(
    text = element_text(size=24), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'none' 
  )

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL G ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

svg("figures/figure2/G_corrmotion.svg", height=4, width=4)


cm <- as.data.frame(read.table("code/tmp.cm", header=T))
clustID <- cm[,1]
cm <- cm[,-1]

h <- heatmap(as.matrix(cm))

cm <- cm[h$rowInd,]
cm <- cm[, h$colInd]
cm[cm<5] <- 0

clustID <- clustID[h$rowInd]

df <- data.frame(i=c(), j=c(), n=c())
for (i in 1:nrow(cm)){
  for (j in 1:ncol(cm)){
    current_n = cm[i,j]
    if (current_n > 0){
      current_df = data.frame(i=i, j=j, n=current_n)
      df <- rbind(df, current_df)
    }
  }
}
df$n <- factor(df$n)


ggplot(df, aes(x=i, y=j)) + 
  geom_point(aes(color=n), shape=20) + 
  theme_classic() + 
  #scale_x_continuous(breaks=c(1:length(clustID)), labels=clustID) + 
  #scale_y_continuous(breaks=c(1:length(clustID)), labels=clustID) + 
  scale_colour_viridis_d(option = "E", begin=0.3, end=1) + 
  xlab("Contact cluster") + 
  ylab("Contact cluster") + 
  theme(
    text=element_text(size=20),
    axis.text = element_blank(),
    #axis.text.x = element_text(size=8, angle=90, hjust=1, vjust=0.5),  
    #axis.text.y = element_text(size=8, hjust=1, vjust=0.5),  
    legend.position = 'right'
    )

dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL H ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

frustp <- as.data.frame(read.table("data/processed/frustration_profiles.txt", header=T))
rmsfp <- as.data.frame(read.table("data/processed/rmsf_profiles.txt", header=T))

ss <- as.data.frame(read.table("data/processed/consensus_aln_ss.txt", header=T))
ss$aln <- ss$aln + 1 # adjust python indexing

sel_b1 <- ss[ss$element=="b1",1]
sel_b2 <- ss[ss$element=="b2",1]
sel_b3 <- ss[ss$element=="b3",1]
sel_b4 <- ss[ss$element=="b4",1]
sel_b5 <- ss[ss$element=="b5",1]
sel_b6 <- ss[ss$element=="b6",1]
sel_b <- c(sel_b1, sel_b2, sel_b3, sel_b4, sel_b5, sel_b6)

sel_a1 <- ss[ss$element=="a1",1]
sel_a2 <- ss[ss$element=="a2",1]
sel_a3 <- ss[ss$element=="a3",1]
sel_a4 <- ss[ss$element=="a4",1]
sel_a5 <- ss[ss$element=="a5",1]
sel_a <- c(sel_a1, sel_a2, sel_a3, sel_a4, sel_a5)


datap <- rbind(
  data.frame(rmsf=colMeans(rmsfp[,sel_a]), frust=colMeans(frustp[,sel_a]), ss=rep("a1-5", length(sel_a)) ),
  data.frame(rmsf=colMeans(rmsfp[,sel_b]), frust=colMeans(frustp[,sel_b]), ss=rep("b1-6", length(sel_b)) )
)


svg("figures/figure2/H_frust_rmsf.svg", height=4, width=5)

ggplot(datap, aes(x=frust, y=rmsf)) + 
  geom_smooth(method=lm, aes(color=ss, fill=ss), alpha=0.2, size=1, linetype='dashed', se=T) + 
  geom_point(aes(color=ss), size=2) +
  labs(x="Frustration", y="Rmsf") + 
  scale_color_manual(values=c("red", "purple")) + 
  scale_fill_manual(values=c("grey70", "grey70")) + 
  theme_classic() +
  theme(
    text = element_text(size=20),
    legend.position = c(0.83, 0.2),
    legend.title = element_blank()
  )

dev.off()

