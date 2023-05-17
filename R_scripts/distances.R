library(tidyverse)
library(rio)
library(RColorBrewer)
library(FSA)
library(ggpubr)
library(here)

load("junction_stats.RData")
distance.list <- junction.stats
names(distance.list) <- gsub("user_","distance_",names(distance.list))

dvg.distances <- bind_rows(distance.list)
dvg.totals <- aggregate(accumulation ~ start_pos + end_pos,data = dvg.distances,
                        FUN = "sum")

#significance testing
dvg.range <- range(dvg.distances$distance)
low.cutoff <- dvg.range[1]+ (0.01*(dvg.range[2]-dvg.range[1])) #first 1% of distance range
high.cutoff <-  dvg.range[2] - (0.15*(dvg.range[2]-dvg.range[1])) #last 15% of distance range
dvg.distances <- dvg.distances %>%
  mutate("Bin"= NA)

for(i in seq_along(dvg.distances$distance)){
  if(0 < dvg.distances$distance[[i]] && dvg.distances$distance[[i]] <= low.cutoff){
    dvg.distances$Bin[[i]] <- "1"
  }
  if(low.cutoff < dvg.distances$distance[[i]] && dvg.distances$distance[[i]] <= high.cutoff){
    dvg.distances$Bin[[i]] <- "2"
  }
  if(high.cutoff < dvg.distances$distance[[i]] && dvg.distances$distance[[i]] <= 30000){
    dvg.distances$Bin[[i]] <- "3"
  }
}
dvg.distances$Bin <- as.factor(dvg.distances$Bin)

dunn.accumulation <- dunnTest(accumulation ~ Bin, data = dvg.distances,
                              method = "bonferroni")
#1 vs 2: p.adj = 0.0324428071
#1 vs 3 p.adj = 0.0561867796
#2 vs 3: p.adj = 0.0005980824

dunn.enrich <- dunnTest(delta_enrich ~ Bin, data = dvg.distances,
                        method = "bonferroni")
#1 vs 2: p.adj = 0.075049122
#1 vs 3 p.adj = 0.008308182
#2 vs 3: p.adj = 1.000000000

#medians
bin1 <- dvg.distances[which(dvg.distances$Bin == "1"),]
bin2 <- dvg.distances[which(dvg.distances$Bin == "2"),]
bin3 <- dvg.distances[which(dvg.distances$Bin == "3"),]
median(bin1$accumulation) #25.74515
median(bin2$accumulation) #19.63858
median(bin3$accumulation) #65.83058

median(bin1$delta_enrich) #0.1046606
median(bin2$delta_enrich) #0.1600003
median(bin3$delta_enrich) #0.1571062

#junction length vs total accumulation
ggplot(data = dvg.distances, aes(x = distance, y = accumulation, color = Bin))+
  geom_point(cex = 3, shape = 1)+
  xlab("Deletion length")+
  ylab("DVG accumulation")+
  scale_x_continuous(limits = c(0,30000), breaks = c(0,10000,20000,30000))+
  scale_y_continuous(limits = c(0,1500), breaks = c(0,500,1000,1500))+
  scale_color_manual(values = c("red","black","blue"))+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19))
ggsave("figs/accumulation.png")

acc.matrix <- as.data.frame(matrix(nrow = 3, ncol = 4, data = NA))
colnames(acc.matrix) <- c("group1","group2","p.adj","p.sig")
acc.matrix$group1 <- c("1","1","2")
acc.matrix$group2 <- c("2","3","3")
acc.matrix$p.adj <- dunn.accumulation$res$P.adj
acc.matrix$p.sig <- c("*","ns","***")

ggplot(data = dvg.distances, aes(x = Bin, y = accumulation, color = Bin))+
  geom_jitter(cex = 3, shape = 1)+
  ylab("Accumulation")+
  scale_y_continuous(limits = c(-0.01,1500), breaks = c(0,500,1000,1500))+
  scale_color_manual(values = c("red","black","blue"))+
  stat_pvalue_manual(acc.matrix, label = "p.sig", y.position = c(1450,1200,950), size = 7)+
  stat_summary(fun = "median", geom = "crossbar", width = .5, color = "black")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        legend.position = "none")
ggsave("figs/accumulation_bar.png")

#junction length vs enrichment over time
ggplot(data = dvg.distances, aes(x = distance, y = delta_enrich, color = Bin))+
  geom_point(cex = 3, shape = 1)+
  xlab("Deletion length")+
  ylab("DVG enrichment")+
  scale_x_continuous(limits = c(0,30000), breaks = c(0,10000,20000,30000))+
  scale_color_manual(values = c("red","black","blue"))+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19))
ggsave("figs/enrichment.png")

en.matrix <- as.data.frame(matrix(nrow = 3, ncol = 4, data = NA))
colnames(en.matrix) <- c("group1","group2","p.adj","p.sig")
en.matrix$group1 <- c("1","1","2")
en.matrix$group2 <- c("2","3","3")
en.matrix$p.adj <- dunn.enrich$res$P.adj
en.matrix$p.sig <- c("ns","**","ns")

ggplot(data = dvg.distances, aes(x = Bin, y = delta_enrich, color = Bin))+
  geom_jitter(cex = 3, shape = 1)+
  ylab("Enrichment")+
  scale_y_continuous(limits = c(-0.01,1.5), breaks = c(0,.25,.5,.75,1))+ 
  scale_color_manual(values = c("red","black","blue"))+
  stat_pvalue_manual(en.matrix, label = "p.sig", y.position = c(1.45,1.3,1.15), size = 7)+
  stat_summary(fun = "median", geom = "crossbar", width = .5, color = "black")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        legend.position = "none")
ggsave("figs/enrich_bar.png")

#junction positions vs TRS positions
positions <- dvg.totals %>%
  arrange(desc(accumulation))
positions <- positions[,-3]
sgmrna <- read.csv("sgmrna.csv", header = T)
sgmrna <- sgmrna %>%
  mutate("segment" = 11:19)
top <- positions[c(1:10),]
top <- top %>%
  mutate("segment" = 10:1)

ggplot(top) + 
  geom_segment(aes(x=start_pos,xend = end_pos, y = segment, yend = segment),
               size=1, linetype = 2)+
  geom_point(aes(x=start_pos,y = segment),color = "blue",size=4)+
  geom_point(aes(x=end_pos, y = segment), color = "red", size=4)+
  geom_segment(data = sgmrna,aes(x=start_pos,xend = end_pos, 
                                 y = segment,yend = segment),
               size=1, color = "dark grey", linetype = 2)+
  geom_point(data = sgmrna,aes(x=start_pos,y = segment),
               size=4, color = "light blue")+
  geom_point(data = sgmrna,aes(x=end_pos,y = segment),
             size=4, color = "light pink")+
  xlab("Junction")+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 19),
        axis.title.x = element_text(size = 22))
ggsave("figs/sgmrna_vs_junctions.png")

#junction start and end distribution (total reads accumulation)
dvg.totals <- arrange(dvg.totals, accumulation, desc = F)
pal <- RColorBrewer::brewer.pal(3, "YlOrRd")
ggplot(dvg.totals, aes(x= start_pos, y = end_pos, color = accumulation))+
  geom_point(cex = 3)+
  xlab("Start position")+
  ylab("End position")+
  scale_color_gradientn(colours = pal,name = "Accumulation")+
  theme_bw()+
  theme(axis.title = element_text(size = 24,
                                  face = "bold"),
        axis.text = element_text(size = 21,
                                 face = "bold"),
        legend.title = element_text(size=24,
                                    face = "bold"),
        legend.text = element_text(size = 21,
                                   face = "bold"),
        plot.title = element_text(size=24,
                                  face = "bold"))
ggsave("figs/junction_reads.png")

#junction start and end distribution (per individual)
dvg.distances <- dvg.distances %>%
  mutate("junction" =  paste0(start_pos,"_",end_pos)) %>%
  relocate(junction,.before = start_pos)
junction.counts <- table(dvg.distances$junction)
junction.counts <- as.data.frame(junction.counts)
colnames(junction.counts) <- c("junction","Frequency")
all.together <- full_join(dvg.distances,junction.counts,by="junction")
all.together <- all.together %>%
  arrange(Frequency, desc = F)

ggplot(all.together, aes(x=start_pos,y=end_pos, color=Frequency))+
  geom_point(cex=3)+
  xlab("Start position")+
  ylab("End position")+
  scale_color_gradientn(colours = pal,name = "Number of \nparticipants")+
  theme_bw()+
  theme(axis.title = element_text(size = 24,
                                  face = "bold"),
        axis.text = element_text(size = 21,
                                 face = "bold"),
        legend.title = element_text(size=24,
                                    face = "bold"),
        legend.text = element_text(size = 21,
                                   face = "bold"),
        plot.title = element_text(size=24,
                                  face = "bold"))
ggsave("figs/junction_occurrences.png")

#frequency histogram
ggplot(all.together, aes(x = Frequency))+
  geom_histogram()+
  xlab("Number of participants \n with junction")+
  ylab("Junction count")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 19, face = "bold"),
        axis.title.x = element_text(size = 22, face = "bold"),
        axis.text.y = element_text(size = 19, face = "bold"),
        axis.title.y = element_text(size = 22, face = "bold"))
ggsave("figs/junction_histogram.png")



