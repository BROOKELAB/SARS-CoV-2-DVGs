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
#bin by percentile of deletion length range
dvg.range <- range(dvg.distances$distance)
cutoff <- dvg.range[1]+ (0.01*(dvg.range[2]-dvg.range[1])) #first 1% of distance range

dvg.distances <- dvg.distances %>%
  mutate("Bin"= NA)
dvg.distances[which(dvg.distances$distance <= cutoff),]$Bin <- "1%"
length(which(dvg.distances$Bin == "1%")) #255 junctions
dvg.distances[which(dvg.distances$distance > cutoff),]$Bin <- "2-100%"
length(which(dvg.distances$Bin == "2-100%")) #54 junctions

dvg.distances$Bin <- as.factor(dvg.distances$Bin)

median(dvg.distances[which(dvg.distances$Bin == "1%"),]$accumulation) #25.74515
median(dvg.distances[which(dvg.distances$Bin == "2-100%"),]$accumulation) #17.77993

median(dvg.distances[which(dvg.distances$Bin == "1%"),]$delta_enrich) #0.1282061
median(dvg.distances[which(dvg.distances$Bin == "2-100%"),]$delta_enrich) #0.2105919

kw.acc <- kruskal.test(dvg.distances$accumulation ~ dvg.distances$Bin) #p-value = 0.001792
kw.enr <- kruskal.test(dvg.distances$delta_enrich ~ dvg.distances$Bin) #p-value = 0.031

#deletion length vs total accumulation
ggplot(data = dvg.distances, aes(x = distance, y = accumulation, color = Bin))+
  geom_point(cex = 3, shape = 1)+
  xlab("Deletion length")+
  ylab("DVG accumulation")+
  scale_x_continuous(limits = c(0,30000), breaks = c(0,10000,20000,30000))+
  scale_y_continuous(limits = c(0,1500), breaks = c(0,500,1000,1500))+
  scale_color_manual(values = c("red","blue"))+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19))
ggsave("figs/accumulation.png")

acc.matrix <- as.data.frame(matrix(nrow = 1, ncol = 4, data = NA))
colnames(acc.matrix) <- c("group1","group2","p.adj","p.sig")
acc.matrix$group1 <- c("1%")
acc.matrix$group2 <- c("2-100%")
acc.matrix$p.adj <- kw.acc$p.value
acc.matrix$p.sig <- c("**")

set.seed(1)
ggplot(data = dvg.distances, aes(x = Bin, y = accumulation, color = Bin))+
  geom_jitter(cex = 3, shape = 1)+
  ylab("Accumulation")+
  scale_y_continuous(limits = c(-0.01,1500), breaks = c(0,500,1000,1500))+
  scale_color_manual(values = c("red","blue"))+
  stat_pvalue_manual(acc.matrix, label = "p.sig", y.position = c(1400), size = 7)+
  stat_summary(fun = "median", geom = "crossbar", width = .5, color = "black")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        legend.position = "none")
ggsave("figs/accumulation_bar.png")

#deletion length vs enrichment over time
ggplot(data = dvg.distances, aes(x = distance, y = delta_enrich, color = Bin))+
  geom_point(cex = 3, shape = 1)+
  xlab("Deletion length")+
  ylab("DVG enrichment")+
  scale_x_continuous(limits = c(0,30000), breaks = c(0,10000,20000,30000))+
  scale_color_manual(values = c("red","blue"))+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19))
ggsave("figs/enrichment.png")

en.matrix <- as.data.frame(matrix(nrow = 1, ncol = 4, data = NA))
colnames(en.matrix) <- c("group1","group2","p.adj","p.sig")
en.matrix$group1 <- c("1%")
en.matrix$group2 <- c("2-100%")
en.matrix$p.adj <- kw.enr$p.value
en.matrix$p.sig <- c("*")

set.seed(1)
ggplot(data = dvg.distances, aes(x = Bin, y = delta_enrich, color = Bin))+
  geom_jitter(cex = 3, shape = 1)+
  ylab("Enrichment")+
  scale_y_continuous(limits = c(-0.01,1.25), breaks = c(0,.25,.5,.75,1))+ 
  scale_color_manual(values = c("red","blue"))+
  stat_pvalue_manual(en.matrix, label = "p.sig", y.position = c(1.2), size = 7)+
  stat_summary(fun = "median", geom = "crossbar", width = .5, color = "black")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        legend.position = "none")
ggsave("figs/enrich_bar.png")

#shared junction positions vs TRS positions
shared.table <- read.csv("shared_table.csv")
split <- strsplit(shared.table$junction,"_")
starts <- as.numeric(sapply(split,"[",1))
ends <- as.numeric(sapply(split,"[",2))
shared.table <- shared.table %>%
  mutate("start_pos" = starts)%>%
  mutate("end_pos" = ends)%>%
  select(start_pos,end_pos)%>%
  arrange(start_pos)%>%
  mutate(segment = 19-as.numeric(rownames(shared.table)))

sgmrna <- read.csv("sgmrna.csv", header = T)
sgmrna <- sgmrna %>%
  mutate("segment" = 20:29)

ggplot(shared.table) + 
  geom_segment(aes(x=start_pos,xend = end_pos, y = segment, yend = segment),
               size=0.5, linetype = 2)+
  geom_segment(aes(x = 0, xend = start_pos, y = segment, yend = segment), 
               size=0.5, linetype = 1)+
  geom_segment(aes(x = end_pos, xend = 29903, y = segment, yend = segment), 
               size=0.5, linetype = 1)+
  geom_point(aes(x=start_pos,y = segment),color = "red",size=1.5)+
  geom_point(aes(x=end_pos, y = segment), color = "red", size=1.5)+
  
  geom_segment(data = sgmrna,aes(x=start_pos,xend = end_pos, y = segment, yend = segment),
               size=0.5, linetype = 2)+
  geom_segment(data = sgmrna, aes(x = 0, xend = start_pos, y = segment, yend = segment), 
               size=0.5, linetype = 1)+
  geom_segment(data = sgmrna, aes(x = end_pos, xend = 29903, y = segment, yend = segment), 
               size=0.5, linetype = 1)+
  geom_point(data = sgmrna,aes(x=start_pos,y = segment),size=1.5, color = "black")+
  geom_point(data = sgmrna,aes(x=end_pos,y = segment),size=1.5, color = "black")+
  xlab("Nucleotide position")+
  theme_classic()+
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



