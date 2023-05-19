library(tidyverse)
library(rio)
library(gdata)
library(here)

load("daily_counts.RData")
counts.list <- daily.counts
names(counts.list) <- gsub("user_", "count_", names(counts.list))

for(i in seq_along(counts.list)){
  length(counts.list[[i]]) <- 10
}

all.counts <- bind_cols(counts.list)
all.counts <- all.counts %>%
  mutate("time" = c(1:10))%>%
  gather(key = "ID", value = "Counts", -time)

all.counts$time <- as.character(all.counts$time)
ggplot(all.counts, aes(x = time, y = Counts)) +
  geom_point(cex = 3)+
  xlab("Time")+
  ylab("Unique junctions")+
  scale_x_discrete(limits = factor(c(1:10)))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 19),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 19),
        axis.title.y = element_text(size = 22))
ggsave("figs/counts_vs_time.png")

all.counts$time <- as.numeric(all.counts$time)
cor.test(all.counts$time, all.counts$Counts, method = "s")
#rho = -0.07390945  #p-value = 0.3639

