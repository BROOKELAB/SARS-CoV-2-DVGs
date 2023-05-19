library(tidyverse)
library(rio)
library(gridExtra)
library(here)

load("daily_sums.RData")
dailies <- daily.sums
names(dailies) <- gsub("user_","daily_",names(dailies))

dailies <- lapply(dailies,as.vector)
ct <- import("naive_ct.xlsx") %>%
  select(-c(sequenced,SNP_count))%>%
  group_by(user_id)%>%
  group_split(user_id)

ct <- lapply(ct, as.data.frame)
dailies <- lapply(dailies, as.data.frame)

for(i in seq_along(ct)){
  ct[[i]] <- bind_cols(ct[[i]],dailies[[i]])
  colnames(ct[[i]])[4] <- "junction_reads"
}

ctplots <- list()
for(i in seq_along(ct)){
  ctplots[[i]] <- ggplot(data = ct[[i]], aes(x = ct, y = junction_reads))+
    geom_point()+
    xlab("Sample Ct")+
    ylab("Total DVG reads")+
    ggtitle(gsub("daily_","",names(dailies)[[i]]))+
    theme_bw()
}

all <- grid.arrange(ctplots[[1]],ctplots[[2]],ctplots[[3]],ctplots[[4]],
                    ctplots[[5]],ctplots[[6]],ctplots[[7]],ctplots[[8]],
                    ctplots[[9]],ctplots[[10]],ctplots[[11]],ctplots[[12]],
                    ctplots[[13]],ctplots[[14]],ctplots[[15]],ctplots[[16]],
                    ctplots[[17]],ctplots[[18]],ctplots[[19]],ctplots[[20]],
                    nrow = 5)
cor.ct <- list()
for(i in seq_along(ct)){
  cor.ct[[i]] <- cor.test(ct[[i]]$ct, ct[[i]]$junction_reads, method  = "k")
}
names(cor.ct) <- gsub("daily_","",names(dailies))

p.ct <- list()
for(i in seq_along(cor.ct)){
  p.ct[[i]] <- cor.ct[[i]]$p.value
}
p.ct <- unlist(p.ct)
which(p.ct < 0.05)

#2  `432870` (p = 0.01591, tau = 0.6479516)
#15 `451152` (p = 0.01414, tau = -0.7142857)



