library(tidyverse)
library(rio)
library(here)

files <- paste0("dvgs/",list.files("dvgs/"))
dvg.tables <- lapply(files, read.csv)

for(i in seq_along(dvg.tables)){
  colnames(dvg.tables[[i]]) <- dvg.tables[[i]][1,]
  dvg.tables[[i]] <- dvg.tables[[i]][-1,]
}

shared.only <- function(dvg){
  for(i in seq_along(colnames(dvg))){
    if(length(which(!is.na(dvg[,i]))) <= 1){
      dvg[,i] <- NA
    }
  }
  not_all_na <- function(x) any(!is.na(x))
  dvg <- dvg %>%
    select(where(not_all_na))
  return(dvg)
}

shared.dvgs <- lapply(dvg.tables, shared.only)
names(shared.dvgs) <- gsub("dvgs/dvg_","user_",files)
names(shared.dvgs) <- gsub(".csv","",names(shared.dvgs))
save(shared.dvgs, file = "shared_dvgs.RData")

shared.names <- list()
for(i in seq_along(shared.dvgs)){
  shared.names[[i]] <- colnames(shared.dvgs[[i]])
}
shared.names <- unlist(shared.names)
shared.names <- shared.names[-which(shared.names == "DI")]
shared.table <- table(shared.names)
shared.table <- as.data.frame(shared.table)
colnames(shared.table) <- c("junction", "Frequency")
write.csv(shared.table, file = "shared_table.csv")

ggplot(shared.table, aes(x = Frequency))+
  geom_bar()+
  xlab("Number of participants \n with junction")+
  ylab("Junction count")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 19, face = "bold"),
        axis.title.x = element_text(size = 22, face = "bold"),
        axis.text.y = element_text(size = 19, face = "bold"),
        axis.title.y = element_text(size = 22, face = "bold"))
ggsave("figs/junction_histogram_shared.png")


