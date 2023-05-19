library(tidyverse)
library(rio)
library(compiler)
library(here)

load("prefilters.RData")
prefilter <- merged
names(prefilter) <- gsub("user_","pref_",names(prefilter))

window <- c(0:20)

filter.window <- function(pre){
  pre <- bind_rows(pre)
  pre <- distinct(pre, start_pos, end_pos, .keep_all = T)
  window <- 0:20
  filter.list <- list()
  filter.count <- NA
  for(i in seq_along(window)){
    filter.list[[i]] <- pre %>%
      filter(!start_pos %in% c((65-i):(65+i)) | !end_pos %in% c((21552-i):(21552+i)))%>% #S
      filter(!start_pos %in% c((69-i):(69+i)) | !end_pos %in% c((25385-i):(25385+i)))%>% #ORF3a
      filter(!start_pos %in% c((69-i):(69+i)) | !end_pos %in% c((26237-i):(26237+i)))%>% #E
      filter(!start_pos %in% c((65-i):(65+i)) | !end_pos %in% c((26469-i):(26469+i)))%>% #M
      filter(!start_pos %in% c((69-i):(69+i)) | !end_pos %in% c((27041-i):(27041+i)))%>% #ORF6
      filter(!start_pos %in% c((69-i):(69+i)) | !end_pos %in% c((27388-i):(27388+i)))%>% #ORF7a
      filter(!start_pos %in% c((65-i):(65+i)) | !end_pos %in% c((27644-i):(27644+i)))%>% #ORF7b
      filter(!start_pos %in% c((65-i):(65+i)) | !end_pos %in% c((27884-i):(27884+i)))%>% #ORF8
      filter(!start_pos %in% c((65-i):(65+i)) | !end_pos %in% c((28256-i):(28256+i)))%>% #N
      filter(!start_pos %in% c((65-i):(65+i)) | !end_pos %in% c((29530-i):(29530+i))) #ORF10
    filter.count[[i]] <- length(filter.list[[i]]$normalized_reads)
  }
  filter.count <- unlist(filter.count)
  return(filter.count)
} #this function is extremely slow but it does work
filter.window <- cmpfun(filter.window)
filtered <- lapply(prefilter, filter.window)

filter.table <- bind_cols(window,filtered)
colnames(filter.table)[1] <- "window"
totals <- rowSums(filter.table[,-1])
plot(window, totals, type = "l") #combined filtering for all participants
filter.gathered <- filter.table %>%
  gather(-window, value = "Count", key = "Participant")

#individual filtering results (plotted on one graph)
allblack <- rep("black",21)
ggplot(data = filter.gathered, aes(x = window, y = Count, color = Participant))+
  geom_line()+
  xlab("Filter window")+
  ylab("Unique junctions")+
  scale_x_continuous(limits = c(0,20),breaks = c(0,5,10,15,20))+
  scale_color_manual(values = allblack)+
  geom_vline(xintercept = 15, linetype = "dashed")+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        legend.position = "none")
ggsave("figs/filter_window.png")




