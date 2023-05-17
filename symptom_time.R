library(tidyverse)
library(rio)
library(here)

load("daily_sums.RData")
dvg.dailies <- daily.sums
names(dvg.dailies) <- gsub("user_", "daily_",names(dvg.dailies))
dvg.names <- gsub("daily_","",names(dvg.dailies))

daily.matrix <- function(daily){
  mat <- as.data.frame(matrix(nrow = length(daily), ncol = 2))
  mat[,1] <- names(daily)
  mat[,2] <- daily
  return(mat)
}
daily.mats <- lapply(dvg.dailies, daily.matrix)

all.daily <- bind_rows(daily.mats)

user.info <- import("all_dvg_user_info.xlsx")
user.info$sample_date <- as.Date(user.info$sample_date)
symptoms <- import("user_symptoms.csv") %>%
  rename(sample_date = labs_vdl_date)%>%
  select(-V1)
symptoms$sample_date <- as.Date(symptoms$sample_date, "%m/%d/%Y")

all.info <- inner_join(user.info, symptoms, by = c("user_id","sample_date"))
all.info <- all.info %>%
  select(-c(sample_barcode,status,ct))%>%
  mutate("dvg_reads" = all.daily$V2)%>%
  relocate(dvg_reads, .after = day_of_infection)
all.info <- all.info %>%
  rowwise()%>%
  mutate("total_symptoms" = sum(c(scratchy_throat,sore_throat,cough,runny_nose,fever_chills,high_temp,
                                  muscle_aches,nausea_vomiting_diarrhea,shortness_of_breath,
                                  unable_to_taste_or_smell,red_eyes)))%>%
  mutate("severe_symptoms" = sum(c(fever_chills,high_temp,muscle_aches,nausea_vomiting_diarrhea,
                                   shortness_of_breath,unable_to_taste_or_smell)))
all.info <- all.info %>%
  ungroup() %>%
  group_split(day_of_infection)

rs <- list()
ps <- list()
rs.severe <- list()
ps.severe <- list()
for(i in seq_along(all.info[1:11])){ #only those dates with > 2 observations
  rs[[i]] <- cor.test(all.info[[i]]$dvg_reads, all.info[[i]]$total_symptoms, 
                      method = "pearson")$estimate #pearson since spearman cannot compute with ties
  ps[[i]] <- cor.test(all.info[[i]]$dvg_reads, all.info[[i]]$total_symptoms, 
                   method = "pearson")$p.value
  rs.severe[[i]] <- cor.test(all.info[[i]]$dvg_reads, all.info[[i]]$severe_symptoms, 
                             method = "pearson")$estimate
  ps.severe[[i]] <- cor.test(all.info[[i]]$dvg_reads, all.info[[i]]$severe_symptoms, 
                             method = "pearson")$p.value
}
rs <- unlist(rs)
ps <- unlist(ps)
rs.severe <- unlist(rs.severe)
ps.severe <- unlist(ps.severe)

cor.mat <- as.data.frame(matrix(data = c(c(1:11),rs,ps), ncol = 3,nrow = 11))
colnames(cor.mat) <- c("Day","R", "P value")
cor.mat <- cor.mat %>%
  gather(key = "Statistic", value = "Value",-Day)

severe.mat <- as.data.frame(matrix(data = c(c(1:11),rs.severe,ps.severe), ncol = 3,nrow = 11))
colnames(severe.mat) <- c("Day","R", "P value")
severe.mat <- severe.mat %>%
  gather(key = "Statistic", value = "Value",-Day)

ggplot(data = cor.mat, aes(x = Day, y = Value, color = Statistic))+
  geom_point(cex = 3)+
  xlab("Day post-enrollment")+
  ylab("")+
  ggtitle("All symptoms")+
  scale_x_discrete(limits = factor(c(1:11)))+
  scale_color_manual(values = c("red","black"))+
  geom_hline(yintercept = 0.05, linetype = "dashed")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19))
ggsave("figs/symptom_cors.png")  

ggplot(data = severe.mat, aes(x = Day, y = Value, color = Statistic))+
  geom_point(cex = 3)+
  xlab("Day post-enrollment")+
  ylab("")+
  ggtitle("Severe symptoms")+
  scale_x_discrete(limits = factor(c(1:11)))+
  scale_color_manual(values = c("red","black"))+
  geom_hline(yintercept = 0.05, linetype = "dashed")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19))
ggsave("figs/symptom_cors_severe.png") 


