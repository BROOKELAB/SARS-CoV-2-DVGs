library(tidyverse)
library(rio)
library(here)

load("total_norm_reads.RData")
dvg.sums <- total.norm.reads
names(dvg.sums) <- gsub("user_","sum_",names(dvg.sums))

dvg.sums <- unlist(dvg.sums)
sample.no <- c("432686"=6,"432870"=9,"433227"=7,"435786"=6,"435805"=9,
               "438577"=6,"442978"=5,"444332"=8,"444446"=9,"444633"=5,
               "445602"=10,"449650"=8,"450241"=6,"450348"=9,"451152"=8,
               "451709"=9,"453058"=7,"459597"=7,"471467"=10,"471588"=9)
dvg.sums.norm <- dvg.sums/sample.no

symptoms <- read.csv("cleaned_data_061021.csv")
symptoms <- symptoms %>%
  filter(user_id == "432686" | user_id == "432870" | user_id == "433227" | 
           user_id == "435786" | user_id == "435805" | user_id == "438577" |
           user_id == "442978" | user_id == "444332" | user_id == "444446" |
           user_id == "444633" | user_id == "445602" | user_id == "449650" |
           user_id == "450241" | user_id == "450348" | user_id == "451152" |
           user_id == "451709" | user_id == "453058" | user_id == "459597" |
           user_id == "471467" | user_id == "471588") %>%
  select(user_id, labs_vdl_date, scratchy_throat, sore_throat, cough, runny_nose, fever_chills, high_temp,
         muscle_aches, nausea_vomiting_diarrhea, shortness_of_breath, 
         unable_to_taste_or_smell, red_eyes)

for(i in seq_along(symptoms[,1])){
  for(j in seq_along(symptoms[1,])){
    if(is.na(symptoms[i,j])){
      symptoms[i,j] <- 0
    }
  }
}

write.csv(symptoms,"user_symptoms.csv")

symptoms <- symptoms %>%
  rowwise() %>%
  mutate("total_symptoms" = sum(c(scratchy_throat, sore_throat, cough, runny_nose, fever_chills, high_temp,
                                muscle_aches, nausea_vomiting_diarrhea, shortness_of_breath, 
                                unable_to_taste_or_smell, red_eyes)))
user.symptoms <- symptoms %>%
  ungroup() %>%
  group_by(user_id) %>%
  group_split(user_id)
names(user.symptoms) <- levels(factor(symptoms$user_id))
user.total.symptoms <- list()
for(i in seq_along(user.symptoms)){
  user.total.symptoms[[i]] <- sum(user.symptoms[[i]]$total_symptoms)
}
user.total.symptoms <- unlist(user.total.symptoms)
corline <- lm(user.total.symptoms ~ dvg.sums.norm) 
summary(corline)
#p = 0.1586

symptom.dvgs <- as_tibble(cbind(dvg.sums.norm, user.total.symptoms))
ggplot(data = symptom.dvgs, aes(x = dvg.sums.norm, y = user.total.symptoms))+
  geom_point(cex=3)+
  xlab("Total DVG reads")+
  ylab("Symptom count")+
  ggtitle("DVG load vs. symptom load")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 22))
ggsave("figs/symptoms_dvgs.png")

#for just anosmia
anosmia <- symptoms %>%
  select(user_id,labs_vdl_date,unable_to_taste_or_smell)
user.anosmia <- anosmia %>%
  group_by(user_id)%>%
  group_split(user_id)
names(user.anosmia) <- levels(factor(anosmia$user_id))
user.total.anosmia <- list()
for(i in seq_along(user.anosmia)){
  user.total.anosmia[[i]] <- sum(user.anosmia[[i]]$unable_to_taste_or_smell)
}
user.total.anosmia <- unlist(user.total.anosmia)
anosmia.cor <- lm(user.total.anosmia ~ dvg.sums.norm) 
summary(anosmia.cor)
#p = 0.2305

anosmia.dvgs <- as_tibble(cbind(dvg.sums.norm, user.total.anosmia))
ggplot(data = anosmia.dvgs, aes(x = dvg.sums.norm, y = user.total.anosmia))+
  geom_point(cex=3)+
  xlab("Total DVG reads")+
  ylab("Symptom count")+
  ggtitle("DVG load vs. anosmia")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 22))
ggsave("figs/anosmia_dvgs.png")

#for just muscle aches
muscle <- symptoms %>%
  select(user_id,labs_vdl_date,muscle_aches)
user.muscle <- muscle %>%
  group_by(user_id)%>%
  group_split(user_id)
names(user.muscle) <- levels(factor(muscle$user_id))
user.total.muscle <- list()
for(i in seq_along(user.muscle)){
  user.total.muscle[[i]] <- sum(user.muscle[[i]]$muscle_aches)
}
user.total.muscle <- unlist(user.total.muscle)
muscle.cor <- lm(user.total.muscle ~ dvg.sums.norm) 
summary(muscle.cor)
#p = 0.8445

muscle.dvgs <- as_tibble(cbind(dvg.sums.norm, user.total.muscle))
ggplot(data = muscle.dvgs, aes(x = dvg.sums.norm, y = user.total.muscle))+
  geom_point(cex=3)+
  xlab("Total DVG reads")+
  ylab("Symptom count")+
  ggtitle("DVG load vs. muscle aches")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 22))
ggsave("figs/muscleache_dvgs.png")

breath <- symptoms %>%
  select(user_id,labs_vdl_date,shortness_of_breath)
user.breath <- breath %>%
  group_by(user_id)%>%
  group_split(user_id)
names(user.breath) <- levels(factor(breath$user_id))
user.total.breath <- list()
for(i in seq_along(user.breath)){
  user.total.breath[[i]] <- sum(user.breath[[i]]$shortness_of_breath)
}
user.total.breath <- unlist(user.total.breath)
breath.cor <- lm(user.total.breath ~ dvg.sums.norm) 
summary(breath.cor)
#p = 0.3892

breath.dvgs <- as_tibble(cbind(dvg.sums.norm, user.total.breath))
ggplot(data = breath.dvgs, aes(x = dvg.sums.norm, y = user.total.breath))+
  geom_point(cex=3)+
  xlab("Total DVG reads")+
  ylab("Symptom count")+
  ggtitle("DVG load vs. shortness of breath")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 22))
ggsave("figs/breath_dvgs.png")

