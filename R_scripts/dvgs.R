library(tidyverse)
library(rio)
library(gdata)
library(here)

#read in ViReMa output
dir.list <- list.dirs("users")[-1]
ext.files <- list()
for(i in seq_along(dir.list)){
  ext.files[[i]] <- paste0(dir.list[[i]],"/", list.files(dir.list[[i]]))
}
for(i in seq_along(ext.files)){
  ext.files[[i]] <- suppressMessages(lapply(ext.files[[i]],read_tsv,col_names=F,comment = "\t"))
}

#split tables by forward and reverse strand junctions
split.library <- function(DItable){
  DItable <- as_tibble(DItable)
  DItable <- DItable %>%
    rename("DI"=X1)
  for(i in seq_along(DItable$DI)){
    if("@NewLibrary: MN985325.1_to_MN985325.1" %in% DItable$DI){
      if(DItable$DI[[i]] == "@NewLibrary: MN985325.1_to_MN985325.1"){
        forwardstart <- i
      }
    } else {
      forwardstart <- NA
    }
  } 
  for(i in seq_along(DItable$DI)){
    if("@NewLibrary: MN985325.1_RevStrand_to_MN985325.1_RevStrand" %in% DItable$DI){
      if(DItable$DI[[i]] == "@NewLibrary: MN985325.1_RevStrand_to_MN985325.1_RevStrand"){
        reversestart <- i
      }
    } else {
      reversestart <- NA
    }
  } 
  for(i in seq_along(DItable$DI)){
    if("@NewLibrary: MN985325.1_to_MN985325.1_RevStrand" %in% DItable$DI){
      if(DItable$DI[[i]] == "@NewLibrary: MN985325.1_to_MN985325.1_RevStrand"){
        forrevstart <- i
      }
    } else {
      forrevstart <- NA
    }
  } 
  for(i in seq_along(DItable$DI)){
    if("@NewLibrary: MN985325.1_RevStrand_to_MN985325.1" %in% DItable$DI){
      if(DItable$DI[[i]] == "@NewLibrary: MN985325.1_RevStrand_to_MN985325.1"){
        revforstart <- i
      }
    } else {
      revforstart <- NA
    }
  } 
  starts <- c(forwardstart,reversestart,forrevstart,revforstart)
  names(starts) <- c("forward","reverse","forrev","revfor")
  starts <- sort(starts)
  table <- list()
  for(i in seq_along(starts)){
    if(i < length(starts)){
      table[[i]] <- DItable[(starts[[i]]+1):(starts[[i+1]]-2),]
    }
    if(i == length(starts)){
      table[[i]] <- DItable[(starts[[i]]+1):(length(DItable$DI)-1),]
    }
  }
  names(table) <- names(starts)
  table <- table[c("forward", "reverse")]
  return(table)
}

tabs <- list()
for(i in seq_along(ext.files)){
  tabs[[i]] <- lapply(ext.files[[i]],split.library)
}

#add temporary placeholder DVG for empty tables
placeholder <- import("placeholder.xlsx")
placeholder <- as_tibble(placeholder)

hold.place <- function(tab){
  for(i in seq_along(tab)){
    if(length(tab[[i]]$forward) < 1){
      tab[[i]]$forward <- placeholder
    }
    if(length(tab[[i]]$reverse) < 1){
      tab[[i]]$reverse <- placeholder
    }
    tab[[i]] <- discard(tab[[i]], ~ is.null(.x) )
  }
  return(tab)
}
tabs <- lapply(tabs, hold.place)

#collate forward and reverse junctions into columns within one table
table.collapse <- function(tab){
  whole <- list()
  for(i in seq_along(tab)){
    whole[[i]] <- gdata::cbindX(tab[[i]][[1]],tab[[i]][[2]])
    colnames(whole[[i]]) <- names(tab[[i]])
  }
  return(whole)
}
whole <- lapply(tabs, table.collapse)

#parse tables
forward.parse <- function(whole.tab){
  forward <- as_tibble(whole.tab$forward)
  forward <- forward %>%
    mutate("start_pos"=NA)%>%
    mutate("end_pos"=NA)%>%
    mutate("min"=NA)%>%
    mutate("max"=NA)%>%
    mutate("reads"=NA)
  for(i in seq_along(forward$value)){
    forward$start_pos[[i]] <- strsplit(forward$value[[i]],"_")[[1]][1]
    forward$end_pos[[i]] <- strsplit(forward$value[[i]],"_")[[1]][3]
    forward$reads[[i]] <- strsplit(forward$value[[i]],"_")[[1]][5]
  }
  forward$start_pos <- as.numeric(forward$start_pos)
  forward$end_pos <- as.numeric(forward$end_pos)
  for(i in seq_along(forward$value)){
    forward$min[[i]] <- min(forward$start_pos[[i]],forward$end_pos[[i]])
    forward$max[[i]] <- max(forward$start_pos[[i]],forward$end_pos[[i]])
  }
  
  forward <- forward %>%
    select(min,max,reads)%>%
    rename("start_pos"=min)%>%
    rename("end_pos"=max)
  return(forward)
}
for.parsed <- list()
for(i in seq_along(whole)){
  for.parsed[[i]] <- lapply(whole[[i]],forward.parse)
}

reverse.parse <- function(whole.tab){
  reverse <- as_tibble(whole.tab$reverse)
  reverse <- reverse %>%
    mutate("start_pos"=NA)%>%
    mutate("end_pos"=NA)%>%
    mutate("min"=NA)%>%
    mutate("max"=NA)%>%
    mutate("reads"=NA)
  for(i in seq_along(reverse$value)){
    reverse$start_pos[[i]] <- strsplit(reverse$value[[i]],"_")[[1]][1]
    reverse$end_pos[[i]] <- strsplit(reverse$value[[i]],"_")[[1]][3]
    reverse$reads[[i]] <- strsplit(reverse$value[[i]],"_")[[1]][5]
  }
  reverse$start_pos <- as.numeric(reverse$start_pos)
  reverse$end_pos <- as.numeric(reverse$end_pos)
  for(i in seq_along(reverse$value)){
    reverse$min[[i]] <- min(reverse$start_pos[[i]],reverse$end_pos[[i]])
    reverse$max[[i]] <- max(reverse$start_pos[[i]],reverse$end_pos[[i]])
  }
  
  reverse <- reverse %>%
    select(min,max,reads)%>%
    rename("start_pos"=min)%>%
    rename("end_pos"=max)
  return(reverse)
}
rev.parsed <- list()
for(i in seq_along(whole)){
  rev.parsed[[i]] <- lapply(whole[[i]],reverse.parse)
}

#merge forward and reverse tables
table.merge <- function(forward,reverse){
  together <- full_join(forward,reverse,by=c("start_pos","end_pos"))
  together <- together %>%
    rename("reads.for"=reads.x)%>%
    rename("reads.rev"=reads.y)
  together$reads.for <- as.numeric(together$reads.for)
  together$reads.rev <- as.numeric(together$reads.rev)
  together <- together[rowSums(is.na(together)) != ncol(together),]
  for(i in seq_along(together$reads.for)){
    if(is.na(together$reads.for[[i]])){
      together$reads.for[[i]] <- 0
    }
    if(is.na(together$reads.rev[[i]])){
      together$reads.rev[[i]] <- 0
    }
  }
  together <- together %>%
    mutate("total_reads"= (reads.for + reads.rev))%>%
    mutate("distance"=end_pos-start_pos)
  return(together)
}
merged <- list()
for(i in seq_along(for.parsed)){
  merged[[i]] <- map2(for.parsed[[i]], rev.parsed[[i]], table.merge)
}

aggregate.reads <- function(merged.tab){
  for(i in seq_along(merged.tab)){
    merged.tab[[i]] <- aggregate(total_reads ~ start_pos + end_pos + distance, data = merged.tab[[i]],
                           FUN = "sum")
  }
  return(merged.tab)
}
merged <- lapply(merged, aggregate.reads)

for(i in seq_along(merged)){
  names(merged[[i]]) <- gsub("_extracted.txt","",(list.files(dir.list[[i]])))
}

names(merged) <- gsub("users/","",dir.list)

#normalize to reads and bare-minimum filter junctions based on length and position
reads <- lapply(paste0("alignments/",list.files("alignments/")), import)
names(reads) <- gsub(".csv","",list.files("alignments/"))

merged.clean <- function(merged.tab, align.tab){
  for(i in seq_along(merged.tab)){
    merged.tab[[i]] <- merged.tab[[i]] %>%
      filter(distance >= 20)%>%
      filter(start_pos < 29870)%>%
      filter(end_pos < 29870)%>%
      mutate("normalized_reads"= (total_reads*1000000)/ ((align.tab$aligned[[i]]) + total_reads)) %>%
      select(-c(total_reads,distance)) %>%
      filter(normalized_reads >= 10)
  }
  return(merged.tab)
}
merged <- map2(merged,reads,merged.clean)

save(merged, file = "prefilters.RData")

#filter junctions within 15 nts of TRS
table.filter <- function(merged.tab){
  merged.tab <- merged.tab %>%
    filter(!start_pos %in% c(55:85) | !end_pos %in% c(21548:21578))%>% #S
    filter(!start_pos %in% c(55:85) | !end_pos %in% c(25378:25408))%>% #ORF3a
    filter(!start_pos %in% c(55:85) | !end_pos %in% c(26230:26260))%>% #E
    filter(!start_pos %in% c(55:85) | !end_pos %in% c(26508:26538))%>% #M
    filter(!start_pos %in% c(55:85) | !end_pos %in% c(27187:27217))%>% #ORF6
    filter(!start_pos %in% c(55:85) | !end_pos %in% c(27379:27409))%>% #ORF7a
    filter(!start_pos %in% c(55:85) | !end_pos %in% c(27741:27771))%>% #ORF7b
    filter(!start_pos %in% c(55:85) | !end_pos %in% c(27879:27909))%>% #ORF8
    filter(!start_pos %in% c(55:85) | !end_pos %in% c(28259:28289))%>% #N
    filter(!start_pos %in% c(55:85) | !end_pos %in% c(29543:29573))%>% #ORF10
    return(merged.tab)
} #sgmrna filtering of 15nts 
together <- list()
for(i in seq_along(merged)){
  together[[i]] <- lapply(merged[[i]], table.filter)
}

#create junction tracking tables for each user
user.di <- lapply(together, purrr::reduce, full_join, by =c("start_pos","end_pos"))

di.vec <- list()
for(i in seq_along(user.di)){
  di.vec[[i]] <- 1:((dim(user.di[[i]])[[2]])-2)
  colnames(user.di[[i]]) <- c("start_pos","end_pos",paste0("normalized_reads_",di.vec[[i]]))
}
junction.di <- user.di #will use for a later analysis

for(i in seq_along(user.di)){
  user.di[[i]] <- user.di[[i]] %>% 
    mutate("DI" = paste0(start_pos,"_",end_pos))%>%
    relocate(DI,.after = end_pos)%>%
    select(-c(start_pos,end_pos))
}
zero.di <- user.di #will use in a second
user.di <- lapply(user.di,t)
names(user.di) <- names(merged)

for(i in seq_along(user.di)){
  filenames <- paste0("dvgs/",gsub("user_","dvg_",names(user.di)),".csv")
  write.csv(user.di[[i]],filenames[[i]])
}

#calculate total normalized reads for each user
for(i in seq_along(zero.di)){
  zero.di[[i]] <- replace(zero.di[[i]],is.na(zero.di[[i]]),0)
}

sum.table <- lapply(zero.di, select, -DI)
total.norm.reads <- lapply(sum.table, sum)

names(total.norm.reads) <- names(user.di)
save(total.norm.reads, file = "total_norm_reads.RData") 

#junction stats per day of infection
daily.sums <- lapply(sum.table, colSums)
names(daily.sums) <- names(user.di)
save(daily.sums, file = "daily_sums.RData") #daily total normalized reads for all junctions

count.di <- sum.table
diversity.count <- function(di.table){
  daily.counts <- list()
  for(i in seq_along(di.table[1,])){
    daily.counts[[i]] <- length(which(di.table[[i]] > 0))
  }
  daily.counts <- unlist(daily.counts)
  return(daily.counts)
}
daily.counts <- lapply(count.di, diversity.count)
names(daily.counts) <- names(user.di)
save(daily.counts, file = "daily_counts.RData") #daily unique junction counts

#deletion accumulation analysis
junction.sums <- lapply(sum.table, rowSums)

#deletion enrichment analysis
for(i in seq_along(zero.di)){
  rownames(zero.di[[i]]) <- zero.di[[i]]$DI
}
zero.di <- lapply(zero.di, select, -DI)
junction.enrich <- zero.di

enrich.calc <- function(di.table,sums){
  for(i in seq_along(di.table[,1])){
    for(j in seq_along(di.table[1,])){
      di.table[i,j] <- di.table[i,j]/(sums[[j]])
    }
  }
  return(di.table)
}
junction.enrich <- map2(junction.enrich, daily.sums, enrich.calc)
junction.enrich <- lapply(junction.enrich, as_tibble)

for(i in seq_along(junction.enrich)){
  junction.enrich[[i]] <- replace(junction.enrich[[i]],is.na(junction.enrich[[i]]),0)
}

junction.enrich <- lapply(junction.enrich, mutate, "delta_enrich" = 0)

delta.enrich.calc <- function(di.table){
  for(i in seq_along(di.table$delta_enrich)){
    di.table$delta_enrich[[i]] <- max(di.table[i,]) - di.table$normalized_reads_1[[i]]
  }
  return(di.table)
}
junction.enrich <- lapply(junction.enrich, delta.enrich.calc)

#collate junction stats
junction.di <- lapply(junction.di, select, start_pos, end_pos)
junction.di <- lapply(junction.di, mutate, "distance" = end_pos - start_pos)

junction.stats <- list()
for(i in seq_along(junction.di)){
  junction.stats[[i]] <- junction.di[[i]] %>%
    mutate("accumulation" = junction.sums[[i]])
  junction.stats[[i]] <- junction.stats[[i]] %>%
    mutate("delta_enrich" = junction.enrich[[i]]$delta_enrich)
}

names(junction.stats) <- names(user.di)
save(junction.stats, file = "junction_stats.RData")




