####################################
# 1_data_process
####################################

## Load packages
	intstall.packages("pacman")
	pacman::p_load(tidyverse, readxl)

## Load data	
	              activ_files <- list.files("./data/Activity/")
	       data_activity_list <- lapply(activ_files, function(x) read_xlsx(paste0("./data/Activity/", x))[-c(1:3),])
	names(data_activity_list) <- activ_files

## Merge
	data_activity <- dplyr::bind_rows(data_activity_list, .id = "id")



# At the end of processing you write the file
	write_csv(data_activity, "./output/data/data_activity.csv")

