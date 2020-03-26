# To be run from data directory where hogwash results are saved
source("../lib/resource_usage_lib.R") # RC wrote this library
library(tidyverse)
project_names <- read_tsv("../data/project_names_and_colors.tsv", 
                          col_names = TRUE)
# Remove projects not included in final prewas draft
project_names <- 
  project_names[!project_names$Project %in%
                  c("Cd_1", "Cd_2", "Kp_2", "Ec_1", "Ece_1", "Lp_1", "Sa_3"), ,
                drop = FALSE]

num_project <- nrow(project_names)
num_row <- num_project * 7
# Quintupled, so we can have: 
# 1. Major allele, no tree given, 1 core
# 2. With AR, no tree given, 1 core
# 3. With AR, no tree given, 4 core
# 4. With AR, no tree given, 10 core
# 5. With AR, given tree, 1 core
# 6. With AR, given tree, 4 core
# 7. With AR, given tree, 10 core

resource_df <- as.data.frame(matrix(NA, nrow = num_row, ncol = 7))
colnames(resource_df) <- c("Time (Hours)",
                           "Memory (GB)", 
                           "Project", 
                           "Method",  # Major Allele or Ancestral Reconstruction
                           "Tree", # Tree Provided or No Tree Provided
                           "Cores (#)", 
                           "Variants (#)")
row_num <- 1
for (i in 1:num_project) {
  project_name <- project_names$Project[i]
  
  # 1. Major allele, no tree given, 1 core
  major_out_file_name <- # ex: Lc_1_major_no_tr_given_1
    paste0("../data/", project_name, "_major_no_tr_given_1.out") 
  if (file.exists(major_out_file_name)) {
    out_file_info <- readLines(con = major_out_file_name, n = 7)
    job_id <- out_file_info[2]
    num_var <- out_file_info[7]
    time_and_mem <- GetResourceUsage(jobId = job_id) # Returns vector of length 2
    time <- time_and_mem[1]
    mem <- time_and_mem[2]
    num_var <- str_squish(num_var)
    num_var <- as.numeric(gsub("variant count: ", "", num_var))
    resource_df[row_num, ] <- 
      c(time, mem, project_name, "Major Allele", "No Tree Provided", 1, num_var)
  }
 
  # 2. With AR, no tree given, 1 core
  AR_no_tr_1_out_file_name <- # ex: Lc_1_w_AR_no_tr_given_1.out
    paste0("../data/", project_name, "_w_AR_no_tr_given_1.out")
  if (file.exists(AR_no_tr_1_out_file_name)) {
    out_file_info <- readLines(con =  AR_no_tr_1_out_file_name, n = 7)
    job_id <- out_file_info[2]
    time_and_mem <- GetResourceUsage(jobId = job_id) # Returns vector of length 2
    time <- time_and_mem[1]
    mem <- time_and_mem[2]
    num_var <- str_squish(num_var)
    num_var <- as.numeric(gsub("variant count: ", "", num_var))
    resource_df[row_num + 1, ] <-  
      c(time, mem, project_name, "Ancestral Reconstruction", "No Tree Provided", 1, num_var)
  }
  
  # 3. With AR, no tree given, 4 core
  AR_no_tr_4_out_file_name <- # ex: Lc_1_w_AR_no_tr_given_4.out
    paste0("../data/", project_name, "_w_AR_no_tr_given_4.out")
  if (file.exists(AR_no_tr_4_out_file_name)) {
    out_file_info <- readLines(con =  AR_no_tr_4_out_file_name, n = 7)
    job_id <- out_file_info[2]
    time_and_mem <- GetResourceUsage(jobId = job_id) # Returns vector of length 2
    time <- time_and_mem[1]
    mem <- time_and_mem[2]
    num_var <- str_squish(num_var)
    num_var <- as.numeric(gsub("variant count: ", "", num_var))
    resource_df[row_num + 2, ] <-  
      c(time, mem, project_name, "Ancestral Reconstruction", "No Tree Provided", 4, num_var)
  }
  
  # 4. With AR, no tree given, 10 core
  AR_no_tr_4_out_file_name <- # ex: Lc_1_w_AR_no_tr_given_4.out
    paste0("../data/", project_name, "_w_AR_no_tr_given_10.out")
  if (file.exists(AR_no_tr_4_out_file_name)) {
    out_file_info <- readLines(con =  AR_no_tr_4_out_file_name, n = 7)
    job_id <- out_file_info[2]
    time_and_mem <- GetResourceUsage(jobId = job_id) # Returns vector of length 2
    time <- time_and_mem[1]
    mem <- time_and_mem[2]
    num_var <- str_squish(num_var)
    num_var <- as.numeric(gsub("variant count: ", "", num_var))
    resource_df[row_num + 3, ] <-  
      c(time, mem, project_name, "Ancestral Reconstruction", "No Tree Provided", 10, num_var)
  }
  
  # 5. With AR, given tree, 1 core
  AR_tr_1_out_file_name <- # ex: Lc_1_w_AR_w_tr_1.out
    paste0("../data/", project_name, "_w_AR_w_tr_1.out")
  if (file.exists(AR_tr_1_out_file_name)) {
    out_file_info <- readLines(con =  AR_tr_1_out_file_name, n = 7)
    job_id <- out_file_info[2]
    time_and_mem <- GetResourceUsage(jobId = job_id) # Returns vector of length 2
    time <- time_and_mem[1]
    mem <- time_and_mem[2]
    num_var <- str_squish(num_var)
    num_var <- as.numeric(gsub("variant count: ", "", num_var))
    resource_df[row_num + 4, ] <-  
      c(time, mem, project_name, "Ancestral Reconstruction", "Tree Provided", 1, num_var)
  }
  
  # 6. With AR, given tree, 1 core
  AR_tr_4_out_file_name <- # ex: Lc_1_w_AR_w_tr_4.out
    paste0("../data/", project_name, "_w_AR_w_tr_4.out")
  if (file.exists(AR_tr_4_out_file_name)) {
    out_file_info <- readLines(con =  AR_tr_4_out_file_name, n = 7)
    job_id <- out_file_info[2]
    time_and_mem <- GetResourceUsage(jobId = job_id) # Returns vector of length 2
    time <- time_and_mem[1]
    mem <- time_and_mem[2]
    num_var <- str_squish(num_var)
    num_var <- as.numeric(gsub("variant count: ", "", num_var))
    resource_df[row_num + 5, ] <-  
      c(time, mem, project_name, "Ancestral Reconstruction", "Tree Provided", 4, num_var)
  }
  
  # 7. With AR, given tree, 1 core
  AR_tr_4_out_file_name <- # ex: Lc_1_w_AR_w_tr_4.out
    paste0("../data/", project_name, "_w_AR_w_tr_10.out")
  if (file.exists(AR_tr_4_out_file_name)) {
    out_file_info <- readLines(con =  AR_tr_4_out_file_name, n = 7)
    job_id <- out_file_info[2]
    time_and_mem <- GetResourceUsage(jobId = job_id) # Returns vector of length 2
    time <- time_and_mem[1]
    mem <- time_and_mem[2]
    num_var <- str_squish(num_var)
    num_var <- as.numeric(gsub("variant count: ", "", num_var))
    resource_df[row_num + 6, ] <-  
      c(time, mem, project_name, "Ancestral Reconstruction", "Tree Provided", 10, num_var)
  }
  row_num <- row_num + 7
}
write.csv(resource_df, 
          file = "../data/prewas_resource_usage.csv", 
          quote = FALSE, 
          row.names = FALSE)
