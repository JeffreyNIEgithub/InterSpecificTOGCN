#!Rscript

.libPaths(c("/home/.opt/sysoft/R-4.0.2/lib64/R/library", "/home/shwzhao/R/x86_64-pc-linux-gnu-library/4.0"))

options(tidyverse.quiet=TRUE)
library(tidyverse)

args <- commandArgs(T)

cat("\n# TPM matrix do not need header.\n")
cat("# Keep sp1 sp2 order right.\n\n")
if(length(args) != 9){
cat("ERROR: Nine location parameters needed.\n")
cat("Rscript TOmatrix.R \\\n\tsp1_sp2.geneid.tsv \\\n\tsp1.TPM.tsv\tsp2.TPM.tsv \\\n\tsp1.DE.TPM.tsv\tsp2.DE.TPM.tsv \\\n\tsp1.TF.tsv\tsp2.TF.tsv \\\n\tsp1.group.tsv\tsp2.group.tsv\n\n")
q()
}

cat("Rscript TOmatrix.R \\\n\t",args[1],"\\\n\t",args[2], args[3], "\\\n\t", args[4], args[5],"\\\n\t", args[6], args[7], "\\\n\t", args[8], args[9], "\n\n")

##################################
# data read

cat("# Data reading...\n\n")

sp1_sp2_geneid <- read_tsv(args[1], col_names = c("sp1", "sp2"), col_types = cols())

sp1_info <- read_tsv(args[8], col_types = cols())
sp2_info <- read_tsv(args[9], col_types = cols())
sp1_info_num <- length(row.names(sp1_info))
sp2_info_num <- length(row.names(sp2_info))

sp1_TPM <- read_tsv(args[2], col_names = c("sp1", sp1_info$Sample), col_types = cols())
sp2_TPM <- read_tsv(args[3], col_names = c("sp2", sp2_info$Sample), col_types = cols())

sp1_TPM_DE <- read_tsv(args[4], col_names = c("sp1", sp1_info$Sample), col_types = cols())
sp2_TPM_DE <- read_tsv(args[5], col_names = c("sp2", sp2_info$Sample), col_types = cols())

sp1_TF <- read_tsv(args[6], col_names = c("sp1", "TF"), col_types = cols())
sp2_TF <- read_tsv(args[7], col_names = c("sp2", "TF"), col_types = cols())

##################################
# 创建文件
time <- str_split(Sys.time(), " ")[[1]]

output <- str_replace_all(str_c("result_", time[1], "_", time[2], sep = ""), ":", "_")
output_sp1 <- str_c(output, "/sp1/", sep = "")
output_sp2 <- str_c(output, "/sp2/", sep = "")
output_sp12 <- str_c(output, "/sp1_sp2/", sep = "")

dir.create(output)
dir.create(output_sp1)
dir.create(output_sp2)
dir.create(output_sp12)

workingdir <- getwd()

##################################
# 各种matix

cat("# Matrix generating...\n\n")

#####
# 根据配对基因id构建matrix
sp1_sp2_TPM <- sp1_sp2_geneid %>%
  inner_join(sp1_TPM, by = "sp1") %>%
  inner_join(sp2_TPM, by = "sp2")

# write_tsv(sp1_sp2_TPM, "sp1_sp2.TPM.tsv")

#####
# 提取双物种TOmatrix
# 更换基因id，软件的id不能过长
sp1_sp2_TPM_TO <- sp1_sp2_TPM %>%
  filter(sp1 %in% sp1_TPM_DE$sp1 | sp2 %in% sp2_TPM_DE$sp2) %>%
  mutate(sp12 = 1, sp12 = str_c("Gene", str_pad(cumsum(sp12), 7, "left", pad = "0")))

# id的对应关系写出
write_tsv(select(sp1_sp2_TPM_TO, sp12, sp1, sp2), str_c(output_sp12, "sp1_sp2.newid.tsv", sep = ""))
write_tsv(select(sp1_sp2_TPM_TO, sp12, everything(), -sp1, -sp2), str_c(output_sp12, "sp1_sp2.TPM.TO.tsv", sep = ""), col_names = F)

#####
# 提取双物种TF TOmatrix
sp1_sp2_TPM_TO_TF <- sp1_sp2_TPM_TO %>%
  filter(sp1 %in% sp1_TF$sp1 | sp2 %in% sp2_TF$sp2) %>%
  select(sp12, everything(), -sp1, -sp2)

write_tsv(sp1_sp2_TPM_TO_TF, str_c(output_sp12, "sp1_sp2.TPM.TO.TF.tsv", sep = ""), col_names = F)

#####
# 提取单物种TOmtrix
# 注意：单物种除了配对id，还有独有id
sp1_TPM_TO <- sp1_TPM %>%
  filter(sp1 %in% sp1_sp2_TPM_TO$sp1 | sp1 %in% sp1_TPM_DE$sp1) %>%
  left_join(sp1_TPM, by = "sp1")

sp2_TPM_TO <- sp2_TPM %>%
  filter(sp2 %in% sp1_sp2_TPM_TO$sp2 | sp2 %in% sp2_TPM_DE$sp2) %>%
  left_join(sp2_TPM, by = "sp2")

write_tsv(sp1_TPM_TO, str_c(output_sp1, "sp1.TPM.TO.tsv", sep = ""), col_names = F)
write_tsv(sp2_TPM_TO, str_c(output_sp2, "sp2.TPM.TO.tsv", sep = ""), col_names = F)

#####
# 提取单物种TF TOmatrix
sp1_TPM_TO_TF <- sp1_TPM_TO %>%
  filter(sp1 %in% sp1_TF$sp1)

sp2_TPM_TO_TF <- sp2_TPM_TO %>%
  filter(sp2 %in% sp2_TF$sp2)

write_tsv(sp1_TPM_TO_TF, str_c(output_sp1, "sp1.TPM.TO.TF.tsv", sep = ""), col_names = F)
write_tsv(sp2_TPM_TO_TF, str_c(output_sp2, "sp2.TPM.TO.TF.tsv", sep = ""), col_names = F)

# Usage: Cutoff #Cond1_samples #Cond2_samples file_of_TF_genes file_of_all_genes

cat("# Command you may need run next:\n")
cat("cd", str_c(workingdir, "/", output_sp1, sep = ""), "\n")
cat("Cutoff", sp1_info_num, sp1_info_num, "sp1.TPM.TO.TF.tsv", "sp1.TPM.TO.tsv", "\n")
cat("GCN", sp1_info_num, sp1_info_num, "sp1.TPM.TO.TF.tsv", "sp1.TPM.TO.tsv", "$Cutoff_pos1", "$Cutoff_pos1", "$Cutoff_neg", "$Cutoff_neg", "\n")
cat("echo $sp1_seed > seed.txt\n")
cat("TO-GCN", sp1_info_num, sp1_info_num, "sp1.TPM.TO.tsv", "sp1.TPM.TO.tsv", "$Cutoff_pos1", "$Cutoff_pos1", "seed.txt 0", "\n\n")

cat("cd", str_c(workingdir, "/", output_sp2, sep = ""), "\n")
cat("Cutoff", sp2_info_num, sp2_info_num, "sp2.TPM.TO.TF.tsv", "sp2.TPM.TO.tsv", "\n")
cat("GCN", sp2_info_num, sp2_info_num, "sp2.TPM.TO.TF.tsv", "sp2.TPM.TO.tsv", "$Cutoff_pos2", "$Cutoff_pos2", "$Cutoff_neg", "$Cutoff_neg", "\n")
cat("echo $sp2_seed > seed.txt\n")
cat("TO-GCN", sp2_info_num, sp2_info_num, "sp2.TPM.TO.tsv", "sp2.TPM.TO.tsv", "$Cutoff_pos2", "$Cutoff_pos2", "seed.txt 0", "\n\n")

cat("cd", str_c(workingdir, "/", output_sp12, sep = ""), "\n")
cat("Cutoff", sp1_info_num, sp2_info_num, "sp1_sp2.TPM.TO.TF.tsv", "sp1_sp2.TPM.TO.tsv", "\n")
cat("Rscript $TOscript/TOcheckseed.R", sp1_info_num, sp2_info_num, "$sp1_seed $sp2_seed sp1_sp2.newid.tsv sp1_sp2.TPM.TO.TF.tsv sp1_sp2.TPM.TO.tsv\n")
cat("GCN", sp1_info_num, sp2_info_num, "sp1_sp2.TPM.TO.TF.2.tsv", "sp1_sp2.TPM.TO.2.tsv", "$Cutoff_pos1", "$Cutoff_pos2", "$Cutoff_neg", "$Cutoff_neg", "\n")
cat("TO-GCN", sp1_info_num, sp2_info_num, "sp1_sp2.TPM.TO.2.tsv", "sp1_sp2.TPM.TO.2.tsv", "$Cutoff_pos1", "$Cutoff_pos2", "seed.txt 0", "\n\n")

cat("# Network you need to combine after all 'GCN' steps finished:\n")
cat("cd", str_c(workingdir, "/", output, sep = ""), "\n")
cat("Rscript $TOscript/TOcombinenet.R sp1_sp2/sp1_sp2.newid.2.tsv sp1/C1+C2+.csv sp2/C1+C2+.csv sp1_sp2/C1+C2+.csv\n\n")
