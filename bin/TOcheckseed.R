#!Rscript

.libPaths(c("/home/.opt/sysoft/R-4.0.2/lib64/R/library", "/home/shwzhao/R/x86_64-pc-linux-gnu-library/4.0"))

options(tidyverse.quiet=TRUE)
library(tidyverse)

args <- commandArgs(T)

if(length(args) != 7){
cat("\nERROR: Seven parameters needed.\n\n")
cat("Rscript TOcheckseed.R \\\n\t$sp1_num $sp2_num \\\n\t$sp1_seed $sp2_seed \\\n\tsp1_sp2.newid.tsv \\\n\tsp1_sp2.TPM.TO.TF.tsv \\\n\tsp1_sp2.TPM.TO.tsv\n\n")
q()
}

cat("Rscript TOcheckseed.R \\\n\t", args[1], args[2], "\\\n\t", args[3], args[4],"\\\n\t", args[5], "\\\n\t", args[6], "\\\n\t", args[7], "\n")

###################################
# 读取数据
sp1_info <- as.numeric(args[1])
sp2_info <- as.numeric(args[2])

sp1_seed <- args[3]
sp2_seed <- args[4]

sp1_sp2_newid <- read_tsv(args[5], col_types = cols())
sp1_sp2_TPM_TO_TF <- read_tsv(args[6], col_names = F, col_types = cols())
sp1_sp2_TPM_TO <- read_tsv(args[7], col_names = F, col_types = cols())

###################################
# 判断是否有合适的seed
seed <- sp1_sp2_newid %>%
  filter(sp1 == sp1_seed & sp2 == sp2_seed) %>%
  .$sp12 %>%
  unique()

if(length(seed) == 0){
# 如果没有

# 提取sp1的表达量
a <- sp1_info + 1
b <- sp1_info + 2
c <- sp1_info + sp2_info + 1

sp1_seed_x <- sp1_sp2_TPM_TO %>%
  filter(X1 == unique(filter(sp1_sp2_newid, sp1 == sp1_seed)$sp12)) %>%
  mutate(X1 = "Gene0000000") %>%
  select(1:all_of(a))

# 提取sp2的表达量
sp2_seed_x <- sp1_sp2_TPM_TO %>%
  filter(X1 == unique(filter(sp1_sp2_newid, sp2 == sp2_seed)$sp12)) %>%
  mutate(X1 = "Gene0000000") %>%
  select(1, all_of(b):all_of(c))

# 创建的seed
seed_x <- sp1_seed_x %>%
  left_join(sp2_seed_x, by = "X1")

# 写入id对应表
sp1_sp2_newid_x <- sp1_sp2_newid %>%
  add_row(sp12 = "Gene0000000", sp1 = sp1_seed, sp2 = sp1_seed)

write_tsv(sp1_sp2_newid_x, "sp1_sp2.newid.2.tsv")

# 写入新的TPM表达矩阵
write_tsv(bind_rows(sp1_sp2_TPM_TO_TF, seed_x), "sp1_sp2.TPM.TO.TF.2.tsv", col_names = F)
write_tsv(bind_rows(sp1_sp2_TPM_TO, seed_x), "sp1_sp2.TPM.TO.2.tsv", col_names = F)

# 写入seed.txt
write("Gene0000000", "seed.txt")
} else{
# 如果有理想的seed

# 直接复制到新文件
write_tsv(sp1_sp2_newid, "sp1_sp2.newid.2.tsv")
write_tsv(sp1_sp2_TPM_TO_TF, "sp1_sp2.TPM.TO.TF.2.tsv", col_names = F)
write_tsv(sp1_sp2_TPM_TO, "sp1_sp2.TPM.TO.2.tsv", col_names = F)
write(seed, "seed.txt")
}

cat("\n# Done.\n")
cat("\n# Check files: sp1_sp2.newid.2.tsv, sp1_sp2.TPM.TO.TF.2.tsv, sp1_sp2.TPM.TO.2.tsv, seed.txt\n\n")
