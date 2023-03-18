#!Rscript

.libPaths(c("/home/.opt/sysoft/R-4.0.2/lib64/R/library", "/home/shwzhao/R/x86_64-pc-linux-gnu-library/4.0"))
cat("# Network combining script.\n\n")

options(tidyverse.quiet=TRUE)
library(tidyverse)

args <- commandArgs(T)

if(length(args) != 4){
cat("ERROR: Four parameters needed.\n\n")
cat("Rscript TOcombinenet.R \\\n\tsp1_sp2/sp1_sp2.newid.2.tsv \\\n\tsp1/C1+C2+.csv \\\n\tsp2/C1+C2+.csv \\\n\tsp1_sp2/C1+C2+.csv\n\n")
q()
}

cat("Rscript TOcombinenet.R \\\n\t", args[1], " \\\n\t", args[2],"\\\n\t", args[3], "\\\n\t", args[4], "\n\n")

################################
# 读取数据
cat("\n# Data reading...\n\n")

sp1_sp2_newid_x <- read_tsv(args[1], col_types = cols())

sp1_net <- read_csv(args[2], col_types = cols()) %>%
  select(source_node = `TF gene ID`, target_node = `gene ID`) %>%
  mutate(type = "sp1")

sp2_net <- read_csv(args[3], col_types = cols()) %>%
  select(source_node = `TF gene ID`, target_node = `gene ID`) %>%
  mutate(type = "sp2")

sp12_net <- read_csv(args[4], col_types = cols()) %>%
  select(source_node = `TF gene ID`, target_node = `gene ID`) %>%
  mutate(type = "sp1_sp2")

################################
# 提取网络之间的连线
between_sp1 <- sp1_sp2_newid_x %>%
  filter(sp12 %in% unique(c(sp12_net$source_node, sp12_net$target_node)) & sp1 %in% unique(c(sp1_net$source_node, sp1_net$target_node))) %>%
  select(source_node = sp12, target_node = sp1) %>%
  mutate(type = "between_sp1")

between_sp2 <- sp1_sp2_newid_x %>%
  filter(sp12 %in% unique(c(sp12_net$source_node, sp12_net$target_node)) & sp2 %in% unique(c(sp2_net$source_node, sp2_net$target_node))) %>%
  select(source_node = sp12, target_node = sp2) %>%
  mutate(type = "between_sp2")

################################
# 合并
all_net <- bind_rows(sp1_net, sp2_net, sp12_net, between_sp1, between_sp2)

write_tsv(all_net, "network.tsv")

cat("# Done.\n\n")
