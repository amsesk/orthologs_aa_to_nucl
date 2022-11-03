library(tidyverse)
library(argparse)

parser = argparse::ArgumentParser()
parser$add_argument('-r', '--prots_to_remove', metavar = "prots_to_remove", action = 'store', help = 'Protein headers to remove from orthotable, because main determined they match multiple orthogroups, before determining LTPs missing orthologs.')
parser$add_argument('-o', '--orthotable', metavar = "orthotable", action = 'store', help = 'Path to the original orthotable.')
parser$add_argument('-p', '--output_path', metavar = "orthotable", action = 'store', help = 'Path to write tsv with LTPs missing orthologs per orthogroup.')
args = parser$parse_args()

duplicated_prots = read_delim(args$prots_to_remove, delim="\n", col_names = F) %>%
  separate(X1, into = c('LTP', 'pid'), sep = '[|]') %>%
  mutate(pid = paste(LTP, pid, sep = '|'))

orthotable = read_delim(args$orthotable, delim = "\t", col_names = T)

# Replace protein names with NA if orthologs_aa_to_nucl determined the protein corresponded
# multiple orthogroups
for (i in 1:length(duplicated_prots$LTP)) {
  col = orthotable %>%
      select(!!sym(duplicated_prots[i,] %>% pull(LTP)))
  colidx = which(colnames(orthotable) == duplicated_prots[i,] %>% pull(LTP))
  #print(orthotable[which(col == duplicated_prots[i,] %>% pull(pid)), colidx])
  orthotable[which(col == duplicated_prots[i,] %>% pull(pid)), colidx] = NA
  #print(orthotable[which(col == duplicated_prots[i,] %>% pull(pid)), colidx])
  #print("------")
}

out = c()
for (i in 1:length(orthotable$marker)) {
  ltp_columns = colnames(orthotable[-which(colnames(orthotable) == "marker")])
  ltp_with_missing = ltp_columns[which(is.na(orthotable[i,-1]))]
  ltp_present = ltp_columns[which(!is.na(orthotable[i,-1]))]
  out = rbind(out, c(orthotable$marker[i], paste(ltp_with_missing,collapse=","), paste(ltp_present,collapse=",")))
}
out = as_tibble(out) %>%
  rename(marker = V1, ltp_missing_ortholog = V2, ltp_ortholog_present = V3)

write_delim(x = out, file = args$output_path, delim = '\t', quote = 'none', escape = 'none', col_names = T)
