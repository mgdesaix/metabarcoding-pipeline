library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

if (length(args)!=1) {
  stop("Need to provide one parameter - blast_summary.txt or equivalent! (input file).n", call.=FALSE)
}

blast <- read_table(args[1],
                    col_names = c("Sample", "OTU", "Depth", "GenBank_ID",
                                  "Kingdom", "Phylum", "Class", "Order", "Family",
                                  "Genus", "Species", "Match", "Query", "Target"),
                    show_col_types = F)

blast_sp <- blast %>%
  group_by(Sample, OTU, Species, Match) %>%
  slice_sample(n = 1)

write_delim(x = blast_sp,
          file = "blast_verbose_species_hits.txt",
           delim="\t")
