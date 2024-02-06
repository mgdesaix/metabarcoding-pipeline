library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

if (length(args)!=1) {
  stop("Need to provide one parameter - blast_summary.txt or equivalent! (input file).n", call.=FALSE)
}

blast <- read_table(args[1],
                    col_names = c("Sample", "OTU", "Depth", "GenBank_ID",
                                  "Kingdom", "Phylum", "Class", "Order", "Family",
                                  "Genus", "Species", "Match"),
                    show_col_types = F)


blast_sp <- blast %>%
  group_by(Sample, Kingdom, Phylum, Class, Order,
           Family, Genus, Species) %>%
  summarize(Depth_sum = sum(Depth),
            # N_OTUs = n(),
            .groups = "drop") %>%
  # filter(Depth > 5) %>%
  arrange(Sample, desc(Depth_sum)) %>%
  pivot_wider(names_from = Sample, values_from = Depth_sum)
write_delim(x = blast_sp,
          file = "blast_species_contingency_table.txt",
           delim="\t")

blast_sample_depth <- blast %>%
  group_by(Sample) %>%
  summarize(Depth = sum(Depth),
            .groups = "drop")
write_delim(x = blast_sample_depth,
          file = "blast_sample_depth.txt",
           delim="\t")


otus <- read_table("./full-sample-otus.txt",
                   col_names = c("Data", "Sample")) %>%
  mutate("Name" = rep(c("OTU", "Sequence"), length.out = n()))

otus_wide <- cbind((otus %>%
                      filter(Name == "OTU") %>%
                      select(Sample, Data) %>%
                      rename("OTU" = "Data")),
                   (otus %>%
                      filter(Name == "Sequence") %>%
                      select(Data) %>%
                      rename("Sequence" = "Data"))) %>%
  tibble()

blast_otus <- blast %>%
  left_join(otus_wide, by = c("Sample", "OTU"))

write_delim(x = blast_otus,
            file = "./blast_OTU_summary.txt",
            delim = "\t")

