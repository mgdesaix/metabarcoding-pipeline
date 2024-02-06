library(tidyverse)

blast <- read_table("./03_results/04_summary/blast_summary.txt",
                    col_names = c("Sample", "OTU", "Depth", "something",
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
write_csv(x = blast_sp,
          file = "./03_results/04_summary/blast_species_summary.txt")

blast_sample_depth <- blast %>%
  group_by(Sample) %>%
  summarize(Depth = sum(Depth),
            .groups = "drop")
write_csv(x = blast_sample_depth,
          file = "./03_results/04_summary/blast_sample_depth.txt")
