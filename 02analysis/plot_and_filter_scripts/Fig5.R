library(tidyverse)
library(janitor)
library(showtext)
library(tidygraph)
library(igraph)
library(UpSetR)
library(dplyr)
library(ggplot2)
library(forcats)
library(svglite)


# R Workflow Author: A. Fernandez-Guerra (Lundbeck Foundation GeoGenetics Centre)
# Modifications By: Louis Kraft (DTU Denmark - Section of Bioinformatics)

showtext_auto()


# Get contigs stats
contigs_files <- list.files("/home/carpedeamManuscript", full.names = TRUE, pattern = ".contig_stats.tsv")

read_contig_files <- function(filename) {
  pattern <- "^(.*)_(.*)\\.contig_stats\\.tsv$"
  matches <- str_match(basename(filename), pattern)
  
  # Extract the matched groups
  sample <- matches[2]
  assembler <- matches[3]
  df <- read_tsv(filename, col_names = FALSE) |>
    setNames(c("contig_id", "length")) |>
    mutate(sample = sample, assembler = assembler)
  if (nrow(df) == 0) {
    return(NULL)
  }
  return(df)
}

contigs <- map_dfr(contigs_files, read_contig_files)

min_contig_length <- 1000

contig_stats <- contigs |>
  filter(length >= min_contig_length) |>
  group_by(sample, assembler) |>
  summarise(
    n_contigs = n(),
    total_length = sum(length),
    mean_length = mean(length),
    median_length = median(length),
    min_length = min(length),
    max_length = max(length),
    .groups = "drop"
  )


# Which is the level of redundancy in the different assemblies based on the ANI?
skani_contigs_files <- list.files("/home/carpedeamManuscript", full.names = TRUE, pattern = "\\.ani.tsv")

read_skani_contig_file <- function(filename) {
  pattern <- "^(.*)_(.*)\\.ani\\.tsv$"
  matches <- str_match(basename(filename), pattern)
  
  # Extract the matched groups
  sample <- matches[2]
  assembler <- matches[3]
  df <- read_tsv(filename) |>
    clean_names() |>
    mutate(sample = sample, assembler = assembler)
  if (nrow(df) == 0) {
    return(NULL)
  }
  return(df)
}

skani_contig <- map_dfr(skani_contigs_files, read_skani_contig_file)

skani_contig_edges <- skani_contig |>
  inner_join(contigs |> rename(ref_name = contig_id, ref_length = length)) |>
  inner_join(contigs |> rename(query_name = contig_id, query_length = length)) |>
  filter(query_length >= min_contig_length, ref_length >= min_contig_length) |>
  filter(ani >= 99) |>
  filter(if_else(ref_length <= query_length, align_fraction_ref, align_fraction_query) >= 90) |>
  select(sample, assembler, ref_name, query_name, ani, align_fraction_ref, align_fraction_query, ref_length, query_length)

s_a <- skani_contig_edges |>
  select(sample, assembler) |>
  distinct()


# Find non-overlapping cliques
find_non_overlapping_cliques <- function(cliques) {
  non_overlapping_cliques <- list()
  used_nodes <- c()
  
  for (clique in cliques) {
    if (length(intersect(clique, used_nodes)) == 0) {
      non_overlapping_cliques <- append(non_overlapping_cliques, list(clique))
      used_nodes <- union(used_nodes, clique)
    }
  }
  
  return(non_overlapping_cliques)
}

# Create a graph for each sample and assembler
g <- s_a |>
  mutate(graph = map2(sample, assembler, ~ {
    # Create the graph
    graph <- skani_contig_edges |>
      filter(sample == .x, assembler == .y) |>
      select(ref_name, query_name, weight = ani, align_fraction_ref, align_fraction_query) |>
      as_tbl_graph(directed = FALSE) |>
      mutate(component = group_components()) |>
      inner_join(contigs |> filter(sample == .x, assembler == .y) |> select(name = contig_id, length))
    
    # Detect cliques and find non-overlapping cliques
    cliques <- max_cliques(graph)
    non_overlapping_cliques <- find_non_overlapping_cliques(cliques)
    
    # Assign clique membership to nodes
    membership <- rep(NA, vcount(graph))
    for (i in seq_along(non_overlapping_cliques)) {
      for (node in non_overlapping_cliques[[i]]) {
        membership[node] <- i
      }
    }
    
    # Add clique membership to the graph
    graph <- graph |>
      activate(nodes) |>
      mutate(clique_membership = membership)
    
    return(graph)
  }))


# Function to get components and their sizes
get_component_info <- function(graph) {
  components <- graph |>
    as_tibble() |>
    group_by(component) |>
    summarise(size = n(), nodes = list(name)) |>
    arrange(desc(size))
  
  num_components <- n_distinct(components$component)
  
  list(components = components, num_components = num_components, num_nodes = igraph::vcount(graph))
}

# Apply the function to each graph in the list and get stats
graph_stats <- g |>
  mutate(
    component_info = map(graph, get_component_info),
    component_size = map(component_info, "components"),
    num_components = map_int(component_info, "num_components"),
    num_nodes = map_int(component_info, "num_nodes"),
  ) |>
  inner_join(contig_stats |> select(sample, assembler, n_contigs), by = c("sample", "assembler"))

# Where do this contigs map in the reference genomes?
skani_search_files <- list.files("/home/carpedeamManuscript", full.names = TRUE, pattern = "\\.search-skani.tsv")
read_skani_search_file <- function(filename) {
  pattern <- "^(.*)_(.*)\\.search-skani\\.tsv$"
  matches <- str_match(basename(filename), pattern)
  
  # Extract the matched groups
  sample <- matches[2]
  assembler <- matches[3]
  df <- read_tsv(filename) |>
    clean_names() |>
    mutate(sample = sample, assembler = assembler)
  if (nrow(df) == 0) {
    return(NULL)
  }
  return(df)
}

# Get the data of mapping the contigs to the original assemblies
skani_search <- map_dfr(skani_search_files, read_skani_search_file) |>
  filter(align_fraction_query >= 95, ani >= 99) |>
  inner_join(contigs |> rename(query_name = contig_id, qlen = length), by = c("sample", "assembler", "query_name")) |>
  filter(qlen > min_contig_length)


# Combine search results with grpah nodes
join_nodes_with_search <- function(graph, sample, assembler) {
  hits <- skani_search |>
    filter(sample == sample, assembler == assembler) |>
    group_by(query_name) |>
    filter(
      ani == max(ani) &
        align_fraction_query == max(align_fraction_query)
    ) |>
    slice(1) |>
    ungroup() |>
    select(query_name, ref_name) |>
    distinct()
  
  
  graph |>
    activate(nodes) |>
    as_tibble() |>
    rename(query_name = name) |>
    left_join(hits) |>
    mutate(is_mapped = ifelse(!is.na(ref_name), "mapped", "unmapped")) |>
    # filter(component == 1, clique_membership == 28) |>
    arrange(clique_membership, desc(length)) |>
    group_by(ref_name, component, clique_membership) |>
    summarise(
      contig_num = n(),
      max_contig_bp = max(length),
      sum_contig_bp = sum(length),
      unmapped_bp = sum(ifelse(is.na(ref_name), length, 0)),
      .groups = "drop"
    ) |>
    mutate(is_mapped = ifelse(!is.na(ref_name), "mapped", "unmapped")) |>
    group_by(clique_membership) |>
    add_count(name = "ref_num") |>
    ungroup() |>
    mutate(type = ifelse(ref_num > 1, "multi", "single"))
}

# Apply the function to each row in graph_stats
graph_stats <- graph_stats |>
  mutate(joined_nodes = pmap(list(graph, sample, assembler), join_nodes_with_search))


# Calculate for each assembly|sample the number of contigs that map to more than one reference genome
# calculate the number of contigs that are duplicated and need to be removed
get_ncontigs_multiple_refs <- function(joined_nodes) {
  df <- joined_nodes |>
    mutate(type = ifelse(is_mapped == "unmapped", "unmapped", type)) |>
    group_by(type) |>
    summarise(
      n = n(),
      contig_num = sum(contig_num),
      contig_unmmaped_bp = sum(unmapped_bp),
      contig_mapped_bp = sum(sum_contig_bp) - contig_unmmaped_bp,
      contig_bp_total = sum(sum_contig_bp),
      contig_bp_keep = sum(max_contig_bp),
      contig_bp_rm = contig_bp_total - contig_bp_keep,
      .groups = "drop"
    ) |>
    mutate(contig2rm = ifelse(type == "multi", 0, contig_num - n))
}

graph_stats <- graph_stats |>
  mutate(
    node_summary = map(joined_nodes, get_ncontigs_multiple_refs),
    n_contigs_multi = map_int(node_summary, ~ {
      # Filter rows where type is "multi"
      multi_rows <- .x |> filter(type == "multi")
      # Extract the first contig_num of the filtered rows (if any)
      if (nrow(multi_rows) > 0) multi_rows$contig_num[1] else 0
    }),
    n_contigs_single = map_int(node_summary, ~ {
      # Filter rows where type is "single"
      single_rows <- .x |> filter(type == "single")
      # Extract the first contig_num of the filtered rows (if any)
      if (nrow(single_rows) > 0) single_rows$contig_num[1] else 0
    }),
    n_contigs_unmapped = map_int(node_summary, ~ {
      # Filter rows where type is "single"
      single_rows <- .x |> filter(type == "unmapped")
      # Extract the first contig_num of the filtered rows (if any)
      if (nrow(single_rows) > 0) single_rows$contig_num[1] else 0
    }),
    n_contigs_single_rm = map_int(node_summary, ~ {
      # Filter rows where type is "single"
      single_rows <- .x |> filter(type == "single")
      # Extract the first contig_num of the filtered rows (if any)
      if (nrow(single_rows) > 0) single_rows$contig2rm[1] else 0
    }),
    n_contigs_multi_rm = map_int(node_summary, ~ {
      # Filter rows where type is "multi"
      single_rows <- .x |> filter(type == "multi")
      # Extract the first contig_num of the filtered rows (if any)
      if (nrow(single_rows) > 0) single_rows$contig2rm[1] else 0
    }),
    n_contigs_unmapped_rm = map_int(node_summary, ~ {
      # Filter rows where type is "single"
      single_rows <- .x |> filter(type == "unmapped")
      # Extract the first contig_num of the filtered rows (if any)
      if (nrow(single_rows) > 0) single_rows$contig2rm[1] else 0
    }),
    contig_bp_single_keep = map_int(node_summary, ~ {
      # Filter rows where type is "single"
      single_rows <- .x |> filter(type == "single")
      # Extract the first contig_num of the filtered rows (if any)
      if (nrow(single_rows) > 0) single_rows$contig_bp_keep[1] else 0
    }),
    contig_bp_single_rm = map_int(node_summary, ~ {
      # Filter rows where type is "single"
      single_rows <- .x |> filter(type == "single")
      # Extract the first contig_num of the filtered rows (if any)
      if (nrow(single_rows) > 0) single_rows$contig_bp_rm[1] else 0
    }),
    contig_bp_multi_keep = map_int(node_summary, ~ {
      # Filter rows where type is "multi"
      single_rows <- .x |> filter(type == "multi")
      # Extract the first contig_num of the filtered rows (if any)
      if (nrow(single_rows) > 0) single_rows$contig_bp_keep[1] else 0
    }),
    contig_bp_multi_rm = map_int(node_summary, ~ {
      # Filter rows where type is "multi"
      single_rows <- .x |> filter(type == "multi")
      # Extract the first contig_num of the filtered rows (if any)
      if (nrow(single_rows) > 0) single_rows$contig_bp_rm[1] else 0
    }),
    contig_bp_unmapped_keep = map_int(node_summary, ~ {
      # Filter rows where type is "single"
      single_rows <- .x |> filter(type == "unmapped")
      # Extract the first contig_num of the filtered rows (if any)
      if (nrow(single_rows) > 0) single_rows$contig_bp_keep[1] else 0
    }),
    contig_bp_unmapped_rm = map_int(node_summary, ~ {
      # Filter rows where type is "single"
      single_rows <- .x |> filter(type == "unmapped")
      # Extract the first contig_num of the filtered rows (if any)
      if (nrow(single_rows) > 0) single_rows$contig_bp_rm[1] else 0
    })
  )

# Get all stats
mapping_stats <- skani_search |>
  select(sample, assembler, query_name, qlen) |>
  distinct() |>
  group_by(sample, assembler) |>
  summarise(mapped_contigs = n(), mapped_bp = sum(qlen), .groups = "drop") |>
  inner_join(contig_stats |> select(sample, assembler, n_contigs, total_length), by = c("sample", "assembler")) |>
  left_join(graph_stats |>
              select(
                sample,
                assembler,
                n_contigs_multi,
                n_contigs_single,
                n_contigs_unmapped,
                n_contigs_single_rm,
                n_contigs_multi_rm,
                n_contigs_unmapped_rm,
                contig_bp_single_keep,
                contig_bp_single_rm,
                contig_bp_multi_keep,
                contig_bp_multi_rm,
                contig_bp_unmapped_keep,
                contig_bp_unmapped_rm
              ), by = c("sample", "assembler")) |>
  replace_na(list(
    n_contigs = 0,
    n_contigs_multi = 0,
    n_contigs_single = 0,
    n_contigs_unmapped = 0,
    n_contigs_single_rm = 0,
    n_contigs_multi_rm = 0,
    n_contigs_unmapped_rm = 0,
    contig_bp_single_keep = 0,
    contig_bp_single_rm = 0,
    contig_bp_multi_keep = 0,
    contig_bp_multi_rm = 0,
    contig_bp_unmapped_keep = 0,
    contig_bp_unmapped_rm = 0
  ))


# Do figure 2
colors <- c(
  "mapped_nondup_bp" = "#2066a8",
  "mapped_dup_var_bp" = "#3594cC",
  "mapped_dup_red_rep_bp" = "#8cc5e3",
  "mapped_dup_red_bp" = "#F09868FF",
  "unmapped_dup_bp" = "#c46666",
  "unmapped_nondup_bp" = "#d8a6a6"
)


# Adjust the base font size for better readability on an A4 page
base_font_size <- 14

mapping_stats |>
  mutate(
    # Total unmapped bp
    unmapped_bp = total_length - mapped_bp,
    # Total unmapped bp that are not duplicates
    unmapped_nondup_bp = unmapped_bp - (contig_bp_unmapped_keep + contig_bp_unmapped_rm),
    # Total unmapped bp that are duplicates
    unmapped_dup_bp = contig_bp_unmapped_keep + contig_bp_unmapped_rm,
    # Total mapped bp that are duplicates and can be due to variation
    mapped_dup_var_bp = contig_bp_multi_keep,
    # Total mapped bp that are duplicates and can be due to redundancy
    mapped_dup_red_bp = contig_bp_single_rm + contig_bp_multi_rm,
    # Total mapped bp that from the duplicates that are representative
    mapped_dup_red_rep_bp = contig_bp_single_keep,
    # Total mapped bp that are duplicates
    mapped_dup_bp = mapped_dup_var_bp + mapped_dup_red_bp + mapped_dup_red_rep_bp,
    # Total mapped bp that are not duplicates
    mapped_nondup_bp = mapped_bp - mapped_dup_bp
  ) |>
  select(
    sample,
    assembler,
    mapped_nondup_bp,
    mapped_dup_var_bp,
    mapped_dup_red_bp,
    mapped_dup_red_rep_bp,
    unmapped_nondup_bp,
    unmapped_dup_bp
  ) |>
  pivot_longer(cols = c(
    mapped_dup_var_bp,
    mapped_dup_red_rep_bp,
    mapped_dup_red_bp,
    mapped_nondup_bp,
    unmapped_dup_bp,
    unmapped_nondup_bp
  ), names_to = "name", values_to = "bp") |>
  mutate(type = ifelse(grepl("mapped", name), "mapped", "unmapped")) |>
  mutate(
    cov = case_when(
      grepl("c3", sample) ~ "3X",
      grepl("c5", sample) ~ "5X",
      grepl("c10", sample) ~ "10X",
      TRUE ~ "unknown"
    ),
    dmg = case_when(
      grepl("high", sample) ~ "high",
      TRUE ~ "unknown"
    ),
    assm = case_when(
      grepl("carpedeam-unsafe", assembler) ~ "Carpedeam\n(unsafe)",
      grepl("carpedeam-safe", assembler) ~ "Carpedeam\n(safe)",
      grepl("megahit", assembler) ~ "Megahit",
      grepl("spades", assembler) ~ "metaSPAdes",
      grepl("penguin", assembler) ~ "PenguiN",
      TRUE ~ "unknown"
    )
  ) |>
  mutate(
    name = fct_relevel(name, c("unmapped_dup_bp", "unmapped_nondup_bp", "mapped_dup_red_bp", "mapped_dup_red_rep_bp", "mapped_dup_var_bp", "mapped_nondup_bp")),
    cov = fct_relevel(cov, c("3X", "5X", "10X")),
    dmg = fct_relevel(dmg, c("high"))
  ) |>
  ggplot(aes(x = assm, y = bp, fill = name)) +
  geom_bar(stat = "identity", position = "stack", width = 1, color = "black", linewidth = 0.3) +
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6)) +
  scale_fill_manual(values = colors, name = NULL) +
  facet_grid(dmg ~ cov, labeller = "label_both", scales = "free_y") +
  theme_bw(base_size = base_font_size) +  # Adjust base_size to change the font size
  theme(
    strip.background = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  ylab("Base pairs") +
  xlab("")

# Save the plot as an SVG file
ggsave("Figure2.png", width = 10, height = 7, dpi = 300)



# Why if there are more contigs in CP we have the same number of proteins?

gff_files <- list.files("/home/carpedeamManuscript", full.names = TRUE, pattern = "gff")

read_gff_files <- function(filename) {
  pattern <- "^(.*)_(.*)\\.gff$"
  matches <- str_match(basename(filename), pattern)
  
  # Extract the matched groups
  sample <- matches[2]
  assembler <- matches[3]
  
  df <- read_tsv(filename, col_names = FALSE, comment = "#") |>
    clean_names() |>
    select(contig_id = x1, start = x4, end = x5) |>
    mutate(sample = sample, assembler = assembler) |>
    mutate(ngene = row_number()) |>
    mutate(assm = case_when(
      assembler == "carpedeam-safe" ~ "cp_s",
      assembler == "carpedeam-unsafe" ~ "cp_u",
      assembler == "megahit" ~ "mh"
    )) |>
    mutate(prot_id = paste(sample, assm, ngene, sep = "-")) |>
    mutate(gene_length = end - start + 1) |>
    select(sample, assembler, assembler, contig_id, prot_id, gene_length)
  
  return(df)
}


# Get prodigal GFF files
gff_data <- map_dfr(gff_files, read_gff_files)

gff_data |>
  group_by(sample, assembler, contig_id) |>
  summarize(n = n(), coding_bp = sum(gene_length), .groups = "drop") |>
  inner_join(contigs, by = c("sample", "assembler", "contig_id")) |>
  mutate(coding_density = coding_bp / length) |>
  filter(coding_density < 1 / 2) |>
  group_by(sample, assembler) |>
  summarise(n = n(), length = sum(length), .groups = "drop") |>
  mutate(
    cov = case_when(
      grepl("c3", sample) ~ "3X",
      grepl("c5", sample) ~ "5X",
      grepl("c10", sample) ~ "10X",
      TRUE ~ "unknown"
    ),
    dmg = case_when(
      grepl("high", sample) ~ "high",
      TRUE ~ "unknown"
    ),
    assm = case_when(
      grepl("carpedeam-unsafe", assembler) ~ "Carpedeam\n(unsafe)",
      grepl("carpedeam-safe", assembler) ~ "Carpedeam\n(safe)",
      grepl("megahit", assembler) ~ "Megahit",
      grepl("spades", assembler) ~ "metaSPAdes",
      grepl("penguin", assembler) ~ "PenguiN",
      TRUE ~ "unknown"
    )
  ) |>
  mutate(
    cov = fct_relevel(cov, c("3X", "5X", "10X")),
    dmg = fct_relevel(dmg, c("high"))
  ) |>
  ggplot(aes(x = n, y = length, fill = assm)) +
  geom_point(shape = 21, size = 3, color = "black") +
  # stat_summary(fun = median, geom = "line", size = 2, color = "#EB554A", group=1) +
  # stat_summary(fun = median, geom = "point", size = 2, fill = "#EB554A", shape = 21) +
  # scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3)) +
  labs(
    x = "# contigs",
    y = "Aggregated bp"
  ) +
  facet_grid(dmg ~ cov, labeller = "label_both", scales = "free") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom"
  ) +
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6)) +
  scale_fill_manual(values = c("#2D2D2D", "#EB554A", "#FFC300", "#F09868", "#91AEB7"), name = NULL)





# Get annotations for the contigs based on Prokka
prokka_gff_files <- list.files("/home/carpedeamManuscript/assembly-annotation-eval/gffProkka/clean_gff", full.names = TRUE, pattern = "gff")

# Sample mappings
sample_mappings <- c(
  "4843e2efab" = "gut_sum_high_c3",
  "5625abfa64" = "gut_sum_high_c5",
  "147b4cc2c4" = "gut_sum_high_c10"
)

read_prokka_gff_files <- function(filename) {
  pattern <- "^(.*).raw-raw.proteins.(.*)\\.gff$"
  matches <- str_match(basename(filename), pattern)
  
  # Extract the matched groups
  sample <- matches[2]
  assembler <- matches[3]
  
  df <- read_tsv(filename, col_names = FALSE, comment = "#") |>
    clean_names() |>
    select(contig_id = x1, feature = x3, start = x4, end = x5) |>
    mutate(
      sample = sample_mappings[sample],  # Map the sample to the corresponding value
      assembler = case_when(
        grepl("config7020", assembler) ~ "carpedeam-safe",
        grepl("config7022", assembler) ~ "carpedeam-unsafe",
        grepl("megahit", assembler) ~ "megahit",
        grepl("penguin", assembler) ~ "penguin",
        grepl("spades", assembler) ~ "spades",
        TRUE ~ NA_character_ # Ensure all cases are covered
      ),
      assm = case_when(
        assembler == "carpedeam-safe" ~ "cp_s",
        assembler == "carpedeam-unsafe" ~ "cp_u",
        assembler == "megahit" ~ "mh",
        assembler == "spades" ~ "ms",
        assembler == "penguin" ~ "pg",
        TRUE ~ NA_character_ # Ensure all cases are covered
      )
    ) |>
    mutate(feature_length = end - start + 1) |>
    select(sample, assembler, contig_id, feature, feature_length)
  
  return(df)
}

# Read the GFFs
prokka_annotations <- map_dfr(prokka_gff_files, read_prokka_gff_files)


prokka_features <- gff_data |>
  group_by(sample, assembler, contig_id) |>
  summarize(n = n(), coding_bp = sum(gene_length), .groups = "drop") |>
  inner_join(contigs) |>
  mutate(coding_density = coding_bp / length) |>
  # filter(coding_density < 1 / 2) |>
  inner_join(skani_search |> select(contig_id = query_name, sample, assembler, ani, align_fraction_query)) |>
  inner_join(prokka_annotations |> select(contig_id, sample, assembler, feature, feature_length))

prokka_features |>
  filter(coding_density < 1 / 2) |>
  group_by(sample, assembler, feature) |>
  summarize(n = n(), length = sum(feature_length), .groups = "drop") |>
  # Replace NA length with 0
  mutate(length = ifelse(is.na(length), 0, length)) |>
  mutate(
    cov = case_when(
      grepl("c3", sample) ~ "3X",
      grepl("c5", sample) ~ "5X",
      grepl("c10", sample) ~ "10X",
      TRUE ~ NA_character_ # Ensure all cases are covered
    ),
    dmg = case_when(
      grepl("high", sample) ~ "high",
      TRUE ~ NA_character_ # Ensure all cases are covered
    ),
    assm = case_when(
      grepl("carpedeam-unsafe", assembler) ~ "Carpedeam\n(unsafe)",
      grepl("carpedeam-safe", assembler) ~ "Carpedeam\n(safe)",
      grepl("megahit", assembler) ~ "Megahit",
      grepl("spades", assembler) ~ "metaSPAdes",
      grepl("penguin", assembler) ~ "PenguiN",
      TRUE ~ NA_character_ # Ensure all cases are covered
    )
  ) |>
  # Check for NA values in cov and dmg
  filter(!is.na(cov) & !is.na(dmg)) |>
  mutate(
    cov = factor(cov, levels = c("3X", "5X", "10X")),
    dmg = factor(dmg, levels = c("high")),
    feature = fct_reorder(feature, length)
  ) |>
  ggplot(aes(x = assm, y = length, fill = feature)) +
  geom_col(alpha = 0.8, color = "black") +
  labs(
    x = "Assembler",
    y = "Base pairs"
  ) +
  facet_grid(dmg ~ cov, labeller = "label_both", scales = "free") +  # Use facet_grid instead of facet_grid2
  theme_bw(base_size = 14) +  # Adjust base_size to change the overall font size
  theme(
    strip.background = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12)  # Adjust size as needed
  ) +
  scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3)) +
  scale_fill_manual(values = c("#2D2D2D", "#EB554A", "#FFC300", "#F09868", "#91AEB7"), name = NULL)
