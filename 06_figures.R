library(tidyverse)
library(glue)
library(patchwork)

# Frame 0 reads comparison ====================================================
samples_ordered <- c(
    "activatedegg_chx", "activatedegg_nochx", "embryo_early", "embryo_mid",
    "embryo_late", "embryo_2h", "embryo_6h", "embryo_12h", "embryo_24h")

logfiles <- list.files('data/psite_logs', full.names = TRUE)
frametest <- read_delim(logfiles, skip = 28, n_max = 3, col_names = FALSE, delim = '  ')

frametest <- frametest |> 
    rename(frame = X1, None = X2, Peak = X3, GBT = X4) |> 
    mutate(lib = rep(sub('.log', '', basename(logfiles)), each = 3)) |> 
    separate(lib, into = c('sample', 'param'), sep = '_rpf.') |> 
    filter(frame == 0) |>
    select(-frame) |> 
    pivot_wider(names_from = 'param', values_from = GBT) |>
    pivot_longer(!sample, values_drop_na = TRUE)  |> 
    distinct(sample, name, .keep_all = TRUE)


plt <- frametest |>
    filter(name != 'psite') |> 
    mutate(name = if_else(name == 'om10', 'PSite', name),
           name = factor(name, levels = c('None', 'Peak', 'PSite')),
           sample = factor(sample, levels = rev(samples_ordered)))

p1 <- ggplot(plt, aes(y = sample, x = value, fill = name)) +
    geom_col(position = 'dodge') +
    geom_vline(xintercept = 0.5, linetype = 2) +
    labs(y = NULL, x = 'Fraction assigned to frame 0', fill = 'Method') +
    scale_fill_brewer(palette = 'Blues') +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_classic() +
    theme(axis.text = element_text(color = 'black'),
          legend.position = c(0.99, 0.01),
          legend.justification = c(1, 0))
ggsave('figures/PSite_peak_none_comparison.pdf', p1, width = 5, height = 5)

# translated ORF prediction comparison ========================================
dmel_txinfo <- read_tsv('data/Drosophila_melanogaster.BDGP6.32.52.gtf.txinfo')

types_fine <- c('CDS', 'N_extension', 'N_truncation', 'uORF', 'uoORF', 'dORF', 'doORF', 'iORF', 'lncRNA-ORF')
types_abnormal <- c('iCDS', 'C_truncation', 'sCDS', 'wORF', 'ncRNA-ORF', 'pseudogene-ORF')

read_orfs <- function(path){
  orfs <- read_tsv(path) |> 
    select(-c(orf_type_ori, gene_biotype, phase_start, phase_end, protein_id:mrna_end_nf)) |> 
    mutate(across(c(tstart, tend, gstart:cds_end), as.integer)) |>
    mutate(orf_type = if_else(
      orf_type == 'ncRNA-ORF' & str_starts(gene_name, 'lncRNA'),
      'lncRNA-ORF', orf_type)) |> 
    filter(orf_type %in% types_fine)
  
  orfs |> 
    mutate(orf_type = factor(orf_type, levels = types_fine),
           # uid_stop = glue('{chrom}:{gend}{strand}'),
           uuid = glue('{chrom}:{gstart}-{gend}{strand}')) |>
    arrange(uuid, orf_type) |>
    distinct(uuid, .keep_all = TRUE)
}

deft_files <- list.files('data/allorfs/default', pattern = 'processed.tsv', full.names = TRUE)
names(deft_files) <- sub('_processed.tsv', '', basename(deft_files))
deft_orfs <- lapply(deft_files, read_orfs)
deft_orfs <- bind_rows(deft_orfs, .id = 'sample')
deft_orfs <- deft_orfs |> separate(sample, c('sample', 'method'), sep = '_rpf.')

psite_files <- list.files('data/allorfs/psite', pattern = 'processed.tsv', full.names = TRUE)
psite_files <- psite_files[str_detect(psite_files, 'om10')]
names(psite_files) <- sub('_processed.tsv', '', basename(psite_files))
psite_orfs <- lapply(psite_files, read_orfs)
psite_orfs <- bind_rows(psite_orfs, .id = 'sample')
psite_orfs <- psite_orfs |> separate(sample, c('sample', 'method'), sep = '_rpf.om10.')

# comparisons 
all_orfs <- bind_rows(Default = deft_orfs, PSite = psite_orfs, .id = 'pipeline')
all_orfs <- all_orfs |>
    mutate(cds_variants = orf_type %in% c('CDS', 'N_extension', 'N_truncation'))

orfs_stat <- all_orfs |> count(pipeline, sample, method, orf_type, cds_variants)

plt <- orfs_stat |>
  group_by(pipeline, sample, method) |> 
  mutate(prop = n/sum(n)) |>
  pivot_wider(names_from = pipeline, values_from = c(n, prop)) |>
  filter(orf_type %in% c('CDS', 'N_extension', 'N_truncation')) |> 
  group_by(sample, method) |> 
  summarize(across(n_Default:prop_PSite, sum), .groups = 'drop')

  plt1 <- orfs_stat |>
    filter(method == 'ribocode') |>
    group_by(pipeline, sample) |> 
    summarise(total = sum(n),
              cds = sum(n[cds_variants]), .groups = 'drop') |> 
    pivot_longer(total:cds) |> 
    mutate(fgrp = str_c(pipeline, '-', if_else(name == 'cds', 'CDS', 'Novel')))

p2 <- ggplot(plt1, aes(y = factor(sample, levels = rev(samples_ordered)),
                x = value, group = pipeline, fill = fgrp)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    scale_fill_brewer(name = NULL, palette = 'Paired',
                      direction = -1) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = 'Number of translated ORFs', y = NULL) +
    theme_classic() +
    theme(axis.text = element_text(color = 'black'),
          legend.position = c(0.99, 0.99),
          legend.justification = c(1, 1))


plt2 <- orfs_stat |>
    filter(method == 'ribotish') |>
    group_by(pipeline, sample) |> 
    summarise(total = sum(n),
              cds = sum(n[cds_variants]), .groups = 'drop') |> 
    pivot_longer(total:cds) |> 
    mutate(fgrp = str_c(pipeline, '-', if_else(name == 'cds', 'CDS', 'Novel')))

p3 <- ggplot(plt2, aes(y = factor(sample, levels = rev(samples_ordered)),
                x = value, group = pipeline, fill = fgrp)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    scale_fill_brewer(name = NULL, palette = 'Paired',
                      direction = -1) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = 'Number of translated ORFs', y = NULL) +
    theme_classic() +
    theme(axis.text = element_text(color = 'black'),
          legend.position = c(0.99, 0.99),
          legend.justification = c(1, 1))

p1 + (p3 + p2 + plot_layout(guides = 'collect')) + plot_layout(widths = c(2, 5))
ggsave('figures/figure_1.pdf', width = 12, height = 3)

x <- plt2 |> filter(name == 'total') |> 
    pivot_wider(names_from = pipeline, values_from = value, id_cols = sample)
wilcox.test(x$Default, x$PSite, paired = TRUE)
