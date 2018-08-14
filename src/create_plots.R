# Function that creates and saves plots from fastqc data
create_plots <- function(df_list) {
  p1 <- df_list$adapter_content %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    gather(key,
           value, -c(ref,
                     Position,
                     group,
                     seqlen)) %>%
    ggplot(aes(factor(
      Position,
      levels = unique(Position),
      ordered = TRUE
    ), value, color = key)) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(outlier.size = 0.5) +
    labs(
      x = "Position in Read",
      y = "Percent (%) Adapter Content",
      color = NULL,
      title = "Adapter content"
    ) +
    scale_colour_jama() +
    scale_y_continuous(limits = c(0, 20)) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom"
    ) +
    facet_wrap(~ seqlen, scales = "free")
  
  p2 <- df_list$per_base_sequence_content %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    gather(key, value, -c(ref, Base, group, seqlen)) %>%
    ggplot(aes(factor(
      Base, levels = unique(Base), ordered = TRUE
    ), value, color = key)) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(outlier.size = 0.5) +
    labs(
      x = "Position in Read",
      y = "Percent (%)",
      color = NULL,
      title = "Per base sequence content"
    ) +
    theme_classic() +
    scale_color_jama() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom"
    ) +
    facet_wrap( ~ seqlen, scales = "free")
  
  p3 <- df_list$sequence_duplication_levels %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    gather(key, value, -c(ref, `Duplication Level`, group, seqlen)) %>%
    ggplot(aes(
      factor(
        `Duplication Level`,
        levels = unique(`Duplication Level`),
        ordered = TRUE
      ),
      value,
      fill = key
    )) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(values = c("#ef8a62",
                                 "#67a9cf")) +
    theme_classic() +
    labs(x = "Duplication Level",
         y = "Percent (%) of Sequences",
         fill = NULL,
         title = "Sequence duplication levels") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.4
          )) +
    facet_wrap( ~ seqlen, scales = "free")
  
  # p4 <- df_list %>%
  #   prepare_seq_len_data() %>%
  #   left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
  #   rename(seqlen = "Sequence length") %>%
  #   ggplot(aes(factor(Length,
  #                     ordered = TRUE,
  #                     levels = unique(Length)),
  #              Count)) +
  #   geom_boxplot(outlier.size = 0.5) +
  #   scale_y_continuous(labels = comma) +
  #   scale_fill_viridis(discrete = TRUE) +
  #   labs(x = "Read Size",
  #        y = "Total sequence length",
  #        title = "Sequence length per read size") +
  #   guides(fill = FALSE) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(
  #     angle = 90,
  #     hjust = 1,
  #     vjust = 0.4
  #   )) +
  #   facet_wrap(~seqlen, scales = "free")
  
  p5 <- df_list$per_sequence_quality_scores %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    ggplot(aes(factor(Quality), Count, fill = factor(Quality))) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(outlier.size = 0.5) +
    scale_y_continuous(labels = comma) +
    scale_fill_viridis(discrete = TRUE) +
    labs(x = "Quality",
         y = "Number of reads",
         title = "Per sequence quality scores") +
    guides(fill = FALSE) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 7)) +
    facet_wrap( ~ seqlen)
  
  p6 <- df_list$per_sequence_gc_content %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    ggplot(aes(factor(`GC Content`), Count)) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(outlier.size = 0.5) +
    labs(x = "GC content (%)",
         y = "Number of reads",
         title = "Per sequence GC content") +
    scale_y_continuous(labels = comma) +
    scale_x_discrete(breaks = as.character(seq(
      from = 0, to = 100, by = 10
    ))) +
    theme_classic() +
    facet_wrap( ~ seqlen, scales = "free")
  
  p7 <- df_list$per_base_n_content %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    ggplot(aes(factor(
      Base, levels = unique(Base), ordered = TRUE
    ), `N-Count`)) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(fill = "#e6e6e6",
                 outlier.size = 0.5) +
    labs(x = "Position in read",
         title = "Per base N content") +
    guides(fill = FALSE) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    facet_wrap( ~ seqlen, scales = "free")
  
  p8 <- df_list$total_deduplicated_percentage %>%
    select(-group) %>%
    gather(key = "ref", value = `1`) %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
    rename(seqlen = "Sequence length",
           perc = `1`) %>%
    mutate(perc = as.numeric(perc)) %>%
    ggplot(aes(factor(seqlen), perc)) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(fill = "#e6e6e6",
                 outlier.size = 0.5) + 
    scale_y_continuous(limits = c(0, 100)) +
    labs(x = "Read size",
         y = "Total percentage of deduplicated reads",
         title = "Total deduplicated percentage") +
    theme_classic()
  
  # This produces a huge figure - excluded for now
  # p9 <- df_list$per_tile_sequence_quality %>%
  #   left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
  #   rename(seqlen = "Sequence length") %>%
  #   ggplot(aes(factor(Base,
  #                     levels = unique(Base),
  #                     ordered = TRUE),
  #              factor(Tile),
  #              fill = Mean)) +
  #   geom_tile() +
  #   scale_fill_viridis(limits = c(-10, 5),
  #                      breaks = c(-10, -8, -6, -4, -2, 0, 2, 4),
  #                      guide=guide_colourbar(ticks=T,nbin=50,barheight=.5,label=T,barwidth=10)) +
  #   labs(y = "Tiles",
  #        x = "Position in read",
  #        title = "Per tile sequence quality (mean)",
  #        fill = NULL) +
  #   theme_classic() +
  #   theme(legend.position="bottom",
  #         legend.justification="center",
  #         legend.direction="horizontal",
  #         legend.text=element_text(color="grey20"),
  #         axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         axis.text.y = element_text(size = 6)) +
  #   facet_wrap(~seqlen, scales = "free")
  
  p10 <- df_list$per_base_sequence_quality %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    mutate(Base = factor(Base,
                         levels = unique(Base),
                         ordered = TRUE)) %>%
    group_by(seqlen) %>%
    mutate(xmax = length(unique(Base)) + 1) %>%
    ungroup() %>%
    ggplot(aes(Base,
               Mean)) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(outlier.size = 0.4,
                 fill = "#7f7f7f") +
    geom_rect(aes(
      ymin = 28,
      ymax = Inf,
      xmin = 0,
      xmax = xmax
    ),
    fill = "#008000",
    alpha = 0.006) +
    geom_rect(aes(
      ymin = 20,
      ymax = 28,
      xmin = 0,
      xmax = xmax
    ),
    fill = "#FFDF00",
    alpha = 0.006) +
    geom_rect(aes(
      ymin = 0,
      ymax = 20,
      xmin = 0,
      xmax = xmax
    ),
    fill = "#ff1919",
    alpha = 0.006) +
    labs(x = "Position in read",
         y = "Sequence quality",
         title = "Per base mean sequence quality") +
    scale_y_continuous(limits = c(0, 42)) +
    theme_classic() +
    theme(axis.text.x = element_text(
      size = 7,
      angle = 90,
      hjust = 1,
      vjust = 0.4
    )) +
    facet_wrap(~ seqlen, scales = "free")
  
  save_plots(p1, "adapter_content", 25, 35)
  save_plots(p2, "per_base_sequence_content", 25, 35)
  save_plots(p3, "sequence_duplication_levels", 25, 30)
  save_plots(p5, "per_sequence_quality_scores", 25, 35)
  save_plots(p6, "per_sequence_gc_content", 25, 30)
  save_plots(p7, "per_base_n_content", 25, 35)
  save_plots(p8, "total_deduplicated_percentage", 20, 20)
  save_plots(p10, "per_base_mean_sequence_quality", 25, 25)
}