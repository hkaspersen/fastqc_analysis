# Function that saves plots
save_plots <- function(plot, title, height, width) {
  ggsave(paste0(output_dir, "/", title, ".svg"),
         plot,
         dpi = 100,
         device = "svg",
         units = "cm",
         height = height,
         width = width)
}