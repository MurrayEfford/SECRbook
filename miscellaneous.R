# utility functions for markdown

# from RMarkdown cookbook Yihui Xie, Christophe Dervieux, Emily Riederer
# 2023-12-30 Chapman & Hall/CRC

colorize <- function(x, color) {
 if (knitr::is_latex_output()) {
  sprintf("\\textcolor{%s}{%s}", color, x)
 } else if (knitr::is_html_output()) {
  sprintf("<span style='color: %s;'>%s</span>", color,
          x)
 } else x
}
