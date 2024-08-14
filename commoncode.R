# code run for each chapter with R chunks

# It is not strictly necessary to run this in every chapter after the first, 
# but it makes debugging easier, as each chapter may be knitted independently

library(secr)
setNumThreads(18)
hareCH6 <- read.capthist("data/hareCH6capt.txt", "data/hareCH6trap.txt", detector = "single")
options(digits = 5, width = 70)  # 85       

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
