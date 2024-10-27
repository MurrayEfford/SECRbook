# from Nick Tierney

# https://www.njtierney.com/post/2022/04/11/rmd-to-qmd/

library(fs)
library(stringr)
.etc

# move chunk options to body
qmd_names <- dir_ls(path = ".", glob = "*.qmd")
sapply(qmd_names, knitr::convert_chunk_header, output=identity)


# edit .Rprofile to include sourcing commoncode.R


knitr::convert_chunk_header("16-troubleshooting.qmd", output=identity)
