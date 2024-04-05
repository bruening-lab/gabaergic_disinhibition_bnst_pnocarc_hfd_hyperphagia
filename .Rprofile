# Check if we are in vscode
is_interactive_and_vscode <-
  interactive() && Sys.getenv("RSTUDIO") == ""
if (is_interactive_and_vscode) {
  # Set vscode as the terminal program
  # This is needed to make it work within tmux sessions
  Sys.setenv(TERM_PROGRAM = "vscode")
}
# Source renv activate script
source("renv/activate.R")
# rVisidata Options
options(rvisidata.tmux = FALSE)
# Radian
options(radian.complete_while_typing = FALSE)
# Attach to vscode
if (is_interactive_and_vscode) {
  source(file.path(Sys.getenv(if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME"), ".vscode-R", "init.R"))
}
