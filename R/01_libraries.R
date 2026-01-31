#===============================================================================
# Library Loading
#===============================================================================

required_packages <- c(
  "tidyverse",
  "nephro",
  "mcr",
  "santoku",
  "rms",
  "yardstick",
  "rsample",
  "marginaleffects",
  "modelsummary",
  "logistf",
  "flextable",
  "gtsummary",
  "gt",
  "glue",
  "scales",
  "Hmisc",
  "rootSolve",
  "table1"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    ". Please install them before running this script."
  )
}

message(
  "Loading required packages: ",
  paste(required_packages, collapse = ", ")
)
invisible(lapply(required_packages, library, character.only = TRUE))
