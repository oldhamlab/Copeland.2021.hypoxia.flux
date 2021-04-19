# utilities.R

src <- function() {
  files <-
    list.files(
    system.file("src/", package = "Copeland.2021.hypoxia.flux"),
    pattern = "\\.R$",
    full.names = TRUE
  )

  invisible(lapply(files, source))
}
