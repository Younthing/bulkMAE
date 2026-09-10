#' Look up the 0.4 name for a pre-0.4 function
#'
#' Version 0.4 removed the old exported aliases so autocomplete shows one API.
#' Use this helper when an old script still calls a 0.3 name.
#'
#' @param old_name A single pre-0.4 function name, with or without `()`.
#'   Use `NULL` to return the full map.
#'
#' @return A replacement string such as `score_gsva()` or
#'   `enrich_go(method = "ora")`. If `old_name` is `NULL`, a data frame with
#'   `old_name` and `new_name` columns.
#'
#' @export
#' @examples
#' bulkmae_rename("run_gsva")
#' bulkmae_rename("run_go_ora")
bulkmae_rename <- function(old_name = NULL) {
  map <- .bulkmae_rename_table()
  if (is.null(old_name)) {
    return(map)
  }
  if (!is.character(old_name) || length(old_name) != 1L || is.na(old_name)) {
    stop("`old_name` must be a single function name.", call. = FALSE)
  }
  key <- .bulkmae_rename_key(old_name)
  matched <- .bulkmae_rename_key(map$old_name) == key
  if (any(matched)) {
    return(map$new_name[matched][[1L]])
  }
  new_hit <- .bulkmae_rename_key(map$new_name) == key
  if (any(new_hit)) {
    old <- map$old_name[new_hit][[1L]]
    stop(
      "`",
      old_name,
      "` is already a 0.4 name (replacement for `",
      old,
      "`).",
      call. = FALSE
    )
  }
  stop(
    "`",
    old_name,
    "` is not a recorded pre-0.4 name. See vignette(\"naming-migration\", \"bulkMAE\").",
    call. = FALSE
  )
}

.bulkmae_rename_key <- function(name) {
  name <- gsub("^`+|`+$", "", name)
  sub("\\(.*$", "", name)
}

.bulkmae_rename_table <- function() {
  path <- system.file("guides", "naming-migration-zh.md", package = "bulkMAE")
  if (!nzchar(path) || !file.exists(path)) {
    path <- file.path("inst", "guides", "naming-migration-zh.md")
  }
  lines <- readLines(path, encoding = "UTF-8")
  old_name <- character()
  new_name <- character()
  for (line in lines) {
    if (!startsWith(line, "|")) {
      next
    }
    cells <- trimws(strsplit(line, "|", fixed = TRUE)[[1L]])
    cells <- cells[nzchar(cells)]
    if (length(cells) < 3L || grepl("0\\.3", cells[[2L]])) {
      next
    }
    if (grepl("^-{3,}$", cells[[1L]])) {
      next
    }
    old <- gsub("^`+|`+$", "", cells[[2L]])
    new <- gsub("^`+|`+$", "", cells[[3L]])
    if (!nzchar(old) || !nzchar(new)) {
      next
    }
    old_name <- c(old_name, old)
    new_name <- c(new_name, new)
  }
  data.frame(old_name = old_name, new_name = new_name, stringsAsFactors = FALSE)
}
