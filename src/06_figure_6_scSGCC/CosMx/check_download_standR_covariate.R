# Check and download standR_covariate repo if missing
# Target path used by CosMXNew scripts:
# /bmbl_data/yuzhou/collaborative/Sizun_lab/INDEPTH/SGCC/SGWT_results/DLBCL_GeoMX/standR_covariate

repo_url <- "https://github.com/BMEngineeR/standR_covariate"
target_dir <- "/bmbl_data/yuzhou/collaborative/Sizun_lab/INDEPTH/SGCC/SGWT_results/DLBCL_GeoMX/standR_covariate"

check_git <- function() {
  nzchar(Sys.which("git"))
}

clone_repo <- function(url, dest) {
  if (dir.exists(dest)) {
    message("Target exists; skipping clone: ", dest)
    return(invisible(TRUE))
  }

  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

  if (check_git()) {
    cmd <- sprintf("git clone %s %s", shQuote(url), shQuote(dest))
    message("Cloning with git: ", cmd)
    status <- system(cmd)
    if (status != 0) stop("git clone failed.")
    return(invisible(TRUE))
  }

  # Fallback: download zip and extract
  zip_url <- paste0(url, "/archive/refs/heads/main.zip")
  tmp_zip <- tempfile(fileext = ".zip")
  message("Downloading zip: ", zip_url)
  utils::download.file(zip_url, tmp_zip, mode = "wb", quiet = FALSE)

  tmp_dir <- tempfile("standR_covariate_")
  dir.create(tmp_dir)
  utils::unzip(tmp_zip, exdir = tmp_dir)

  # Locate extracted folder (repo-main)
  extracted <- list.dirs(tmp_dir, full.names = TRUE, recursive = FALSE)
  if (length(extracted) == 0) stop("No extracted folder found in zip.")

  file.rename(extracted[1], dest)
  unlink(tmp_zip)
  unlink(tmp_dir, recursive = TRUE)
  invisible(TRUE)
}

if (!dir.exists(target_dir)) {
  message("standR_covariate not found at: ", target_dir)
  clone_repo(repo_url, target_dir)
} else {
  message("standR_covariate already present: ", target_dir)
}

# Optional: show repo status if git is available
if (dir.exists(target_dir) && check_git()) {
  message("Repo status:")
  system(sprintf("git -C %s status -s", shQuote(target_dir)))
}
