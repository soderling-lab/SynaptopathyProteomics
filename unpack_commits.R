#!/usr/bin/env Rscript

# 1. Seach git commit history for a file:
# $ git log --all --full-history -- "**/thefile.*"

# 2. Identify the correct commit by searching throught their files:
# $ git diff-tree --no-commit-id --name-only -r <commit_id>"

# 3. Once you find the file, you can get it with:
# git checkout <SHA>^ -- <path-to-file>

# Loop through commit history.
git_hist <- readLines("history.txt")

# Get commits. Remove "lazy git commit"
idx <- grepl("commit",git_hist)
out <- grepl("lazy git commit",git_hist)
commits <- git_hist[idx & !out]

# Get commit dates.
dates <- git_hist[grep("Date",git_hist)]
date_list <- lapply(strsplit(dates,":"),"[",c(2,3,4))
commit_dates <- sapply(date_list,function(x) trimws(paste(x,collapse=" ")))

# Get commit ids.
ids <- sapply(strsplit(commits,"\\ "),"[",2)

# Get files associated with each commit.
cmd <- "git diff-tree --no-commit-id --name-only -r"
myfiles <- lapply(ids,function(x) system(paste(cmd,x),intern=TRUE))

# 4th commit may contain my file.
mypath <- "code/3_WPCNA/ARCHIVE/wgcna-optimization.py"
ids[4]
id <- "082e71a30cf10f028273fd170348d5bee4492c4f"

# git checkout 082e71a30cf10f028273fd170348d5bee4492c4f -- code/3_WPCNA/ARCHIVE/wgcna-optimization.py

