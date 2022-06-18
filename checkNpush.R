#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("methods")) # Rscript for CMD doesn't load this automatically
suppressPackageStartupMessages(library("optparse"))

# ---- inputs ----

parser <- OptionParser()
option_list <- list( 
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false", 
              dest="verbose", help="Print little output"),
  make_option(c("-m", "--message"), type="character", default='Minor fixes',
              help="Message for git commit -m [default %default]")
#   make_option("--generator", default="rnorm",
#               help = "Function to generate random deviates [default \"%default\"]"),
#   make_option("--mean", default=0,
#               help="Mean if generator == \"rnorm\" [default %default]"),
#   make_option("--sd", default=1, metavar="standard deviation",
#               help="Standard deviation if generator == \"rnorm\" [default %default]"),
#   make_option(c("-o", "--out_dir"), type = "character",
#               help="Directory for output files")
)
# wd = getwd()
inputs = commandArgs(trailingOnly=FALSE)
script_dir = dirname(sub("--file=", "", inputs[grep("--file=",inputs)]))
if(!length(script_dir)){
  script_dir = getwd()
}

# get real args
parser = OptionParser(
  usage =
    "%prog [options] <PACKAGE_DIR>
Description: 
  This script aims for the goal to check and build the package and push to github.",
  option_list = option_list
) # donot change the format of this doc.

itfs = parse_args(parser, positional_arguments = TRUE)
opts = itfs$options
args = list(
  input_dir = itfs$args[1]
)

# if none arguments input
if(is.na(args$input_dir)){
  args$input_dir = getwd()
}


setwd(args$input_dir)
# ---- loading dependencies ----
suppressPackageStartupMessages(library("devtools"))
suppressPackageStartupMessages(library("roxygen2"))
suppressPackageStartupMessages(library("usethis"))

# ---- main ----
roxygen2::roxygenize()
devtools::document()
devtools::build_readme()

# devtools::load_all()
# devtools::check()

# system('git add .')
# system(paste0("git commit -m '", opts$message, "'"))


