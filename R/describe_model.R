#!/usr/bin/Rscript
library("optparse")
library("base64enc")
option_list <- list(
    make_option(c("-t","--target"), action = "store", help = "Target Dataset to classify"),
    make_option(c("-o", "--output"), action = "store", help = "Name of output file (defaults to <target>-described.json")
)

parser <- OptionParser(usage = "usage: %prog -t <target> [-o <output>]", option_list = option_list)

args <- parse_args(parser)

# http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
interpolate_path <- function(){
    initial.options <- commandArgs(trailingOnly = FALSE)
    file.arg.name <- "--file="
    script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
    script.basename <- dirname(script.name)
    return(script.basename)
}

base_path <- interpolate_path()

suppressPackageStartupMessages(source(file.path(base_path, "label_and_disambiguate.R")))

model_data <- prepareAnnotatedModel(args$target)
model_forest <- fitModel(model_data)

temp_image <- tempfile()

png(filename)
require(randomForest)
varImpPlot(model_forest)
dev.off()

img_data <- base64encode(temp_image)

