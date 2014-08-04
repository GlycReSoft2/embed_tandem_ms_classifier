#!/usr/bin/Rscript
library("optparse")

option_list <- list(
    make_option(c("-g","--gold-standard"), action = "store", help = "Gold Standard Dataset to train the model on"),
    make_option(c("-t","--target"), action = "store", help = "Target Dataset to classify"),
    make_option(c("-o", "--output"), action = "store", help = "Name of output file (defaults to <target>-scored.csv")
)

parser <- OptionParser(usage = "usage: %prog -g <gold-standard> -t <target> [-o <output>]", option_list = option_list)

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

suppressPackageStartupMessages(source(file.path(base_path, "label-and-disambiguate.R")))

if(!('target' %in% names(args)) || !("gold-standard" %in% names(args))){
    stop("Must provide both gold-standard and target datasets")
}

gold_standard <- prepareAnnotatedModel(args[['gold-standard']])

target <- prepareModel(args[["target"]])

model_fit <- fitModel(gold_standard)

target$MS2_Score <- predict(model_fit, target, "prob")


labeled_target <- labelAmbiguity(target)

output_file <- args$output

if(is.null(output_file)){
    library("tools")
    output_file <- paste(file_path_sans_ext(args[["target"]]), ".scored.csv", sep = "")
}

write.csv(labeled_target, output_file)
cat(output_file)