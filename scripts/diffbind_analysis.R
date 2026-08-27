#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(DiffBind); library(rtracklayer)})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 9L) stop("Usage: diffbind_analysis.R <manifest> <cohort> <consensus> <output_root> <outdir> <min_members> <alpha> <block_columns_csv> <subtract_control:true|false>")
manifest_file <- args[[1]]; cohort_id <- args[[2]]; consensus_file <- args[[3]]; output_root <- args[[4]]
outdir <- args[[5]]; min_members <- as.integer(args[[6]]); alpha <- as.numeric(args[[7]])
block_text <- args[[8]]; subtract_control <- tolower(args[[9]]) == "true"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
metadata <- read.delim(manifest_file, stringsAsFactors=FALSE, check.names=FALSE)
targets <- metadata[metadata$cohort_id == cohort_id & metadata$is_control == "FALSE", , drop=FALSE]
if (!nrow(targets)) stop("cohort absent")
condition_counts <- table(targets$condition)
eligible <- names(condition_counts[condition_counts >= min_members])
if (length(eligible) < 2L) {writeLines('{"status":"SKIPPED","reason":"insufficient replicated conditions"}',file.path(outdir,"SKIPPED.json"));quit(save="no",status=0)}
targets <- targets[targets$condition %in% eligible,,drop=FALSE]
bam_path <- function(key) file.path(output_root,"03_alignment","analysis",paste0(key,".host.analysis.bam"))
sheet <- data.frame(SampleID=targets$sample_key,Tissue=targets$cell_type,Factor=targets$factor,
    Condition=targets$condition,Treatment=targets$treatment,Replicate=as.integer(targets$replicate),
    bamReads=vapply(targets$sample_key,bam_path,character(1)),
    Peaks=consensus_file,PeakCaller="bed",stringsAsFactors=FALSE)
has_controls <- all(targets$control_key != "." & nzchar(targets$control_key))
if (subtract_control && !has_controls) {
    writeLines('{"status":"SKIPPED","reason":"matched controls unavailable"}',file.path(outdir,"SKIPPED.json"))
    quit(save="no",status=0)
}
if (has_controls) {
    sheet$bamControl <- vapply(targets$control_key,bam_path,character(1))
    sheet$ControlID <- targets$control_key
}
write.csv(sheet,file.path(outdir,"diffbind_samplesheet.csv"),row.names=FALSE,quote=TRUE)
db <- dba(sampleSheet=sheet)
db <- dba.blacklist(db, blacklist=FALSE, greylist=has_controls)
db <- dba.count(db, peaks=import(consensus_file), bSubControl=subtract_control,
                bScaleControl=subtract_control, bUseSummarizeOverlaps=TRUE)
blocks <- character(); if(block_text != "." && nzchar(block_text)) blocks <- trimws(strsplit(block_text,",",fixed=TRUE)[[1]])
if(length(blocks)) {writeLines('{"status":"SKIPPED","reason":"DiffBind v0.1 does not map arbitrary block columns; use primary DESeq2Enrichment"}',file.path(outdir,"SKIPPED.json"));quit(save="no",status=0)}
design <- "~Condition"
db <- dba.contrast(db, categories=DBA_CONDITION, minMembers=min_members)
db <- dba.analyze(db, method=DBA_DESEQ2)
contrast_metadata <- as.data.frame(dba.show(db,bContrasts=TRUE),stringsAsFactors=FALSE)
contrast_count <- nrow(contrast_metadata)
label_value <- function(row,candidates,fallback) {
    for(candidate in candidates) if(candidate %in% names(row)) {
        value <- as.character(row[[candidate]][[1]])
        if(!is.na(value) && nzchar(value)) return(value)
    }
    fallback
}
safe_label <- function(value) gsub("[^A-Za-z0-9._-]+","_",value)
comparison_rows <- list()
for (index in seq_len(contrast_count)) {
    metadata_row <- contrast_metadata[index,,drop=FALSE]
    numerator <- label_value(metadata_row,c("Group1","Name1","Condition1"),paste0("group1_contrast_",index))
    reference <- label_value(metadata_row,c("Group2","Name2","Condition2"),paste0("group2_contrast_",index))
    comparison_id <- paste0(safe_label(numerator),"_vs_",safe_label(reference))
    result <- dba.report(db, contrast=index, method=DBA_DESEQ2, th=1)
    result_table <- as.data.frame(result)
    all_name <- paste0("contrast_",index,"_all.tsv")
    significant_name <- paste0("contrast_",index,"_significant.tsv")
    write.table(result_table,file.path(outdir,all_name),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
    significant <- dba.report(db, contrast=index, method=DBA_DESEQ2, th=alpha)
    significant_table <- as.data.frame(significant)
    write.table(significant_table,file.path(outdir,significant_name),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
    if("Fold" %in% names(significant_table)) {
        higher_in_numerator <- sum(significant_table$Fold > 0,na.rm=TRUE)
        higher_in_reference <- sum(significant_table$Fold < 0,na.rm=TRUE)
    } else {
        higher_in_numerator <- NA_integer_; higher_in_reference <- NA_integer_
    }
    comparison_rows[[index]] <- data.frame(comparison_id=comparison_id,numerator=numerator,reference=reference,
        tested=nrow(result_table),significant=nrow(significant_table),higher_in_numerator=higher_in_numerator,
        higher_in_reference=higher_in_reference,status="SUCCESS",all_results=all_name,
        significant_results=significant_name,stringsAsFactors=FALSE)
}
if(length(comparison_rows)) write.table(do.call(rbind,comparison_rows),file.path(outdir,"comparison_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
saveRDS(db,file.path(outdir,"diffbind_object.rds"))
writeLines(c(paste("status: SUCCESS"),paste("design:",design),paste("subtract_control:",subtract_control),paste("contrasts:",contrast_count)),file.path(outdir,"summary.txt"))
writeLines(capture.output(sessionInfo()),file.path(outdir,"session_info.txt"))
