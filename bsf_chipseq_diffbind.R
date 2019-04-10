#!/usr/bin/env Rscript
#
# Copyright 2013 - 2019 Michael K. Schuster
#
# Biomedical Sequencing Facility (BSF), part of the genomics core facility
# of the Research Center for Molecular Medicine (CeMM) of the
# Austrian Academy of Sciences and the Medical University of Vienna (MUW).
#
#
# This file is part of BSF R.
#
# BSF R is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BSF R is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with BSF R.  If not, see <http://www.gnu.org/licenses/>.

suppressPackageStartupMessages(expr = library(package = "optparse"))

argument_list <- parse_args(object = OptionParser(
  option_list = list(
    make_option(
      opt_str = c("-v", "--verbose"),
      action = "store_true",
      default = TRUE,
      help = "Print extra output [default]",
      type = "logical"
    ),
    make_option(
      opt_str = c("-q", "--quietly"),
      action = "store_false",
      default = FALSE,
      dest = "verbose",
      help = "Print little output",
      type = "logical"
    ),
    make_option(
      opt_str = c("--comparison"),
      dest = "comparison",
      help = "Comparison name",
      type = "character"
    ),
    make_option(
      opt_str = c("--factor"),
      dest = "factor",
      help = "ChIP factor",
      type = "character"
    ),
    make_option(
      opt_str = c("--sample-annotation"),
      dest = "sample_annotation",
      help = "Sample annotation sheet",
      type = "character"
    ),
    make_option(
      opt_str = c("--threads"),
      default = 1L,
      dest = "threads",
      help = "Number of parallel processing threads [1]",
      type = "integer"
    )
  )
))

# Start of main script ----------------------------------------------------


# suppressPackageStartupMessages(expr = library(package = "BiocParallel"))
suppressPackageStartupMessages(expr = library(package = "DiffBind"))

prefix <-
  paste("chipseq", "diff", "bind", argument_list$comparison, argument_list$factor, sep = "_")

output_directory <- prefix

if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

# Initialise a DBA object -------------------------------------------------


diffbind_dba <- NULL
file_path <- file.path(prefix, paste0(prefix, '_dba.Rdata'))
if (file.exists(file_path) &&
    file.info(file_path)$size > 0L) {
  message("Loading a DiffBind DBA object")
  load(file = file_path)
} else {
  # Create a DBA object ---------------------------------------------------
  message("Creating a DiffBind DBA object")
  diffbind_dba <-
    dba(sampleSheet = argument_list$sample_annotation,
        bCorPlot = FALSE)
  
  # Count via the BiocParallel package that can be controlled more easily.
  # Unfortunately, this does not work because DBA_PARALLEL_BIOC does not seem to be exported.
  # diffbind_dba$config$parallelPackage <- DBA_PARALLEL_BIOC
  # Set the number of parallel threads in the MulticoreParam instance.
  # BiocParallel::register(BPPARAM = MulticoreParam(workers = argument_list$threads))
  
  # Plot heatmap on peak caller scores ------------------------------------
  message("Plotting a correlation heatmap based on peak caller score data")
  grDevices::png(filename = file.path(
    prefix,
    paste(prefix, "correlation_peak_caller_score.png", sep = "_")
  ))
  return_value <-
    DiffBind::dba.plotHeatmap(DBA = diffbind_dba, margin = 25)
  base::invisible(x = grDevices::dev.off())
  
  # Count reads -----------------------------------------------------------
  message("Counting reads")
  diffbind_dba <-
    DiffBind::dba.count(DBA = diffbind_dba, bCorPlot = FALSE, bParallel = FALSE)
  
  # Plot heatmap on read counts -------------------------------------------
  message("Plotting a correlation heatmap based on read counts")
  grDevices::png(filename = file.path(
    prefix,
    paste(prefix, "correlation_read_counts.png", sep = "_")
  ))
  return_value <-
    DiffBind::dba.plotHeatmap(DBA = diffbind_dba, margin = 25)
  base::invisible(x = grDevices::dev.off())
  
  # Establish contrasts -----------------------------------------------------
  message("Establishing contrasts by tissue")
  # The categories default to DiffBind::DBA_TISSUE, DiffBind::DBA_FACTOR, DiffBind::DBA_CONDITION and DiffBind::DBA_TREATMENT.
  diffbind_dba <-
    DiffBind::dba.contrast(DBA = diffbind_dba, minMembers = 2)
  # Check if setting contrasts was successful. It may not be, if less than two replicates were available.
  if (is.null(x = diffbind_dba$contrasts)) {
    # Set the mask manually. For the moment this only works for two samples.
    if (nrow(x = diffbind_dba$samples) == 2) {
      message("In lack of replicates, setting contrasts on the basis of the first two conditions")
      diffbind_conditions <-
        unique(x = diffbind_dba$class[DiffBind::DBA_CONDITION,])
      diffbind_dba <- DiffBind::dba.contrast(
        DBA = diffbind_dba,
        group1 = DiffBind::dba.mask(
          DBA = diffbind_dba,
          attribute = DiffBind::DBA_CONDITION,
          value = diffbind_conditions[1]
        ),
        group2 = DiffBind::dba.mask(
          DBA = diffbind_dba,
          attribute = DiffBind::DBA_CONDITION,
          value = diffbind_conditions[2]
        ),
        name1 = diffbind_conditions[1],
        name2 = diffbind_conditions[2]
      )
      rm(diffbind_conditions)
    }
  }
  
  # Run differential binding affinity analysis ----------------------------
  message("Running differential binding affinity analysis")
  diffbind_dba <-
    DiffBind::dba.analyze(DBA = diffbind_dba, bCorPlot = FALSE)
  
  # Plot heatmap on differential binding affinity -------------------------
  message("Plotting correlation heatmap based on differential binding affinity analysis")
  grDevices::png(filename = file.path(prefix, paste(
    prefix, "correlation_analysis.png", sep = "_"
  )))
  return_value <-
    DiffBind::dba.plotHeatmap(DBA = diffbind_dba, margin = 25)
  base::invisible(x = grDevices::dev.off())
  
  # Save DBA object -------------------------------------------------------
  message("Saving DBA object to disk")
  save(diffbind_dba, file = file_path)
}
rm(file_path)

# Create score-based PCA plot ---------------------------------------------
# Create a PCA plot irrespective of contrasts on the basis of scores in the main binding matrix.
message(sprintf("Creating PCA plot for comparison %s and factor %s", argument_list$comparison, argument_list$factor))
grDevices::png(filename = file.path(prefix, paste(prefix, "pca_plot.png", sep = "_")))
DiffBind::dba.plotPCA(DBA = diffbind_dba, attributes = DiffBind::DBA_CONDITION)
base::invisible(x = dev.off())

process_per_contrast <-
  function(contrast, group1, group2, db_number) {
    # Process per row of a contrasts data frame obtained via DiffBind::dba.show()
    # contrast the row.names() string of the data frame indicating the contrast number
    # group1 Group1 value
    # group2 Group2 value
    
    # The working directory has been set to the output_directory.
    
    # Write differentially bound sites ------------------------------------
    # These are the significantly differentially bound sites on the basis of the
    # configured FDR threshold.
    if (db_number == 0L) {
      message(
        sprintf(
          "Skipping differentially bound sites for comparison %s, factor %s and contrast %s versus %s",
          argument_list$comparison,
          argument_list$factor,
          group1,
          group2
        )
      )
    } else {
      message(
        sprintf(
          "Writing differentially bound sites for comparison %s, factor %s and contrast %s versus %s to disk",
          argument_list$comparison,
          argument_list$factor,
          group1,
          group2
        )
      )
      # To annotate peaks as differentially bound or not, export all sites as a
      # GRanges object and write it as a BED file to disk. All sites can be obtained by
      # setting the FDR threshold (th) to 1.0.
      
      granges_object <- DiffBind::dba.report(
        DBA = diffbind_dba,
        contrast = as.integer(x = contrast),
        th = 1.0,
        bNormalized = TRUE,
        bCalled = TRUE,
        bCounts = TRUE,
        bCalledDetail = TRUE,
        DataType = DiffBind::DBA_DATA_GRANGES
      )
      # The GenomicRanges::GRanges object returned by DiffBind::dba.report() does not have a valid GenomeInfoDb::Seqinfo object assigned.
      # Since the GenomeInfoDb::genomeStyles() function provides only mapping information about chromosomes,
      # but not on extra-chromosomal contigs, a clean assignment is impossible. The GenomeInfoDb::seqlevels(x) <- value assignment
      # requires a named character vector with a (complete) mapping.
      
      # Set the BED score on the basis of the FDR value, scaled and centered to fit
      # UCSC Genome Browser conventions.
      # The score should be an integer and range from 0 (white) to 1000 (black).
      GenomicRanges::score(x = granges_object) <-
        1000L - as.integer(x = round(
          x = scale(
            x = granges_object$FDR,
            center = min(granges_object$FDR),
            scale = diff(x = range(granges_object$FDR))
          ) * 1000.0
        ))
      
      rtracklayer::export.bed(
        object = granges_object,
        con = sprintf("%s_peaks_%s__%s.bed", prefix, group1, group2)
      )
      rm(granges_object)
      
      # The report function is quite peculiar in that it insists on a prefix DBA_
      # for the file name.
      file_path <-
        sprintf("%s_%s_report_%s__%s", argument_list$comparison, argument_list$factor, group1, group2)
      tryCatch(
        expr = {
          base::invisible(
            x = DiffBind::dba.report(
              DBA = diffbind_dba,
              contrast = as.integer(x = contrast),
              bNormalized = TRUE,
              bCalled = TRUE,
              bCounts = TRUE,
              bCalledDetail = TRUE,
              file = file_path,
              DataType = DiffBind::DBA_DATA_FRAME
            )
          )
        },
        error = function(cond) {
          message("DiffBind::dba.report failed with message:\n", cond, appendLF = TRUE)
        }
      )
      
      # Link differentially bound sites -----------------------------------
      # Create a symbolic link from the akward report file name to standard file names,
      # used by this script.
      file_path <-
        sprintf("DBA_%s_%s_report_%s__%s.csv",
                argument_list$comparison,
                argument_list$factor,
                group1,
                group2)
      link_path <-
        sprintf("%s_report_%s__%s.csv",
                prefix,
                group1,
                group2)
      if (file.exists(file_path) && !file.exists(link_path)) {
        if (!file.symlink(from = file_path,
                          to = link_path)) {
          warning(paste0(
            "Encountered an error linking DBA report file: ",
            file_path
          ))
        }
      }
      rm(link_path, file_path)
    }
    
    # Create MA plot ------------------------------------------------------
    message(
      sprintf(
        "Creating MA plot for comparison %s, factor %s and contrast %s versus %s",
        argument_list$comparison,
        argument_list$factor,
        group1,
        group2
      )
    )
    grDevices::png(filename = sprintf("%s_ma_plot_%s__%s.png", prefix, group1, group2))
    DiffBind::dba.plotMA(
      DBA = diffbind_dba,
      bNormalized = TRUE,
      bXY = FALSE,
      # FALSE for a MA plot.
      contrast = as.integer(x = contrast)
    )
    base::invisible(x = grDevices::dev.off())
    
    # Create Scatter plot -------------------------------------------------
    message(
      sprintf(
        "Creating scatter plot for comparison %s, factor %s and contrast %s versus %s",
        argument_list$comparison,
        argument_list$factor,
        group1,
        group2
      )
    )
    grDevices::png(filename = sprintf("%s_scatter_plot_%s__%s.png", prefix, group1, group2))
    DiffBind::dba.plotMA(
      DBA = diffbind_dba,
      bNormalized = TRUE,
      bXY = TRUE,
      # TRUE for a scatter plot.
      contrast = as.integer(x = contrast)
    )
    base::invisible(x = grDevices::dev.off())
    
    # Create PCA plot -----------------------------------------------------
    # This PCA plot is based upon the differential binding affinity analysis for the contrast.
    if (db_number == 0L) {
      message(
        sprintf(
          "Skipping PCA plot for comparison %s, factor %s and contrast %s versus %s",
          argument_list$comparison,
          argument_list$factor,
          group1,
          group2
        )
      )
    } else {
      message(
        sprintf(
          "Creating PCA plot for comparison %s, factor %s and contrast %s versus %s",
          argument_list$comparison,
          argument_list$factor,
          group1,
          group2
        )
      )
      grDevices::png(filename = sprintf("%s_pca_plot_%s__%s.png", prefix, group1, group2))
      tryCatch(expr = {
        DiffBind::dba.plotPCA(
          DBA = diffbind_dba,
          attributes = DiffBind::DBA_CONDITION,
          contrast = as.integer(x = contrast)
        )
      },
      error = function(cond) {
        message("DiffBind::dba.plotPCA failed with message:\n", cond, "\n", appendLF = TRUE)
      })
      base::invisible(x = dev.off())
    }
    
    # Create Box plot -----------------------------------------------------
    if (db_number == 0L) {
      message(
        sprintf(
          "Skipping Box plot for comparison %s, factor %s and contrast %s versus %s",
          argument_list$comparison,
          argument_list$factor,
          group1,
          group2
        )
      )
    } else {
      message(
        sprintf(
          "Creating Box plot for comparison %s, factor %s and contrast %s versus %s",
          argument_list$comparison,
          argument_list$factor,
          group1,
          group2
        )
      )
      grDevices::png(filename = sprintf("%s_box_plot_%s__%s.png", prefix, group1, group2))
      tryCatch(expr = {
        DiffBind::dba.plotBox(
          DBA = diffbind_dba,
          bNormalized = TRUE,
          contrast = as.integer(x = contrast)
        )
      },
      error = function(cond) {
        message("DiffBind::dba.plotBox failed with message:\n", cond, "\n", appendLF = TRUE)
      })
      base::invisible(x = grDevices::dev.off())
    }
    
    return(TRUE)
  }

# Get a data frame with all contrasts to apply the above function to each row.
contrast_frame <-
  DiffBind::dba.show(DBA = diffbind_dba, bContrasts = TRUE)

# Replace '!' characters with 'not_'.
contrast_frame$Group1 <-
  gsub(pattern = "!",
       replacement = "not_",
       x = contrast_frame$Group1)
contrast_frame$Group2 <-
  gsub(pattern = "!",
       replacement = "not_",
       x = contrast_frame$Group2)

# Write the contrasts data frame to disk.
write.csv(x = contrast_frame,
          file = file.path(prefix, paste(prefix, "contrasts.csv", sep = "_")),
          row.names = FALSE)

# Since DiffBind is quite peculiar in writing report files,
# set the output directory as new working directory. Sigh.
original_directory <- setwd(dir = output_directory)

# Filter out rows that have no significantly differentially occupied binding sites.
# FIXME: The analysis in the 5th column depends on the algorithm used.
# It seems available from DBA$config$AnalysisMethod.
# Apply the function to each row and discard the return value.
# print(x = "Contrast frame:")
# print(x = head(x = contrast_frame))
# print(x = str(object = contrast_frame))
# print(x = "db_number:")
# print(x = as.integer(x = as.character(x = contrast_frame[, 5L])))
return_value <-
  mapply(
    FUN = process_per_contrast,
    row.names(contrast_frame),
    contrast_frame$Group1,
    contrast_frame$Group2,
    # Since column 5 is a factor it needs converting into a character,
    # before converting into an integer.
    as.integer(x = as.character(x = contrast_frame[, 5L]))
  )
rm(return_value, contrast_frame)

# DiffBind::dba.overlap(DBA = diffbind_dba, mode = DiffBind::DBA_OLAP_RATE)

# message("Creating a Venn diagram")
# grDevices::png(filename = paste(prefix, "box_plot.png", sep = "_"))
# DiffBind::dba.plotVenn(DBA = diffbind_dba)
# base::invisible(x = grDevices::dev.off())

rm(
  diffbind_dba,
  process_per_contrast,
  prefix,
  output_directory,
  original_directory,
  argument_list
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())
