#
# Library module of bsfR GRC to UCSC genome mapping functions.
#
# Copyright 2013 - 2020 Michael K. Schuster
#
# Biomedical Sequencing Facility (BSF), part of the genomics core facility of
# the Research Center for Molecular Medicine (CeMM) of the Austrian Academy of
# Sciences and the Medical University of Vienna (MUW).
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

#' Get a genome information tibble.
#'
#' @param resource_directory A \code{character} scalar with the resource (i.e.
#'   genome and transcriptome) directory path.
#' @param ensembl_version A \code{character} scalar with the Ensembl version.
#'
#' @return A \code{tibble} with genome information.
#' @export
#' @importFrom rlang .data .env :=
#'
#' @examples
#' \dontrun{
#'  genome_tibble <-
#'    bsfR::bsfg_get_genome_tibble(
#'      resource_directory = "/path/to/resources/genomes/",
#'      ensembl_version = "96"
#'    )
#' }
bsfg_get_genome_tibble <-
  function(resource_directory, ensembl_version) {
    # Start with the bsfg_genome_tibble packaged in the sysdata.R file.
    genome_tibble <- bsfg_genome_tibble

    # Set genome and transcriptome directory names
    for (file_type in c("ncbi", "ucsc")) {
      # "genome_path_ncbi" = "resource_directory/GRCh38"
      # "genome_path_ucsc" = "resource_directory/hg38"

      # "transcriptome_path_ncbi" = "resource_directory/GRCh38_e87"
      # "transcriptome_path_ucsc" = "resource_directory/hg38_e87"

      assembly_version <-
        paste("assembly", "version", file_type, sep = "_")

      genome_path <- paste("genome", "path", file_type, sep = "_")

      transcriptome_path <-
        paste("transcriptome", "path", file_type, sep = "_")

      genome_tibble <-
        dplyr::mutate(
          .data = genome_tibble,
          !!genome_path := file.path(.env$resource_directory,
                                     .data[[assembly_version]]),
          !!transcriptome_path := file.path(.env$resource_directory,
                                            paste(
                                              .data[[assembly_version]],
                                              paste0("e", .env$ensembl_version),
                                              sep = "_"
                                            ))
        )

      rm(transcriptome_path, genome_path, assembly_version)
    }
    rm(file_type)

    # Set assembly report paths for NCBI and UCSC directories.
    genome_tibble <- dplyr::mutate(
      .data = genome_tibble,
      # Set the assembly report name, which is the last element of the URL split
      # by "/" characters.
      "assembly_report_name" = purrr::map_chr(
        .x = stringr::str_split(
          string = .data$assembly_report_url,
          pattern = stringr::fixed(pattern = "/")
        ),
        .f = ~ .[length(x = .)]
      ),
      # Set assembly report paths for NCBI directories.
      "assembly_report_path_ncbi" = file.path(
        .env$resource_directory,
        .data$assembly_version_ncbi,
        .data$assembly_report_name
      ),
      # Set assembly report paths for UCSC directories.
      "assembly_report_path_ucsc" = file.path(
        .env$resource_directory,
        .data$assembly_version_ucsc,
        .data$assembly_report_name
      ),
      # Set a species prefix.
      "species_prefix" = base::gsub(
        pattern = " ",
        replacement = "_",
        x = .data$scientific_name
      )
    )

    return(genome_tibble)
  }

#' Get a genome information list.
#'
#' @param resource_directory A \code{character} scalar with the resource (i.e.
#'   genome and transcriptome) directory path.
#' @param ensembl_version A \code{character} scalar with the Ensembl version.
#' @param ncbi_version A \code{character} vector matched via the "in" operator.
#' @param ucsc_version A \code{character} vector matched via the "in" operator.
#'
#' @return A \code{list} with genome information.
#' @export
#' @importFrom rlang .data .env
#'
#' @examples
#' \dontrun{
#'  genome_list <-
#'    bsfR::bsfg_get_genome_list(
#'      resource_directory = "/path/to/resources/genomes/",
#'      ensembl_version = "96",
#'      ncbi_version = "GRCh38",
#'      ucsc_version = "hg38"
#'    )
#' }
bsfg_get_genome_list <-
  function(resource_directory,
           ensembl_version,
           ncbi_version,
           ucsc_version) {
    # The dplyr::filter() function automatically discards cases that evaluate to NA.
    genome_tibble <- dplyr::filter(
      .data = bsfR::bsfg_get_genome_tibble(resource_directory = resource_directory, ensembl_version = ensembl_version),
      .data$assembly_version_ncbi %in% .env$ncbi_version &
        .data$assembly_version_ucsc %in% .env$ucsc_version
    )
    if (nrow(x = genome_tibble) != 1L) {
      # This should only retrieve one row of the tibble.
      stop("Could not get a single record for source and target assembly version.")
    }
    genome_list <- as.list(genome_tibble)
    rm(genome_tibble)

    return(genome_list)
  }

#' Create NCBI and UCSC "genome" or "transcriptome" resource directories.
#'
#' @param genome_list A \code{list} of genome assembly annotation.
#' @param biological_type A \code{character} scalar of "genome" or "transcriptome".
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  bsfR::bsfg_create_resource_directory(
#'    genome_list = genome_list,
#'    biological_type = "genome",
#'    verbose = verbose
#'  )
#' }
bsfg_create_resource_directory <-
  function(genome_list, biological_type, verbose = FALSE) {
    stopifnot(biological_type %in% c("genome", "transcriptome"))

    # Create NCBI and UCSC-style genome and transcriptome directories, if they
    # do not exist already and if the assembly version is defined.
    for (assembly_type in c("ncbi", "ucsc")) {
      # for (biological_type in c("genome", "transcriptome")) {
      # Check if the assembly_version is not NA.
      if (is.na(x = genome_list[[paste("assembly_version", assembly_type, sep = "_")]])) {
        next()
      }
      file_path <-
        genome_list[[paste(biological_type, "path", assembly_type, sep = "_")]]
      if (!file.exists(file_path)) {
        if (verbose) {
          message("Creating directory: ", file_path)
        }
        dir.create(path = file_path,
                   showWarnings = TRUE,
                   recursive = FALSE)
      }
      rm(file_path)
      # }
      # rm(biological_type)
    }
    rm(assembly_type)

    return()
  }

#' Read an NCBI assembly report from the NCBI FTP site.
#'
#' The assembly report provides a conversion table from NCBI-style to UCSC-style
#' sequence region names and gets imported into a \code{tibble}. Adds in
#' conversion information for the mitochondrion, which is not part of the
#' assembly report.
#'
#' @param genome_list A \code{list} of genome assembly annotation.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{tibble} with NCBI assembly information.
#' @export
#'
#' @examples
#' \dontrun{
#'  assembly_report_tibble <-
#'    bsfR::bsfg_get_assembly_report(
#'      genome_list,
#'      verbose = FALSE
#'    )
#' }
bsfg_get_assembly_report <- function(genome_list, verbose = FALSE) {
  # Create NCBI and UCSC-style genome and transcriptome directories, if they do
  # not exist already and if the assembly version is defined.
  bsfR::bsfg_create_resource_directory(
    genome_list = genome_list,
    biological_type = "genome",
    verbose = verbose
  )

  if (!file.exists(genome_list$assembly_report_path_ucsc)) {
    if (verbose) {
      message(
        "Downloading assembly report: ",
        genome_list$assembly_report_url,
        " ... to destination: ",
        genome_list$assembly_report_path_ucsc
      )
    }

    if (utils::download.file(
      url = genome_list$assembly_report_url,
      destfile = genome_list$assembly_report_path_ucsc
    ) > 0L) {
      stop("Failed to download: ", genome_list$assembly_report_url)
    }
  }

  if (verbose) {
    message("Reading assembly report: ",
            genome_list$assembly_report_path_ucsc)
  }

  assembly_tibble <-
    readr::read_tsv(
      file = genome_list$assembly_report_path_ucsc,
      col_names = c(
        "sequence_name",
        "sequence_role",
        "assigned_molecule",
        "assigned_molecule_location_type",
        "genbank_accn",
        "relationship",
        "refseq_accn",
        "assembly_unit",
        "sequence_length",
        "ucsc_style_name"
      ),
      col_types = readr::cols(
        sequence_name = readr::col_character(),
        sequence_role = readr::col_character(),
        assigned_molecule = readr::col_character(),
        assigned_molecule_location_type = readr::col_character(),
        genbank_accn = readr::col_character(),
        relationship = readr::col_character(),
        refseq_accn = readr::col_character(),
        assembly_unit = readr::col_character(),
        sequence_length = readr::col_integer(),
        ucsc_style_name = readr::col_character()
      ),
      comment = "#"
    )

  # The NCBI assembly report does not contain information about the
  # mitochondrion.
  if (!("chrM" %in% assembly_tibble$ucsc_style_name)) {
    assembly_tibble <-
      tibble::add_row(
        .data = assembly_tibble,
        "ucsc_style_name" = "chrM",
        "sequence_name" = "MT"
      )
  }

  return(assembly_tibble)
}

#' Get an Ensembl transcriptome annotation file from the Ensembl FTP site.
#'
#' @param genome_list A \code{list} of genome assembly annotation.
#' @param ensembl_version A \code{character} scalar with the Ensembl version
#'   number.
#' @param ensembl_ftp A \code{character} scalar with the Ensembl FTP URL prefix.
#'   Defaults to "ftp://ftp.ensembl.org".
#' @param file_type A \code{character} scalar with the file type i.e. "gff3" of
#'   "gtf".
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{character} scalar with the Ensembl transcriptome annotation
#'   file path.
#' @export
#'
#' @examples
#' \dontrun{
#'  ensembl_transcriptome_gtf <-
#'    bsfR::bsfg_get_ensembl_transcriptome(
#'      genome_list = genome_list,
#'      ensembl_version = "100",
#'      ensembl_ftp = "ftp://ftp.ensembl.org",
#'      file_type = "gtf",
#'      verbose = FALSE
#'    )
#' }
bsfg_get_ensembl_transcriptome <-
  function(genome_list,
           ensembl_version,
           ensembl_ftp = "ftp://ftp.ensembl.org",
           file_type = "gtf",
           verbose = FALSE) {
    stopifnot(file_type %in% c("gff3", "gtf"))

    # Create NCBI and UCSC-style genome and transcriptome directories, if they
    # do not exist already and if the assembly version is defined.
    bsfR::bsfg_create_resource_directory(
      genome_list = genome_list,
      biological_type = "transcriptome",
      verbose = verbose
    )

    # Source name (e.g. Homo_sapiens.GRCh38.87.gff3.gz or Homo_sapiens.GRCh38.87.gtf.gz)
    file_name <- paste(
      genome_list$species_prefix,
      genome_list$assembly_version_ncbi,
      ensembl_version,
      file_type,
      "gz",
      sep = "."
    )
    file_path <-
      file.path(genome_list$transcriptome_path_ncbi, file_name)

    if (!file.exists(file_path)) {
      # Source path (e.g. ftp://ftp.ensembl.org/pub/release-87/gff3/homo_sapiens/Homo_sapiens.GRCh38.87.gff3.gz)
      # Source path (e.g. ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz)
      # If the source file does not exist, download it first.
      file_url <- paste(
        ensembl_ftp,
        "pub",
        paste("release", ensembl_version, sep = "-"),
        file_type,
        tolower(x = genome_list$species_prefix),
        file_name,
        sep = "/"
      )
      if (verbose) {
        message(
          "Downloading Ensembl transcriptome: ",
          file_url,
          " ... to destination: ",
          file_path
        )
      }
      if (utils::download.file(url = file_url, destfile = file_path) < 0L) {
        stop("Failed to download: ", file_url)
      }
      rm(file_url)
    }
    rm(file_name)

    return(file_path)
  }

#' Convert sequence levels between NCBI and UCSC.
#'
#' Either a NCBI or UCSC \code{GenomicRanges::GRanges} object needs specifying. Based on the
#' preselected genome list the corresponding NCBI assembly report will be loaded.
#' The mapping procedure is complicated by the fact that UCSC-style names map to
#' either NCBI sequence names (e.g. 1, 2, HSCHR1_CTG1_UNLOCALIZED, ...) or NCBI
#' GenBank accession numbers (e.g. CM000663.2, CM000664.2, KI270706.1, ...).
#' Therefore, mapping orientation-specific maps need building depending on the
#' NCBI names that are used.
#'
#' @param genome_list A \code{list} of genome assembly annotation.
#' @param source_granges A \code{GenomicRanges::GRanges} object.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{GenomicRanges::GRanges} object with converted sequence levels.
#' @export
#'
#' @examples
#' \dontrun{
#'  target_granges <-
#'    bsfR::bsfg_convert_seqlevels(
#'      genome_list,
#'      source_granges = source_granges,
#'      verbose = FALSE
#'    )
#' }
bsfg_convert_seqlevels <-
  function(genome_list,
           source_granges,
           verbose = FALSE) {
    target_granges <- NULL
    ncbi_granges <- NULL
    ucsc_granges <- NULL
    ncbi_seqinfo <- NULL
    ucsc_seqinfo <- NULL

    # Is the source GenomicRanges::GRanges object associated with a GenomeInfoDb::Seqinfo object?
    source_seqinfo <- GenomicRanges::seqinfo(x = source_granges)
    if (!is.null(x = source_seqinfo)) {
      # Get a GenomeInfoDb::Seqinfo object for NCBI.
      if (!is.na(x = genome_list$bsgenome_ncbi)) {
        ncbi_seqinfo <-
          rtracklayer::SeqinfoForBSGenome(genome = genome_list$bsgenome_ncbi)
      }
      # Get a GenomeInfoDb::Seqinfo object for UCSC.
      if (!is.na(x = genome_list$bsgenome_ucsc)) {
        ucsc_seqinfo <-
          rtracklayer::SeqinfoForBSGenome(genome = genome_list$bsgenome_ucsc)
      }

      if (GenomeInfoDb::genome(x = source_seqinfo)[1L] == GenomeInfoDb::genome(x = ncbi_seqinfo)[1L]) {
        ncbi_granges <- source_granges
        ncbi_seqinfo <- source_seqinfo
      }
      if (GenomeInfoDb::genome(x = source_seqinfo)[1L] == GenomeInfoDb::genome(x = ucsc_seqinfo)[1L]) {
        ucsc_granges <- source_granges
        ucsc_seqinfo <- source_seqinfo
      }
    }

    if (is.null(x = ncbi_seqinfo)) {
      stop("Could not retrieve a NCBI Seqinfo object.")
    }
    if (verbose) {
      print(x = paste(
        "NCBI sequence levels:",
        length(x = GenomeInfoDb::seqlevels(x = ncbi_seqinfo))
      ))
      print(x = GenomeInfoDb::seqlevels(x = ncbi_seqinfo))
    }

    if (is.null(x = ucsc_seqinfo)) {
      stop("Could not retrieve a UCSC Seqinfo object.")
    }
    if (verbose) {
      print(x = paste(
        "UCSC sequence levels:",
        length(x = GenomeInfoDb::seqlevels(x = ucsc_seqinfo))
      ))
      print(x = GenomeInfoDb::seqlevels(x = ucsc_seqinfo))
    }

    assembly_report_tibble <-
      bsfR::bsfg_get_assembly_report(genome_list = genome_list, verbose = verbose)

    if (!is.null(x = ucsc_granges)) {
      # Map from UCSC to NCBI sequence levels.
      #
      # Create a vector of UCSC sequence levels named by NCBI sequence levels.
      # Initialise the vector with as many NA values as the NCBI
      # GenomeInfoDb::Seqinfo object has seqlevels. Then map to NCBI assembly
      # report sequence names and to NCBI assembly report GenBank accession
      # numbers.

      ucsc_levels <-
        base::rep(NA_character_, times = length(x = GenomeInfoDb::seqlevels(x = ncbi_seqinfo)))
      names(x = ucsc_levels) <-
        GenomeInfoDb::seqlevels(x = ncbi_seqinfo)

      ucsc_names <-
        assembly_report_tibble$ucsc_style_name

      # Test for NCBI sequence names.
      # (e.g. 1, 2, HSCHR1_CTG1_UNLOCALIZED, ...)
      names(x = ucsc_names) <-
        assembly_report_tibble$sequence_name
      ucsc_map <-
        ucsc_names[names(x = ucsc_levels)]

      # Merge only those sequence names into the UCSC seqlevels that could be
      # resolved.
      ucsc_levels[!is.na(x = ucsc_map)] <-
        ucsc_map[!is.na(x = ucsc_map)]

      # Test for NCBI GenBank accession numbers.
      # (e.g. CM000663.2, CM000664.2, KI270706.1, ...)
      names(x = ucsc_names) <-
        assembly_report_tibble$genbank_accn
      ucsc_map <-
        ucsc_names[names(x = ucsc_levels)]

      # Merge only those GenBank accessions into the UCSC seqlevels that could
      # be resolved.
      ucsc_levels[!is.na(x = ucsc_map)] <-
        ucsc_map[!is.na(x = ucsc_map)]

      # Replace the remaining UCSC NA values with their NCBI names.
      ucsc_levels[is.na(x = ucsc_levels)] <-
        names(x = ucsc_levels[is.na(x = ucsc_levels)])

      rm(ucsc_map, ucsc_names)
      if (verbose) {
        print(x = paste(
          "UCSC sequence levels with NCBI names:",
          length(x = ucsc_levels)
        ))
        print(x = ucsc_levels)
      }

      match_integer <-
        base::match(x = ucsc_levels,
                    table = GenomeInfoDb::seqlevels(x = ucsc_granges))
      if (verbose) {
        message("Map UCSC to NCBI sequence levels")
        # Print the integer match vector.
        print(x = paste("Match integer vector:", length(x = match_integer)))
        print(x = match_integer)
        # Create a character match vector.
        match_character <-
          ucsc_levels
        names(x = match_character) <-
          GenomeInfoDb::seqlevels(x = ucsc_granges)[match_integer]
        print(x = paste("Match character vector:", length(x = match_character)))
        print(x = match_character)
        rm(match_character)
      }
      target_granges <- ucsc_granges
      GenomicRanges::seqinfo(x = target_granges, new2old = match_integer) <-
        ncbi_seqinfo
      rm(match_integer, ucsc_levels)
    }

    if (!is.null(x = ncbi_granges)) {
      # Map from NCBI to UCSC sequence levels.
      #
      # Create a vector of NCBI sequence levels named by the UCSC
      # GenomeInfoDb::Seqinfo object. Initialise the vector with as many NA
      # values as the UCSC GenomeInfoDb::Seqinfo object has seqlevels. Then map
      # to NCBI assembly report sequence names and to NCBI assembly report
      # GenBank accession numbers.

      ncbi_levels <-
        base::rep(NA_character_, times = length(x = GenomeInfoDb::seqlevels(x = ucsc_seqinfo)))
      names(x = ncbi_levels) <-
        GenomeInfoDb::seqlevels(x = ucsc_seqinfo)
      # if (verbose) {
      #   print(x = paste("NCBI levels:", length(x = ncbi_levels)))
      #   print(x = ncbi_levels)
      # }

      # Test for NCBI assembly report sequence names.
      # (e.g. 1, 2, HSCHR1_CTG1_UNLOCALIZED, ...)
      ncbi_names <-
        assembly_report_tibble$sequence_name
      names(x = ncbi_names) <-
        assembly_report_tibble$ucsc_style_name
      ncbi_map <-
        ncbi_names[names(x = ncbi_levels)]
      # if (verbose) {
      #   print(x = paste("NCBI map sequence names:", length(ncbi_map)))
      #   print(x = ncbi_map)
      # }

      # Merge only those NCBI sequence names that match to NCBI
      # GenomicRanges::GRanges.
      ncbi_levels[ncbi_map %in% GenomeInfoDb::seqlevels(x = ncbi_granges)] <-
        ncbi_map[ncbi_map %in% GenomeInfoDb::seqlevels(x = ncbi_granges)]
      # if (verbose) {
      #   print(x = paste("NCBI levels:", length(x = ncbi_levels)))
      #   print(x = ncbi_levels)
      # }

      # Test for NCBI assembly report GenBank accession numbers.
      # (e.g. CM000663.2, CM000664.2, KI270706.1, ...)
      ncbi_names <-
        assembly_report_tibble$genbank_accn
      names(x = ncbi_names) <-
        assembly_report_tibble$ucsc_style_name
      ncbi_map <-
        ncbi_names[names(x = ncbi_levels)]
      # if (verbose) {
      #   print(x = paste("NCBI map GenBank accessions:", length(x = ncbi_map)))
      #   print(x = ncbi_map)
      # }

      # Merge only those NCBI GenBank accession that match to NCBI
      # GenomicRanges::GRanges.
      ncbi_levels[ncbi_map %in% GenomeInfoDb::seqlevels(x = ncbi_granges)] <-
        ncbi_map[ncbi_map %in% GenomeInfoDb::seqlevels(x = ncbi_granges)]
      # if (verbose) {
      #   print(x = paste("NCBI levels:", length(x = ncbi_levels)))
      #   print(x = ncbi_levels)
      # }

      # Replace the remaining NCBI NA values with their UCSC names.
      ncbi_levels[is.na(x = ncbi_levels)] <-
        names(x = ncbi_levels[is.na(x = ncbi_levels)])
      rm(ncbi_map, ncbi_names)

      if (verbose) {
        print(x = paste(
          "NCBI sequence levels with UCSC names:",
          length(x = ncbi_levels)
        ))
        print(x = ncbi_levels)
      }

      match_integer <-
        base::match(x = ncbi_levels,
                    table = GenomeInfoDb::seqlevels(x = ncbi_granges))
      if (verbose) {
        message("Map NCBI to UCSC sequence levels")
        print(x = paste("Match integer vector:", length(x = match_integer)))
        print(x = match_integer)
        match_character <-
          ncbi_levels
        names(x = match_character) <-
          GenomeInfoDb::seqlevels(x = ncbi_granges)[match_integer]
        print(x = paste("Match character vector:", length(x = match_character)))
        print(x = match_character)
        rm(match_character)
      }
      target_granges <- ncbi_granges
      GenomicRanges::seqinfo(x = target_granges, new2old = match_integer) <-
        ucsc_seqinfo
      rm(match_integer)
    }
    rm(assembly_report_tibble,
       ncbi_granges,
       ncbi_seqinfo,
       ucsc_granges,
       ucsc_seqinfo)

    return(target_granges)
  }
