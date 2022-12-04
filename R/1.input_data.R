#' Test eps and minPts combinations for DBSCAN analysis
#'
#' This function tests all combinations of eps and minPts for DBSCAN analysis
#' of reference samples indicated in refID. The results are represented in
#' scatterplots exported to a pdf file.
#' @inheritParams dPCP
#' @param  refID a string or a character vector of chipID (Thermo Fisher) or
#'   the complete file name with the extension (Bio-Rad) of reference sample(s)
#'   to  be analysed.
#' @param  reference.quality numeric. Between 0 and 1. Quality threshold
#'   to subset the data (just for Thermo Fisher). If different thresholds have
#'   to be applied to various reference samples, a vectror of the same length
#'   of \code{refID} has to be provided.
#' @param  eps a numeric vector of values to be tested. Maximum distance
#'   between elements within a cluster in a DBSCAN analysis.
#'   See also \code{\link[dbscan]{dbscan}}.
#' @param  minPts a numeric vector of values to be tested. Number of minimum
#'   elements to assemble a cluster in a DBSCAN analysis.
#'   See also \code{\link[dbscan]{dbscan}}.
#' @return A pdf file containing the scatterplots of DBSCAN analysis performed
#' with all combinations of eps and minPts.
#' Each reference generates a different pdf file.
#' @examples
#' \donttest{
#' library(dPCP)
#'
#' #Find path of sample table and location of reference and input files
#' sampleTable <- system.file("extdata", "Template_sampleTable.csv",
#'                      package = "dPCP")
#'
#' fileLoc <- system.file("extdata", package = "dPCP")
#'
#' dbscan_combination("dilution20200313_B01_Amplitude.csv",
#'                    file.location = fileLoc, system = "bio-rad",
#'                    eps = c(150, 160, 180, 190), minPts = c(80, 100, 120))
#'
#' unlink("dilution20200313_B01_Amplitude.pdf")
#' }
#' @export

dbscan_combination <- function(refID, system = NULL, file.location = ".",
                               reference.quality = 0.5,
                               eps = c(120, 150, 180, 200),
                               minPts = c(20, 50, 80, 100)) {

  if (!is.character(refID))
    stop("'refID' must be a character string or vector")

  if (any(c("Thermo Fisher", "ThermoFisher", "thermo fisher", "thermofisher",
            "Thermo", "thermo", "t", "T") == system)) {
    system <- "Thermo Fisher"

  } else if (any(c("Bio-Rad", "BioRad", "Bio-rad", "Biorad", "bio-rad",
                   "biorad", "Bio", "bio", "B", "b") == system)) {
    system <- "Bio-Rad"
    reference.quality <- "Defined by Bio-Rad"

  } else if (any(c("Other", "other", "O", "o") == system)) {
    system <- "Other"
    reference.quality <- "Defined by the supplier"

  } else {stop("system must be either Thermo Fisher, Bio-Rad, or other")}

  if ((system == "Thermo Fisher") & (
    !is.numeric(reference.quality) || any(reference.quality > 1) ||
    any(reference.quality < 0)))
    stop("Invalid value for 'reference.quality'. It must be between 0 and 1")

  if ((system == "Thermo Fisher") & (
    length(reference.quality) > 1 &
    length(reference.quality) != length(unique(refID))))
    stop("Invalid value for 'reference.quality'. It must be a numeric value or
         a numeric vector of lenght equal to the number of reference samples")

  if (!is.numeric(c(eps, minPts))) stop("'eps' and 'minPts' must be numeric")

  #Calculate all combinations
  combinations <- expand.grid(eps, minPts)

  comb <- lapply(seq(length(refID)), function(x) {

    if (system == "Thermo Fisher") {
      #List reference experiment file
      eds.file <- list.files(file.location,
                             pattern = paste(refID[x], ".eds", sep = ""))

      if (length(eds.file) == 0)
        stop(paste0("An experiment file of reference", " '", refID[x],
                    "'", " can not be found."))

      if (length(eds.file) > 1)
        stop(paste0("Multiple experiment files for reference", " '",
                    refID[x], "'."))

      #Extract fluoresce and quality data
      quality.value <- utils::read.csv(
        unz(paste0(file.location, "/", eds.file),
            "apldbio/sds/quality_metrics/QVhOverallScores.CSV"))

      fluorescence <- utils::read.csv(
        unz(paste0(file.location, "/", eds.file),
            "apldbio/sds/quants_dye/ComponentData.CSV"))

      reference.data <- matrix(c(
        quality.value[, 1], fluorescence[, 3], fluorescence[, 1]), ncol = 3,
        dimnames = list(1:nrow(quality.value), c("Quality", "Vic", "Fam")))

      if (length(reference.quality) > 1) {
        reference.quality <- reference.quality[x]
      }

      #Keep data with quality score higher than reference.quality
      submatrix <- subset.matrix(
        reference.data[, c(2, 3)], reference.data[, 1] >= reference.quality)

      if (nrow(submatrix) < 10000)
        stop(paste0("The number of wells qualified by reference.quality is less
                    than 10000 for reference ", "'", refID[x], "'.", "
                    Try to reduce the quality threshold or run another chip"))

    } else {

      ref.file <- list.files(file.location, pattern = refID[x])

      if (length(ref.file) == 0)
        stop(paste0("No file of reference", " '", refID[x], "'",
                    " can be found."))

      if (length(ref.file) > 1)
        stop(paste0("Multiple files for reference", " '", refID[x], "'."))

      if (stringr::str_sub(ref.file, start = -4) == ".csv") {

        #If it is a `.csv` file, extract fluorescence and quality data
        fluorescence <- utils::read.csv(
          paste0(file.location, "/", ref.file), na.strings = c("NA", ""),
          colClasses = c("numeric", "numeric", "NULL"))

        submatrix <- matrix(c(fluorescence[,2],fluorescence[,1]), ncol = 2)
        colnames(submatrix)  <-  c("Vic", "Fam")
      }
    }



    #Try all combinations and plot the results
    graph.comb <- lapply(seq(nrow(combinations)), function(y) {

      dbclus <- dbscan::dbscan(submatrix, combinations[y, 1],
                               combinations[y, 2])

      p <- ggplot(as.data.frame(submatrix), aes_string(x = "Vic", y = "Fam")) +

        geom_point(size = 0.5, aes(color = factor(dbclus$cluster))) +

        geom_point(size = 0.5,
                   data = as.data.frame(submatrix)[dbclus$cluster == 0, ],
                   aes_string(x = "Vic", y = "Fam"), colour = "grey40") +

        labs(title = paste0("eps = ", dbclus$eps, ", minPts = ",
                            dbclus$minPts)) +

        theme(
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.position = "none",
          plot.title = element_text(size = 14),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text = element_text(size = 9)
        )
    })
    #Draw 4 plots per page and export to pdf
    tot <- ggpubr::ggarrange(plotlist = graph.comb, ncol = 2, nrow = 2)

    if (system == "Thermo Fisher") {
      ggpubr::ggexport(tot, filename = paste0(refID[x], ".pdf"))

    } else {
      ggpubr::ggexport(
        tot,
        filename = paste0(stringr::str_sub(refID[x], end = -5), ".pdf"))
    }
  })
}


#' Read sample table
#'
#' This function reads a file containing the essential information about the
#' samples and experimental settings. The file has to be filled out by the user
#' and formatted as described in the vignette.
#' @inheritParams dPCP
#' @return An object of class \code{sample_table}.
#' @examples
#' library(dPCP)
#'
#' #Find path of sample table and location of reference and input files
#' sampleTable <- system.file("extdata", "Template_sampleTable.csv",
#'                      package = "dPCP")
#'
#' fileLoc <- system.file("extdata", package = "dPCP")
#'
#' #Read sample table file
#' sample.table <- read_sampleTable(sampleTable, system = "bio-rad",
#'                                  file.location = fileLoc)
#' @export

read_sampleTable <- function(file, system = NULL, file.location = ".") {

  if (!is.character(file)) stop("'file' must be a character string")

  if (stringr::str_sub(file, start = -4) != ".csv")
    stop("'file' must be a '.csv' file. Please use the template provided.")


  if (any(c("Thermo Fisher", "ThermoFisher", "thermo fisher", "thermofisher",
            "Thermo", "thermo", "t", "T") == system)) {
    system <- "Thermo Fisher"

  } else if (any(c("Bio-Rad", "BioRad", "Bio-rad", "Biorad", "bio-rad",
                   "biorad", "Bio", "bio", "B", "b") == system)) {
    system <- "Bio-Rad"

  } else if (any(c("Other", "other", "O", "o") == system)) {
    system <- "Other"
    reference.quality <- "Defined by the supplier"

  } else {stop("system must be either Thermo Fisher, Bio-Rad, or other")}

  if (!is.character(file.location))
    stop("'file.location' must be a path name indicating reference and
         sample files location")


  sample.table <- utils::read.csv(file,
                                   stringsAsFactors = FALSE,
                                   na.strings = c("NA", ""))

  if (ncol(sample.table) != 9)
    stop("Columns number is not correct. Please use the template provided.")

  if (any(colnames(sample.table) !=
          c("Sample.name", "Chip.ID.Well.ID", "No.of.targets", "FAM.target",
            "Target.3", "Target.4", "VIC.HEX.target", "Reference",
            "Dilution")) ||
      any(is.na(colnames(sample.table))))
    stop("Columns names are not correct. Please use the template provided.")

  if (any(grepl(",", sample.table))) {
    new.col <- lapply(sample.table[, c(1, 4:7)], function(x) {
      x <- gsub(",", ".", x)
    })
    sample.table[, c(1, 4:7)] <- new.col
  }

  if (any(is.na(sample.table$Sample.name))) {
    na.names <- sapply(which(is.na(sample.table$Sample.name)), function(x) {
      paste0("Sample", x)
    })
    sample.table$Sample.name[which(is.na(sample.table$Sample.name))] <-
      na.names
  }

  if (any(is.na(sample.table$Chip.ID.Well.ID)))
    stop("Sample Chip or Well ID must be provided")

  if (!is.integer(sample.table$No.of.targets) ||
      any(is.na(sample.table$No.of.targets)) ||
      any(sample.table$No.of.targets > 4) ||
      any(sample.table$No.of.targets < 1))
    stop("Invalid value for No of targets")

  if (any(duplicated(sample.table$Sample.name)) &&
      any(duplicated(sample.table$Sample.name) !=
          duplicated(paste(sample.table$Sample.name,
                           sample.table$No.of.targets,"_"))))
    stop("Replicates must have the same number of targets")

  if (any(is.na(sample.table$Reference))) {

    ref <- lapply(which(is.na(sample.table$Reference)), function(x) {

      if (system == "Thermo Fisher") {
        ref.file <- sample.table$Chip.ID.Well.ID[x]
      } else {
        ref.file <- list.files(
          file.location,
          pattern = paste0(as.character(sample.table$Chip.ID.Well.ID[x]),
                           "_Amplitude.csv"))
      }

      if (length(ref.file) > 1) {

        stop(paste0("More than one file has been found for the sample ",
                    sample.table$Chip.ID.Well.ID[x]))

      } else if (length(ref.file) == 0)
        stop(paste0("No file has been found for the sample ",
                    sample.table$Chip.ID.Well.ID[x]))

      ref.file
    })

    sample.table$Reference[which(is.na(sample.table$Reference))] <-
      unlist(ref)
  }

  if (!is.numeric(sample.table$Dilution) ||
      any(is.na(sample.table$Dilution)) ||
      any(sample.table$Dilution > 1) ||
      any(sample.table$Dilution <= 0))
    stop("Invalid value for Dilution.")

  class(sample.table) <- "sample_table"
  return(sample.table)
}


#' Read reference files
#'
#' This function reads the results files of reference samples listed in the
#' sample table. Fluoresce intensity and quality value (just for Thermo Fisher)
#' are collected.
#' If a \code{\link{reference_dbscan}} template file with the same input
#' paramters (reference ID, eps, minPts) is available, fluorescence data,
#' quality value and dbscan analysis results are retrived from the template
#' file.
#' @param  sample.table object of class \code{sample_table}, inherited from
#'   \code{\link{read_sampleTable}}.
#' @param  eps,minPts numeric. Input parameters for the DBSCAN algorithm. If
#'   they match the paramters of \code{\link{reference_dbscan}} template file,
#'   the data are retrived from the template.
#' @inheritParams dPCP
#' @return An object of class \code{read_reference} containing a sublist for
#'   each reference. Each sublist has the following components:
#'   \item{quality}{value of the \code{reference.quality} parameter.}
#'   \item{data}{a matrix with the fluorescence intensities and quality
#'   values.}
#'   \item{dbscan}{an object of class \code{dbscan_fast}, inherited from
#'     \code{\link[dbscan]{dbscan}}. This component is available only if a
#'     \code{\link{reference_dbscan}} template file is used to retrive the
#'     data.}
#' @examples
#' library(dPCP)
#'
#' #Find path of sample table and location of reference and input files
#' sampleTable <- system.file("extdata", "Template_sampleTable.csv",
#'                      package = "dPCP")
#'
#' fileLoc <- system.file("extdata", package = "dPCP")
#'
#' #Read sample table file
#' sample.table <- read_sampleTable(sampleTable, system = "bio-rad",
#'                                  file.location = fileLoc)
#'
#' #Read reference files
#' ref <- read_reference(sample.table, system = "bio-rad",
#'                       file.location = fileLoc)
#' @export

read_reference <- function(sample.table, system = NULL, file.location = ".",
                           reference.quality = 0.5, eps = NULL,
                           minPts = NULL) {

  if (!inherits(sample.table, "sample_table"))
    stop("'sample.table' must be an object of class sample_table")


  if (any(c("Thermo Fisher", "ThermoFisher", "thermo fisher", "thermofisher",
            "Thermo", "thermo", "t", "T") == system)) {
    system <- "Thermo Fisher"

  } else if (any(c("Bio-Rad", "BioRad", "Bio-rad", "Biorad", "bio-rad",
                   "biorad", "Bio", "bio", "B", "b") == system)) {
    system <- "Bio-Rad"
    reference.quality <- "Defined by Bio-Rad"

  } else if (any(c("Other", "other", "O", "o") == system)) {
    system <- "Other"
    reference.quality <- "Defined by the supplier"

  } else {stop("system must be either Thermo Fisher, Bio-Rad, or other")}


  if ((system == "Thermo Fisher") & (
    !is.numeric(reference.quality) || any(reference.quality > 1) ||
    any(reference.quality < 0)))
    stop("Invalid value for 'reference.quality'. It must be between 0 and 1")

  if ((system == "Thermo Fisher") & (
    length(reference.quality) > 1 &
    length(reference.quality) != length(unique(sample.table$Reference))))
    stop("Invalid value for 'reference.quality'. It must be a numeric value or
         a numeric vector of lenght equal to the number of reference samples")

  if (!is.character(file.location))
    stop("'file.location' must be a path name indicating reference and
         sample files location")

  if (!is.null(eps) & !is.numeric(eps)) stop("'eps' must be numeric")

  if (!is.null(minPts) & !is.numeric(minPts)) stop("'minPts' must be numeric")

  if (length(eps) > 1 & length(eps) != length(unique(sample.table$Reference)))
    stop("Invalid value for 'eps'. It must be a numeric value or a numeric
         vector of lenght equal to the number of reference samples")

  if (length(minPts) > 1 & length(minPts) !=
      length(unique(sample.table$Reference)))
    stop("Invalid value for 'minPts'. It must be a numeric value or a numeric
         vector of lenght equal to the number of reference samples")

  if (system == "Thermo Fisher") {
    #List reference .eds files and DBSCAN reference files (.rts)
    reference.filesnames <- sapply(unique(sample.table$Reference), function(x){

      eds.files <- list.files(file.location,
                              pattern = paste0(x, ".eds"))

      if (!is.null(eps) & !is.null(minPts)) {

        if (length(eps) == 1) {
          new.eps <- eps
        } else {
          new.eps <- eps[match(x, unique(sample.table$Reference))]
        }

        if (length(minPts) == 1) {
          new.minPts <- minPts
        } else {
          new.minPts <- minPts[match(x, unique(sample.table$Reference))]
        }

        if (length(reference.quality) == 1) {
          DB.files <- list.files(
            file.location, pattern = paste0(x, "_", "qlty",
                                            reference.quality, "_",
                                            "eps", new.eps, "_",
                                            "minPts", new.minPts, "_DB",
                                            ".rds"))
        } else {
          DB.files <- list.files(
            file.location, pattern = paste0(
              x, "_", "qlty",
              reference.quality[match(x, unique(sample.table$Reference))],
              "_", "eps", new.eps, "_", "minPts", new.minPts, "_DB", ".rds"))
        }
      } else {
        DB.files <- NULL
      }

      if (length(eds.files) == 0 & length(DB.files) == 0)
        stop(paste0("Neither experiment file nor DBSCAN prediction file of
                    reference", " '", x, "'", " can be found."))

      if (length(eds.files) > 1)
        stop(paste0("Multiple experiment files for reference", " '", x, "'."))

      if (length(DB.files) > 1)
        stop(paste0("Multiple DBSCAN prediction files for reference", " '", x,
                    "'."))

      if (length(DB.files) == 1) {
        eds.files <- DB.files
      }

      if (length(eds.files) == 1 & length(DB.files) == 1) {
        print(paste0(
          "Both experiment and DBSCAN prediction files of reference", " '", x,
          "'", " have been found. DBSCAN prediction file will be used for the
          analysis."))}

      eds.files
    })

    #Extract data from .eds or .rts files
    subquality <- lapply(seq_along(reference.filesnames), function(y) {

      #If it is a `.rds` file, extract dbscan results
      if (stringr::str_sub(reference.filesnames[y], start = -4) == ".rds") {

        quality <- readRDS(paste0(file.location, "/",
                                  reference.filesnames[y]))$quality
        data <- readRDS(paste0(file.location, "/",
                               reference.filesnames[y]))$data
        submatrix <- readRDS(paste0(file.location, "/",
                                    reference.filesnames[y]))$dbscan

        refdata <- list("quality" = quality, "data" = data,
                        "dbscan" = submatrix)

      } else if (stringr::str_sub(reference.filesnames[y], start = -4) ==
                 ".eds") {

        #If it is a `.eds` file, extract fluorescence and quality data
        quality.value <- utils::read.csv(
          unz(paste0(file.location, "/", reference.filesnames[y]),
              "apldbio/sds/quality_metrics/QVhOverallScores.CSV"))

        fluorescence <- utils::read.csv(
          unz(paste0(file.location, "/", reference.filesnames[y]),
              "apldbio/sds/quants_dye/ComponentData.CSV"))

        reference.data <- matrix(
          c(quality.value[, 1], fluorescence[, 3], fluorescence[, 1]),
          ncol = 3, dimnames = list(1:nrow(quality.value),
                                    c("Quality", "Vic", "Fam")))

        #Keep data with quality score higher than reference.quality
        if (length(reference.quality) == 1) {
          submatrix <- subset.matrix(reference.data[, c(2, 3)],
                                     reference.data[, 1] >= reference.quality)
          refdata <- list("quality" = reference.quality, "data" = submatrix)
        } else {
          submatrix <- subset.matrix(
            reference.data[, c(2, 3)], reference.data[, 1] >=
              reference.quality[y])
          refdata <- list("quality" = reference.quality[y], "data" = submatrix)
        }
        if (nrow(submatrix) < 10000)
          stop(paste0("The number of wells qualified by 'reference.quality' is
                    less than 10000 for reference ", "'",
                      stringr::str_sub(reference.filesnames[y], start = -10,
                                       end = -5), "'.", "
                    Try to reduce the quality threshold or run another chip"))

        refdata
      }
    })
  } else {

    #Extract data from .csv or .rts files
    subquality <- lapply(unique(sample.table$Reference), function(y) {

      ref.file <- list.files(file.location, pattern = y)

      if (length(ref.file) == 0)
        stop(paste0("No file of reference", " '", y, "'", " can be found."))

      if (length(ref.file) > 1)
        stop(paste0("Multiple files for reference", " '", y, "'."))

      #If it is a `.rds` file, extract dbscan results
      if (stringr::str_sub(ref.file, start = -4) == ".rds") {

        quality <- readRDS(paste0(file.location, "/",
                                  ref.file))$quality
        data <- readRDS(paste0(file.location, "/",
                               ref.file))$data
        submatrix <- readRDS(paste0(file.location, "/",
                                    ref.file))$dbscan

        refdata <- list("quality" = quality, "data" = data,
                        "dbscan" = submatrix)

      } else if (stringr::str_sub(ref.file, start = -4) == ".csv") {

        #If it is a `.csv` file, extract fluorescence and quality data
        fluorescence <- utils::read.csv(
          paste0(file.location, "/", ref.file),
          na.strings = c("NA", ""),
          colClasses = c("numeric", "numeric"))

        data <- matrix(c(fluorescence[,2], fluorescence[,1]), ncol = 2)
        colnames(data)  <-  c("Vic", "Fam")

        refdata <- list("quality" = reference.quality, "data" = data)
      }
    })
  }

  names(subquality) <- unique(sample.table$Reference)

  class(subquality) <- "read_reference"

  return(subquality)
}


#' Read sample files
#'
#' This function reads the results files of samples listed in the sample table.
#' Fluoresce intensity and quality value (just for Thermo Fisher) are
#' collected.
#' @inheritParams read_reference
#' @inheritParams dPCP
#' @return An object of class \code{read_sample} containing a sublist for each
#'   sample. Each sublist has the following components:
#'   \item{quality}{value of the \code{sample.quality} parameter.}
#'   \item{data}{a matrix with the fluorescence intensities and quality
#'   values.}
#' @examples
#' library(dPCP)
#'
#' #Find path of sample table and location of reference and input files
#' sampleTable <- system.file("extdata", "Template_sampleTable.csv",
#'                      package = "dPCP")
#'
#' fileLoc <- system.file("extdata", package = "dPCP")
#'
#' #Read sample table file
#' sample.table <- read_sampleTable(sampleTable, system = "bio-rad",
#'                                  file.location = fileLoc)
#'
#' #Read reference files
#' ref <- read_reference(sample.table, system = "bio-rad",
#'                       file.location = fileLoc)
#'
#' #Read samples files
#' samp <- read_sample(sample.table, system = "bio-rad",
#'                     file.location = fileLoc)
#' @export

read_sample <- function(sample.table, system = NULL, file.location = ".",
                        sample.quality = 0.5, partition.volume = NULL) {

  if (!inherits(sample.table, "sample_table"))
    stop("'sample.table' must be an object of class sample_table")

  if (!is.character(file.location))
    stop("'file.location' must be a path name indicating reference and
         sample files location")

  if (any(c("Thermo Fisher", "ThermoFisher", "thermo fisher", "thermofisher",
            "Thermo", "thermo", "t", "T") == system)) {
    system <- "Thermo Fisher"

  } else if (any(c("Bio-Rad", "BioRad", "Bio-rad", "Biorad", "bio-rad",
                   "biorad", "Bio", "bio", "B", "b") == system)) {
    system <- "Bio-Rad"
    sample.quality <- "Defined by Bio-Rad"

  } else if (any(c("Other", "other", "O", "o") == system)) {

    if (is.numeric(partition.volume) & length(partition.volume) == 1) {
      system <- "Other"
      sample.quality <- paste0(
        "Defined by the supplier. Partition volume: ", partition.volume)
    } else {
      stop("partition.volume must be numeric")
    }

  } else {stop("system must be either Thermo Fisher, Bio-Rad, or other")}


  if ((system == "Thermo Fisher") & (
    !is.numeric(sample.quality) || any(sample.quality > 1) ||
    any(sample.quality < 0)))
    stop("Invalid value for 'sample.quality'. It must be between 0 and 1")

  if ((system == "Thermo Fisher") & (
    length(sample.quality) > 1 &
    length(sample.quality) != length(sample.table$Chip.ID.Well.ID)))
    stop("Invalid value for 'sample.quality'. It must be a numeric value or
         a numeric vector of lenght equal to the number of samples")

  #Assign name to the samples (samplename_chipID)
  name <- paste(sample.table$Sample.name, sample.table$Chip.ID, sep = "_")

  if (system == "Thermo Fisher") {
    #List .eds sample files
    filesnames <- sapply(sample.table$Chip.ID, function(x) {

      eds.files <- list.files(file.location,
                              pattern = paste0(as.character(x), ".eds"))

      if (length(eds.files) == 0)
        stop(paste0("Cannot find the file of sample ", x))

      if (length(eds.files) > 1)
        stop(paste0("Multiple files for sample", x))

      eds.files
    })

    #Extract quality and fluorescence data from .eds files
    subquality <- lapply(seq_along(filesnames), function(y) {

      quality <- utils::read.csv(
        unz(paste0(file.location, "/", filesnames[[y]]),
            "apldbio/sds/quality_metrics/QVhOverallScores.CSV"))

      fluorescence <- utils::read.csv(
        unz(paste0(file.location, "/",filesnames[[y]]),
            "apldbio/sds/quants_dye/ComponentData.CSV"))

      alldata <- matrix(
        c(quality[, 1], fluorescence[, 3], fluorescence[, 1]), ncol = 3,
        dimnames = list(1:nrow(quality), c("Quality", "Vic", "Fam")))

      #Keep data with quality score higher than reference.quality
      if (length(sample.quality) > 1) {
        sample.quality <- sample.quality[y]
      }

      submatrix <- subset.matrix(
        alldata[, c(2, 3)], alldata[, 1] >= sample.quality)
      submatrix <- list("quality" = sample.quality, "data" = submatrix)
    })




  } else {

    #List .csv sample files
    filesnames <- sapply(sample.table$Chip.ID.Well.ID, function(x) {

      csv.files <- list.files(
        file.location, pattern = paste0(as.character(x), "_Amplitude.csv"))

      if (length(csv.files) == 0)
        stop(paste0("Cannot find the file of sample ", x))

      if (length(csv.files) > 1)
        stop(paste0("Multiple files for sample", x))

      csv.files
    })

    subquality <- lapply(filesnames, function(y) {

      fluorescence <- utils::read.csv(paste0(file.location, "/", y),
                                      na.strings = c("NA", ""),
                                      colClasses = c("numeric", "numeric"))

      data <- matrix(c(fluorescence[,2],fluorescence[,1]), ncol = 2)
      colnames(data)  <-  c("Vic", "Fam")

      list("quality" = sample.quality, "data" = data)
    })
  }

  names(subquality) <- name

  class(subquality) <- "read_sample"

  return(subquality)
}
