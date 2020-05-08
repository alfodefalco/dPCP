#' Find the empty partitions and single target clusters in the reference sample
#'
#' This function computes a DBSCAN analysis to identify single target clusters
#' in the reference samples listed in the sample table.
#' If a \code{\link{reference_dbscan}} template file with the same input
#' paramters (reference ID, eps, minPts) is available, data are retrived
#' from the template file.
#' @param reference.subquality an object of class \code{read_reference},
#'   inherited from \code{\link{read_reference}}.
#' @inheritParams read_reference
#' @inheritParams dPCP
#' @rdname plot.reference_dbscan
#' @return An object of class \code{reference_dbscan} containing a sublist for
#'   each reference. Each sublist has the following components:
#'   \item{quality}{quality threshold used in \code{\link{read_reference}}.}
#'   \item{data}{a matrix with the fluorescence intensities and quality
#'   values.}
#'   \item{dbscan}{an object of class \code{dbscan_fast}, inherited from
#'   \code{\link[dbscan]{dbscan}}.}
#' @examples
#' \dontrun{
#' library(dPCP)
#'
#' #Find path of sample table and location of reference and input files
#' sampleTable <- system.file("extdata", "Template_sampleTable.csv",
#'                      package = "dPCP")
#'
#' fileLoc <- system.file("extdata",package = "dPCP")
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
#'
#' #Reference DBSCAN clustering
#' dbref <- reference_dbscan(ref, sample.table, save.template = FALSE)
#'
#' plot(dbref, reference = "all", save.plot = FALSE)
#' }
#' @export

reference_dbscan <- function(reference.subquality, sample.table, eps = 200,
                             minPts = 50, save.template = FALSE) {

  if (class(reference.subquality) != "read_reference")
    stop("reference.subquality must be an object of class read_reference")

  if (class(sample.table) != "sample_table")
    stop("'sample.table' must be an object of class sample_table")

  if (!is.numeric(c(eps, minPts))) stop("eps and minPts must be numeric")

  if (length(eps) > 1 & length(eps) != length(unique(sample.table$Reference)))
    stop("Invalid value for eps. It must be a numeric value or a numeric vector
         of lenght equal to the number of reference samples")

  if (length(minPts) > 1 & length(minPts) !=
      length(unique(sample.table$Reference)))
    stop("Invalid value for minPts. It must be a numeric value or a numeric
         vector of lenght equal to the number of reference samples")

  if (!is.logical(save.template) && !is.character(save.template))
    stop("save.template must be logical or for Thermo Fisher system a character
         vector with chip ID of reference samples to be saved")

  if (is.character(save.template) &&
      all(is.na(match(save.template, unique(sample.table$Reference)))))
    stop("save.template must be logical or for Thermo Fisher system a character
         vector with chip ID of reference samples to be saved")

  dbclus.ref <- lapply(seq_along(reference.subquality), function(x) {

    if (length(eps) == 1) {
      new.eps <- eps
    } else {
      new.eps <- eps[x]
    }

    if (length(minPts) == 1) {
      new.minPts <- minPts
    } else {
      new.minPts <- minPts[x]
    }

    #Check if data are dbscan results. If so, collect data; if not, calculate
    #dbscan results
    if (all(
      class(reference.subquality[[x]][[length(reference.subquality[[x]])]]) ==
      c("dbscan_fast", "dbscan"))) {

      dbclus <- reference.subquality[[x]]$dbscan

      ref.clus <- list(
        "quality" = reference.subquality[[x]]$quality,
        "data" = reference.subquality[[x]]$data, "dbscan" = dbclus)

      if (ref.clus$dbscan$eps != eps) {
        print(
          paste0(
            "The value of eps chosen for reference ",
            names(reference.subquality)[x],
            " does not match with eps in DBSCAN prediction file. Your choice
            will be discarded and DBSCAN prediction file data will be used.
            If you wish to use new settings, move DBSCAN prediction file to
            another folder"))
      }

      if (ref.clus$dbscan$minPts != minPts) {
        print(
          paste0(
            "The value of minPts chosen for reference ",
            names(reference.subquality)[x],
            " does not match with minPts in DBSCAN prediction file. Your choice
            will be discarded and DBSCAN prediction file data will be used.
            If you wish to use new settings, move DBSCAN prediction file to
            another folder"))
      }

    } else {
      dbclus <- dbscan::dbscan(reference.subquality[[x]]$data, new.eps,
                               new.minPts)

      ref.clus <- list("quality" = reference.subquality[[x]]$quality,
                       "data" = reference.subquality[[x]]$data,
                       "dbscan" = dbclus)
    }

    if (length(
      unique(ref.clus$dbscan$cluster)[unique(ref.clus$dbscan$cluster) != 0]) <
      sample.table$No.of.targets[
        match(names(reference.subquality)[x], sample.table$Reference)] + 1)
      stop(paste0("Number of clusters predicted for ", "'",
                  names(reference.subquality)[x], "'",
                  " are less than expected. Try to change eps and/or minPts
                  values."))

    if (any(!is.na(
      match(save.template, c(TRUE, names(reference.subquality)[x]))))) {

      not.allowed1 <- paste("\\\\", "/", ":", "\\*", "\\?", "\"", ">", "<",
                            "\\|", sep = "|")
      not.allowed2 <- c("CON", "PRN", "AUX", "NUL", "COM1", "COM2", "COM3",
                        "COM4", "COM5", "COM6", "COM7", "COM8", "COM9", "LPT1",
                        "LPT2", "LPT3", "LPT4", "LPT5", "LPT6", "LPT7", "LPT8",
                        "LPT9")

      if (grepl(not.allowed1, names(reference.subquality)[x])) {

        new.name <- (gsub(not.allowed1, ".", names(reference.subquality)[x]))

      } else if (any(names(reference.subquality)[x] == not.allowed2)) {

        new.name <- paste0("Reference", x)

      } else {
        new.name <- names(reference.subquality)[x]
      }

      if (is.character(reference.subquality[[x]]$quality) &
          length(reference.subquality[[x]]) == 2) {
        saveRDS(ref.clus, file =  paste0(
          stringr::str_sub(names(reference.subquality)[x], end = -5),
          "_reference_", "eps", new.eps, "_", "minPts", new.minPts, "_DB",
          ".rds"))

      } else if (is.numeric(reference.subquality[[x]]$quality) &
                 (length(reference.subquality[[x]]) == 2)) {
        saveRDS(ref.clus, file =  paste0(
          new.name, "_", "qlty", reference.subquality[[x]]$quality, "_", "eps",
          new.eps, "_", "minPts", new.minPts, "_DB", ".rds"))
      }
    }

    ref.clus
  })
  names(dbclus.ref) <- names(reference.subquality)

  class(dbclus.ref) <- "reference_dbscan"

  return(dbclus.ref)
}


#' Plot results of DBSCAN analysis for reference samples
#'
#' @param x an object of class \code{reference_dbscan}
#' @param ... Arguments to be passed to methods
#' @inheritParams dPCP
#' @export
#' @import ggplot2

plot.reference_dbscan <- function(x, ..., reference = "all",
                                  save.plot = FALSE, format = "png",
                                  dpi = 300) {

  if (any(reference != "all") & any(is.na(match(reference, names(x)))))
    stop("'reference' must be `all` or a character vector with chip ID (Thermo
         Fisher) or the file name (Bio-rad) of reference samples to be showed")

  if (!is.logical(save.plot))
    stop("'save.plot' must be logical")

  if (all(format != c("eps", "ps", "tex", "pdf", "jpeg", "tiff", "png",
                      "bmp", "svg", "wmf")))
    stop("Uknown format")

  if (!is.numeric(dpi) || length(dpi) != 1)
    stop("'dpi' must be a numeric value")

  if (all(reference == "all")) {
    plotreference <- 1:length(x)
  } else {
    plotreference <- sort(match(reference, names(x)))
  }

  refDBplot <- lapply(plotreference, function(y) {

    p <- ggplot(as.data.frame(x[[y]]$data),
                aes_string(x = "Vic", y = "Fam")) +

      geom_point(
        size = 0.9, aes(color = factor(x[[y]]$dbscan$cluster))) +

      geom_point(
        size = 0.9,
        data = as.data.frame(
          x[[y]]$data)[x[[y]]$dbscan$cluster == 0, ],
        aes_string(x = "Vic", y = "Fam"), colour = "grey40") +

      labs(title = paste0("DBSCAN reference '", names(x)[y], "'",
                          "\neps = ", x[[y]]$dbscan$eps,
                          ", minPts = ", x[[y]]$dbscan$minPts)) +

      theme(
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none",
        plot.title = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 9)
      )

    if (isTRUE(save.plot)) {

      not.allowed1 <- paste("\\\\", "/", ":", "\\*", "\\?", "\"", ">", "<",
                            "\\|", sep = "|")

      if (grepl(not.allowed1, names(x)[y])) {
        new.name <- (gsub(not.allowed1, ".", names(x)[y]))
      } else {
        new.name <- names(x)[y]
      }

      ggsave(paste0(new.name, "_refDBplot", ".", format), p, dpi = dpi)
    }

    graphics::plot(p)
  })
  names(refDBplot) <- names(x)[plotreference]

  return(refDBplot)
}


#' Predict empty partitions and single target clusters in the samples (deep
#' workflow)
#'
#' This function predicts the single target clusters in the samples using the
#' DBSCAN analysis of a reference as guide.
#' @param sample.subquality an object of class \code{read_sample}, inherited
#'   from \code{\link{read_sample}}.
#' @inheritParams read_reference
#' @param  referenceDB an object of class \code{reference_dbscan}, inherited
#'   from \code{\link{reference_dbscan}}
#' @rdname plot.sample_dbscan
#' @return An object of class \code{sample_dbscan} containing a sublist for
#'   each sample. Each sublist has the following components:
#'   \item{quality}{quality threshold used in \code{\link{read_sample}}.}
#'   \item{reference}{reference ID.}
#'   \item{dbscan}{an object of class \code{dbscan_fast}, inherited from
#'     \code{\link[dbscan]{dbscan}}.}
#' @examples
#' \dontrun{
#' library(dPCP)
#'
#' #Find path of sample table and location of reference and input files
#' sampleTable <- system.file("extdata", "Template_sampleTable.csv",
#'                      package = "dPCP")
#'
#' fileLoc <- system.file("extdata",package = "dPCP")
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
#'
#' #Reference DBSCAN clustering
#' dbref <- reference_dbscan(ref, sample.table, save.template = FALSE)
#'
#' #Sample DBSCAN prediction
#' dbsam <- sample_dbscan(samp, sample.table, dbref, file.location = fileLoc)
#'
#' plot(dbsam, sample = "all", save.plot = FALSE)
#' }
#' @export

sample_dbscan <- function(sample.subquality, sample.table, referenceDB,
                          file.location = ".") {

  if (class(sample.subquality) != "read_sample")
    stop("sample.subquality must be an object of class read_sample")

  if (class(sample.table) != "sample_table")
    stop("'sample.table' must be an object of class sample_table")

  if (class(referenceDB) != "reference_dbscan")
    stop("referenceDB must be an object of class reference_dbscan")

  DB <- lapply(seq_along(sample.subquality), function(x) {

    #Select the appropriate reference
    reference <- match(sample.table$Reference[x], names(referenceDB))

    if (is.character(sample.subquality[[x]]$quality)) {
      #Check if sample ID = reference ID.
      #If so,take dbscan results from previous step; if not, do the prediction
      if (list.files(
        file.location,
        pattern = paste0(as.character(sample.table$Chip.ID.Well.ID[x]),
                         "_Amplitude.csv")) == sample.table$Reference[x]) {

        datacluster <- cbind.data.frame(
          sample.subquality[[x]]$data,
          "cluster" = referenceDB[[reference]]$dbscan$cluster)
      } else {
        dbclus <- stats::predict(referenceDB[[reference]]$dbscan,
                                 data = referenceDB[[reference]]$data,
                                 newdata = sample.subquality[[x]]$data)

        datacluster <- cbind.data.frame(sample.subquality[[x]]$data,
                                        "cluster" = dbclus)
      }
    } else {
      #Check if sample ID = reference ID.
      #If so,take dbscan results from previous step; if not, do the prediction
      if (sample.table$Chip.ID[x] == sample.table$Reference[x] &
          sample.subquality[[x]]$quality == referenceDB[[reference]]$quality) {

        datacluster <- cbind.data.frame(
          sample.subquality[[x]]$data,
          "cluster" = referenceDB[[reference]]$dbscan$cluster)
      } else {
        dbclus <- stats::predict(referenceDB[[reference]]$dbscan,
                                 data = referenceDB[[reference]]$data,
                                 newdata = sample.subquality[[x]]$data)

        datacluster <- cbind.data.frame(sample.subquality[[x]]$data,
                                        "cluster" = dbclus)
      }
    }

    list("quality" = sample.subquality[[x]]$quality,
         "reference" = sample.table$Reference[x],
         "data" = datacluster)
  })
  names(DB) <- names(sample.subquality)

  class(DB) <- "sample_dbscan"

  return(DB)
}


#' Plot results of DBSCAN analysis of samples
#'
#' @param x an object of class \code{sample_dbscan}
#' @param ... Arguments to be passed to methods
#' @inheritParams dPCP
#' @export

plot.sample_dbscan <- function(x, ..., sample = "all", save.plot = FALSE,
                               format = "png", dpi = 300) {

  if (all(sample != "all") & any(sample > length(x)))
    stop("sample must be `all` or a numeric vector indicating the row number of
         samples in sample.table")

  if (!is.logical(save.plot))
    stop("save.plot must be logical")

  if (all(format != c("eps", "ps", "tex", "pdf", "jpeg", "tiff", "png", "bmp",
                      "svg", "wmf")))
    stop("Uknown format")

  if (!is.numeric(dpi) || length(dpi) != 1) stop("dpi must be a numeric value")


  if (all(sample == "all")) {
    plotsample <- 1:length(x)
  } else {
    plotsample <- sort(sample)
  }

  sampleDBplot <- lapply(plotsample, function(y) {

    p <- ggplot(as.data.frame(x[[y]]$data),
                aes_string(x = "Vic", y = "Fam")) +

      geom_point(size = 1, aes(color = factor(x[[y]]$data[, 3]))) +

      geom_point(
        size = 1,
        data = as.data.frame(
          x[[y]]$data)[x[[y]]$data[, 3] == 0, ],
        aes_string(x = "Vic", y = "Fam"), colour = "grey40") +

      labs(title = names(x)[y]) +

      theme(
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none",
        plot.title = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 9)
      )

    if (isTRUE(save.plot)) {

      not.allowed1 <- paste("\\\\", "/", ":", "\\*", "\\?", "\"", ">", "<",
                            "\\|", sep = "|")

      if (grepl(not.allowed1, names(x)[y])) {
        new.name <- (gsub(not.allowed1, ".", names(x)[y]))
      } else {
        new.name <- names(x)[y]
      }

      ggsave(paste0(new.name, "_sampleDBplot", ".", format), p, dpi = dpi)
    }

    graphics::plot(p)
  })
  names(sampleDBplot) <- names(x)[plotsample]

  return(sampleDBplot)
}


#' Prediction of clusters centroid position
#'
#' This function calculates the centroids position of empty partitions and
#' single target clusters, then predicts the centroids position of multi target
#' clusters.
#' While \code{centers_data} (default workflow) calculates the centroids
#' position of empty partitions and single target clusters from the DBSCAN
#' analysis of reference, \code{centers_deep} (deep workflow) takes as input
#' the DBSCAN analysis of samples.
#' @inheritParams sample_dbscan
#' @rdname plot.centers_data
#' @return An object of class \code{centers_data} containing a sublist for
#'   each sample. Each sublist has the following components:
#'   \item{quality}{quality threshold used in \code{\link{read_sample}}.}
#'   \item{reference}{reference ID.}
#'   \item{centers}{a data frame with the centroids coordinates.}
#'   \item{data}{a data frame with the fluorescence intensities.}
#' @examples
#' \dontrun{
#' library(dPCP)
#'
#' #Find path of sample table and location of reference and input files
#' sampleTable <- system.file("extdata", "Template_sampleTable.csv",
#'                      package = "dPCP")
#'
#' fileLoc <- system.file("extdata",package = "dPCP")
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
#'
#' #Reference DBSCAN clustering
#' dbref <- reference_dbscan(ref, sample.table, save.template = FALSE)
#'
#' #Predict position of clusters centroid from reference DBSCAN results
#' cent <- centers_data(samp, sample.table,dbref)
#'
#' plot(cent, sample = "all", save.plot = FALSE)
#' }
#' @export

centers_data <- function(sample.subquality, sample.table,
                         referenceDB) {

  if (class(sample.subquality) != "read_sample")
    stop("sample.subquality must be an object of class read_sample")

  if (class(sample.table) != "sample_table")
    stop("'sample.table' must be an object of class sample_table")

  class(sample.table) <- "data.frame"

  if (class(referenceDB) != "reference_dbscan")
    stop("referenceDB must be an object of class reference_dbscan")

  centersCoor <- lapply(seq_along(sample.subquality), function(x) {

    mrks <- sample.table[x, 3]

    reference <- which(names(referenceDB) == sample.table[x, 8])

    notnoise.refdata <- subset(referenceDB[[reference]]$data,
                               referenceDB[[reference]]$dbscan$cluster > 0)

    notnoise.refcluster <- subset(referenceDB[[reference]]$dbscan$cluster,
                                  referenceDB[[reference]]$dbscan$cluster > 0)

    notnoiseDB <- cbind.data.frame(
      notnoise.refdata, "cluster" = notnoise.refcluster)

    #Calculate position of reference sample centers
    cent.ref <- lapply(unique(notnoiseDB$cluster), function(y) {
      vic.centers <- mean(subset(notnoiseDB[, 1], notnoiseDB$cluster == y))
      fam.centers <- mean(subset(notnoiseDB[, 2], notnoiseDB$cluster == y))

      cbind.data.frame("Vic" = vic.centers, "Fam" = fam.centers)
    })
    centersDB <- do.call(rbind.data.frame, cent.ref)

    centersDB$dist <- raster::pointDistance(centersDB, c(0, 0), lonlat = FALSE)

    empty <- subset(centersDB, centersDB$dist == min(centersDB$dist))
    row.names(empty) <- 1
    noEmpty <- subset(centersDB, centersDB$dist != min(centersDB$dist))
    row.names(noEmpty) <- 1:nrow(noEmpty)

    if (mrks == 1) {
      allcenters <- centersDB[, c(1, 2)][order(centersDB$dist), ]

      if (length(which(!is.na(sample.table[x, 4:7]))) == 1) {
        row.names(allcenters) <- c(
          "Empty",
          sample.table[x, which(!is.na(sample.table[x, 4:7])) + 3])
      } else {
        row.names(allcenters) <- c("Empty", "Target1")
      }

      cent.data <- list("quality" = sample.subquality[[x]]$quality,
                        "reference" = names(referenceDB)[reference],
                        "centers" = allcenters,
                        "data" = data.frame(sample.subquality[[x]]$data))
    } else if (mrks > 1) {
      targetFam <- which.min(noEmpty$Vic)
      targetVic <- which.min(noEmpty$Fam)

      if (mrks == 2) {
        centMat <- rbind(empty[, c(1, 2)], noEmpty[targetFam, c(1, 2)],
                         noEmpty[targetVic, c(1, 2)])

        if (length(which(!is.na(sample.table[x, 4:7]))) == 2) {
          row.names(centMat) <- c(
            "Empty",
            sample.table[x, which(!is.na(sample.table[x, 4:7])) + 3])
        } else {
          row.names(centMat) <- c("Empty", "Target1", "Target2")
        }

        doublePos <- centMat[2, ] + centMat[3, ] - centMat[1, ]

        allcenters <- rbind(centMat, doublePos)
        row.names(allcenters) <- c(
          row.names(centMat),
          paste(row.names(centMat)[2], row.names(centMat)[3], sep = " + "))

        cent.data <- list("quality" = sample.subquality[[x]]$quality,
                          "reference" = names(referenceDB)[reference],
                          "centers" = allcenters,
                          "data" = data.frame(sample.subquality[[x]]$data))

      } else if (mrks == 3) {
        notVicFam <- noEmpty[-c(targetFam, targetVic), ]

        targetmix <- subset(notVicFam, notVicFam$dist == min(notVicFam$dist))

        centMat <- rbind(empty[, c(1, 2)],
                         noEmpty[targetFam, c(1, 2)],
                         targetmix[, c(1, 2)],
                         noEmpty[targetVic, c(1, 2)])

        if (length(which(!is.na(sample.table[x, 4:7]))) == 3) {
          row.names(centMat) <- c(
            "Empty",
            sample.table[x, which(!is.na(sample.table[x, 4:7])) + 3])
        } else {
          row.names(centMat) <- c("Empty", "Fam", "Mix1", "Vic")
        }

        FamMix1 <- centMat[2, ] + centMat[3, ] - centMat[1, ]

        FamVic <- centMat[2, ] + centMat[4, ] - centMat[1, ]

        Mix1Vic <- centMat[3, ] + centMat[4, ] - centMat[1, ]

        FamMix1Vic <- FamMix1 + centMat[4, ] - centMat[1, ]

        allcenters <- rbind(centMat, FamMix1, FamVic, Mix1Vic, FamMix1Vic)

        row.names(allcenters) <- c(row.names(centMat),
                                   paste(row.names(centMat)[2],
                                         row.names(centMat)[3], sep = " + "),
                                   paste(row.names(centMat)[2],
                                         row.names(centMat)[4], sep = " + "),
                                   paste(row.names(centMat)[3],
                                         row.names(centMat)[4], sep = " + "),
                                   paste(row.names(centMat)[2],
                                         row.names(centMat)[3],
                                         row.names(centMat)[4], sep = " + "))

        cent.data <- list("quality" = sample.subquality[[x]]$quality,
                          "reference" = names(referenceDB)[reference],
                          "centers" = allcenters,
                          "data" = data.frame(sample.subquality[[x]]$data))

      } else if (mrks == 4) {

        notVicFam <- noEmpty[-c(targetFam, targetVic), ]
        notVicFam <- notVicFam[, c(1, 2)][order(notVicFam$dist), ]
        targetsmix <- notVicFam[c(1, 2), ]
        targetsmix <- targetsmix[order(targetsmix$Vic), ]

        centMat <- rbind(empty[, c(1, 2)],
                         noEmpty[targetFam, c(1, 2)],
                         targetsmix[1, c(1, 2)],
                         targetsmix[2, c(1, 2)],
                         noEmpty[targetVic, c(1, 2)])

        if (length(which(!is.na(sample.table[x, 4:7]))) == 4) {
          row.names(centMat) <- c(
            "Empty",
            sample.table[x, which(!is.na(sample.table[x, 4:7])) + 3])
        } else {
          row.names(centMat) <- c("Empty", "Fam", "Target3",
                                  "Target4", "Vic")
        }


        FamTarget3 <- centMat[2, ] + centMat[3, ] - centMat[1, ]

        FamTarget4 <- centMat[2, ] + centMat[4, ] - centMat[1, ]

        Target3Target4 <- centMat[3, ] + centMat[4, ] - centMat[1, ]

        FamVic <- centMat[2, ] + centMat[5, ] - centMat[1, ]

        Target3Vic <- centMat[3, ] + centMat[5, ] - centMat[1, ]

        Target4Vic <- centMat[4, ] + centMat[5, ] - centMat[1, ]

        FamTarget3Target4 <- FamTarget3 + centMat[4, ] - centMat[1, ]

        FamTarget3Vic <- FamTarget3 + centMat[5, ] - centMat[1, ]

        FamTarget4Vic <- FamTarget4 + centMat[5, ] - centMat[1, ]

        Target3Target4Vic <- Target3Target4 + centMat[5, ] - centMat[1, ]

        FamTarget3Target4Vic <- FamTarget3Target4 + centMat[5, ] - centMat[1, ]


        allcenters <- rbind(centMat, FamTarget3, FamTarget4, Target3Target4,
                            FamVic, Target3Vic, Target4Vic, FamTarget3Target4,
                            FamTarget3Vic, FamTarget4Vic, Target3Target4Vic,
                            FamTarget3Target4Vic)

        row.names(allcenters) <- c(
          row.names(centMat),
          paste(row.names(centMat)[2], row.names(centMat)[3], sep = " + "),
          paste(row.names(centMat)[2], row.names(centMat)[4], sep = " + "),
          paste(row.names(centMat)[3], row.names(centMat)[4], sep = " + "),
          paste(row.names(centMat)[2], row.names(centMat)[5], sep = " + "),
          paste(row.names(centMat)[3], row.names(centMat)[5], sep = " + "),
          paste(row.names(centMat)[4], row.names(centMat)[5], sep = " + "),
          paste(row.names(centMat)[2], row.names(centMat)[3],
                row.names(centMat)[4], sep = " + "),
          paste(row.names(centMat)[2], row.names(centMat)[3],
                row.names(centMat)[5], sep = " + "),
          paste(row.names(centMat)[2], row.names(centMat)[4],
                row.names(centMat)[5], sep = " + "),
          paste(row.names(centMat)[3], row.names(centMat)[4],
                row.names(centMat)[5], sep = " + "),
          paste(row.names(centMat)[2], row.names(centMat)[3],
                row.names(centMat)[4], row.names(centMat)[5], sep = " + "))

        cent.data <- list("quality" = sample.subquality[[x]]$quality,
                          "reference" = names(referenceDB)[reference],
                          "centers" = allcenters,
                          "data" = data.frame(sample.subquality[[x]]$data))
      }
    }

    cent.data
  })
  names(centersCoor) <- names(sample.subquality)

  class(centersCoor) <- "centers_data"

  return(centersCoor)
}


#' Prediction of clusters centroid position
#'
#' @param sample.dbscan an object of class \code{sample_dbscan}, inherited
#'   from \code{\link{sample_dbscan}}.
#' @inheritParams sample_dbscan
#' @rdname plot.centers_data
#' @examples
#' \dontrun{
#' library(dPCP)
#'
#' #Find path of sample table and location of reference and input files
#' sampleTable <- system.file("extdata", "Template_sampleTable.csv",
#'                      package = "dPCP")
#'
#' fileLoc <- system.file("extdata",package = "dPCP")
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
#'
#' #Reference DBSCAN clustering
#' dbref <- reference_dbscan(ref, sample.table, save.template = FALSE)
#'
#' #Sample DBSCAN prediction
#' dbsam <- sample_dbscan(samp, sample.table, dbref, file.location = fileLoc)
#'
#' #Predict clusters centroid position from DBSCAN prediction
#' cent <- centers_deep(dbsam, sample.table,dbref)
#'
#' plot(cent, sample = "all", save.plot = FALSE)
#' }
#' @export

centers_deep <- function(sample.dbscan, sample.table, referenceDB) {

  if (class(sample.dbscan) != "sample_dbscan")
    stop("sample.dbscan must be an object of class sample_dbscan")

  if (class(sample.table) != "sample_table")
    stop("'sample.table' must be an object of class sample_table")

  class(sample.table) <- "data.frame"

  if (class(referenceDB) != "reference_dbscan")
    stop("referenceDB must be an object of class reference_dbscan")

  centersCoor <- lapply(seq_along(sample.dbscan), function(x) {

    mrks <- sample.table[x, 3]

    notnoise <- subset(sample.dbscan[[x]]$data,
                       sample.dbscan[[x]]$data$cluster > 0)

    reference <- which(names(referenceDB) == sample.table[x, 8])

    notnoise.refdata <- subset(referenceDB[[reference]]$data,
                               referenceDB[[reference]]$dbscan$cluster > 0)

    notnoise.refcluster <- subset(referenceDB[[reference]]$dbscan$cluster,
                                  referenceDB[[reference]]$dbscan$cluster > 0)

    notnoiseDB <- cbind.data.frame(
      notnoise.refdata, "cluster" = notnoise.refcluster)

    #if DBSCAN clusters identified in the samples have less than 10 elements,
    #the center of that cluster is predicted using reference data.
    if (any(table(notnoise$cluster) < 10) ||
        any(!unique(notnoiseDB$cluster) %in% unique(notnoise$cluster))) {

      missed1 <- as.numeric(
        subset(names(table(notnoise$cluster)), table(notnoise$cluster) < 10))

      missed2 <- unique(notnoiseDB$cluster)[!unique(notnoiseDB$cluster) %in%
                                              unique(notnoise$cluster)]

      allmissed <- c(missed1, missed2)

      notmissed <- unique(
        notnoiseDB$cluster)[!unique(notnoiseDB$cluster) %in% allmissed]

      missedcent.ref <- lapply(allmissed, function(y) {
        vic.centers <- mean(subset(notnoiseDB[, 1], notnoiseDB$cluster == y))
        fam.centers <- mean(subset(notnoiseDB[, 2], notnoiseDB$cluster == y))

        cbind.data.frame("Vic" = vic.centers, "Fam" = fam.centers)
      })
      missedcent.ref <- do.call(rbind.data.frame, missedcent.ref)

      notmissedcent.ref <- lapply(notmissed, function(z) {
        vic.centers <- mean(subset(
          referenceDB[[reference]]$data[, 1],
          referenceDB[[reference]]$dbscan$cluster == z))

        fam.centers <- mean(subset(
          referenceDB[[reference]]$data[, 2],
          referenceDB[[reference]]$dbscan$cluster == z))

        centersDB.ref <- cbind.data.frame("Vic" = vic.centers,
                                          "Fam" = fam.centers)
      })

      notmissedcent.ref <- do.call(rbind.data.frame, notmissedcent.ref)

      cent.sam <- lapply(notmissed, function(k) {
        vic.centers <- mean(subset(
          sample.dbscan[[x]]$data[, 1], sample.dbscan[[x]]$data$cluster == k))

        fam.centers <- mean(subset(
          sample.dbscan[[x]]$data[, 2], sample.dbscan[[x]]$data$cluster == k))

        centersDB.sam <- cbind.data.frame("Vic" = vic.centers,
                                          "Fam" = fam.centers)
      })

      cent.sam <- do.call(rbind.data.frame, cent.sam)
      dif <- colMeans(cent.sam - notmissedcent.ref)

      cent.missed.vic <- missedcent.ref[, 1] + dif[1]
      cent.missed.fam <- missedcent.ref[, 2] + dif[2]
      cent.missed <- cbind.data.frame("Vic" = cent.missed.vic,
                                      "Fam" = cent.missed.fam)

      centersDB <- rbind.data.frame(cent.missed, cent.sam)
      centersDB$dist <- raster::pointDistance(centersDB, c(0, 0),
                                              lonlat = FALSE)
    }
    else {
      cent.dist <- lapply(sort(unique(notnoise$cluster)), function(j) {
        vic.centers <- mean(subset(
          sample.dbscan[[x]]$data[, 1], sample.dbscan[[x]]$data$cluster == j))

        fam.centers <- mean(subset(
          sample.dbscan[[x]]$data[, 2], sample.dbscan[[x]]$data$cluster == j))

        centersDB <- cbind.data.frame("Vic" = vic.centers, "Fam" = fam.centers)
        centersDB$dist <- raster::pointDistance(centersDB, c(0, 0),
                                                lonlat = FALSE)
        centersDB
      })
      centersDB <- do.call(rbind.data.frame, cent.dist)
    }

    empty <- subset(centersDB, centersDB$dist == min(centersDB$dist))
    row.names(empty) <- 1
    noEmpty <- subset(centersDB, centersDB$dist != min(centersDB$dist))
    row.names(noEmpty) <- 1:nrow(noEmpty)

    if (mrks == 1) {
      allcenters <- centersDB[, c(1, 2)][order(centersDB$dist), ]

      if (length(which(!is.na(sample.table[x, 4:7]))) == 1) {
        row.names(allcenters) <- c(
          "Empty",
          sample.table[x, which(!is.na(sample.table[x, 4:7])) + 3])
      } else {
        row.names(allcenters) <- c("Empty", "Target1")
      }

      cent.data <- list("quality" = sample.dbscan[[x]]$quality,
                        "reference" = sample.dbscan[[x]]$reference,
                        "centers" = allcenters,
                        "data" = sample.dbscan[[x]]$data[, c(1, 2)])
    } else if (mrks > 1) {
      targetFam <- which.min(noEmpty$Vic)
      targetVic <- which.min(noEmpty$Fam)

      if (mrks == 2) {
        centMat <- rbind(empty[, c(1, 2)], noEmpty[targetFam, c(1, 2)],
                         noEmpty[targetVic, c(1, 2)])

        if (length(which(!is.na(sample.table[x, 4:7]))) == 2) {
          row.names(centMat) <- c(
            "Empty",
            sample.table[x, which(!is.na(sample.table[x, 4:7])) + 3])
        } else {
          row.names(centMat) <- c("Empty", "Target1", "Target2")
        }

        doublePos <- centMat[2, ] + centMat[3, ] - centMat[1, ]

        allcenters <- rbind(centMat, doublePos)
        row.names(allcenters) <- c(
          row.names(centMat),
          paste(row.names(centMat)[2], row.names(centMat)[3], sep = " + "))

        cent.data <- list("quality" = sample.dbscan[[x]]$quality,
                          "reference" = sample.dbscan[[x]]$reference,
                          "centers" = allcenters,
                          "data" = sample.dbscan[[x]]$data[, c(1, 2)])

      } else if (mrks == 3) {
        notVicFam <- noEmpty[-c(targetFam, targetVic), ]

        targetmix <- subset(notVicFam, notVicFam$dist == min(notVicFam$dist))

        centMat <- rbind(empty[, c(1, 2)],
                         noEmpty[targetFam, c(1, 2)],
                         targetmix[, c(1, 2)],
                         noEmpty[targetVic, c(1, 2)])

        if (length(which(!is.na(sample.table[x, 4:7]))) == 3) {
          row.names(centMat) <- c(
            "Empty",
            sample.table[x, which(!is.na(sample.table[x, 4:7])) + 3])
        } else {
          row.names(centMat) <- c("Empty", "Fam", "Mix1", "Vic")
        }

        FamMix1 <- centMat[2, ] + centMat[3, ] - centMat[1, ]

        FamVic <- centMat[2, ] + centMat[4, ] - centMat[1, ]

        Mix1Vic <- centMat[3, ] + centMat[4, ] - centMat[1, ]

        FamMix1Vic <- FamMix1 + centMat[4, ] - centMat[1, ]

        allcenters <- rbind(centMat, FamMix1, FamVic, Mix1Vic, FamMix1Vic)

        row.names(allcenters) <- c(row.names(centMat),
                                   paste(row.names(centMat)[2],
                                         row.names(centMat)[3], sep = " + "),
                                   paste(row.names(centMat)[2],
                                         row.names(centMat)[4], sep = " + "),
                                   paste(row.names(centMat)[3],
                                         row.names(centMat)[4], sep = " + "),
                                   paste(row.names(centMat)[2],
                                         row.names(centMat)[3],
                                         row.names(centMat)[4], sep = " + "))

        cent.data <- list("quality" = sample.dbscan[[x]]$quality,
                          "reference" = sample.dbscan[[x]]$reference,
                          "centers" = allcenters,
                          "data" = sample.dbscan[[x]]$data[, c(1, 2)])

      } else if (mrks == 4) {

        notVicFam <- noEmpty[-c(targetFam, targetVic), ]
        notVicFam <- notVicFam[, c(1, 2)][order(notVicFam$dist), ]
        targetsmix <- notVicFam[c(1, 2), ]
        targetsmix <- targetsmix[order(targetsmix$Vic), ]

        centMat <- rbind(empty[, c(1, 2)],
                         noEmpty[targetFam, c(1, 2)],
                         targetsmix[1, c(1, 2)],
                         targetsmix[2, c(1, 2)],
                         noEmpty[targetVic, c(1, 2)])

        if (length(which(!is.na(sample.table[x, 4:7]))) == 4) {
          row.names(centMat) <- c(
            "Empty",
            sample.table[x, which(!is.na(sample.table[x, 4:7])) + 3])
        } else {
          row.names(centMat) <- c("Empty", "Fam", "Target3",
                                  "Target4", "Vic")
        }


        FamTarget3 <- centMat[2, ] + centMat[3, ] - centMat[1, ]

        FamTarget4 <- centMat[2, ] + centMat[4, ] - centMat[1, ]

        Target3Target4 <- centMat[3, ] + centMat[4, ] - centMat[1, ]

        FamVic <- centMat[2, ] + centMat[5, ] - centMat[1, ]

        Target3Vic <- centMat[3, ] + centMat[5, ] - centMat[1, ]

        Target4Vic <- centMat[4, ] + centMat[5, ] - centMat[1, ]

        FamTarget3Target4 <- FamTarget3 + centMat[4, ] - centMat[1, ]

        FamTarget3Vic <- FamTarget3 + centMat[5, ] - centMat[1, ]

        FamTarget4Vic <- FamTarget4 + centMat[5, ] - centMat[1, ]

        Target3Target4Vic <- Target3Target4 + centMat[5, ] - centMat[1, ]

        FamTarget3Target4Vic <- FamTarget3Target4 + centMat[5, ] - centMat[1, ]


        allcenters <- rbind(centMat, FamTarget3, FamTarget4, Target3Target4,
                            FamVic, Target3Vic, Target4Vic, FamTarget3Target4,
                            FamTarget3Vic, FamTarget4Vic, Target3Target4Vic,
                            FamTarget3Target4Vic)

        row.names(allcenters) <- c(
          row.names(centMat),
          paste(row.names(centMat)[2], row.names(centMat)[3], sep = " + "),
          paste(row.names(centMat)[2], row.names(centMat)[4], sep = " + "),
          paste(row.names(centMat)[3], row.names(centMat)[4], sep = " + "),
          paste(row.names(centMat)[2], row.names(centMat)[5], sep = " + "),
          paste(row.names(centMat)[3], row.names(centMat)[5], sep = " + "),
          paste(row.names(centMat)[4], row.names(centMat)[5], sep = " + "),
          paste(row.names(centMat)[2], row.names(centMat)[3],
                row.names(centMat)[4], sep = " + "),
          paste(row.names(centMat)[2], row.names(centMat)[3],
                row.names(centMat)[5], sep = " + "),
          paste(row.names(centMat)[2], row.names(centMat)[4],
                row.names(centMat)[5], sep = " + "),
          paste(row.names(centMat)[3], row.names(centMat)[4],
                row.names(centMat)[5], sep = " + "),
          paste(row.names(centMat)[2], row.names(centMat)[3],
                row.names(centMat)[4], row.names(centMat)[5], sep = " + "))

        cent.data <- list("quality" = sample.dbscan[[x]]$quality,
                          "reference" = sample.dbscan[[x]]$reference,
                          "centers" = allcenters,
                          "data" = sample.dbscan[[x]]$data[, c(1, 2)])
      }
    }

    cent.data
  })
  names(centersCoor) <- names(sample.dbscan)

  class(centersCoor) <- "centers_data"

  return(centersCoor)
}


#' Plot clusters centroid
#'
#' @param x an object of class \code{centers_data}
#' @param ... Arguments to be passed to methods
#' @inheritParams dPCP
#' @examples
#' \dontrun{
#' plot(cent, sample = c(2,5), save.plot = TRUE, format = "jpg",
#'       dpi = 300)
#' }
#' @export

plot.centers_data <- function(x, ..., sample = "all", save.plot = FALSE,
                              format = "png", dpi = 300) {

  if (all(sample != "all") & any(sample > length(x)))
    stop("sample must be `all` or a numeric vector indicating the row number of
         samples in sample.table")

  if (!is.logical(save.plot))
    stop("save.plot must be logical")

  if (all(format != c("eps", "ps", "tex", "pdf", "jpeg", "tiff", "png", "bmp",
                      "svg", "wmf")))
    stop("Uknown format")

  if (!is.numeric(dpi) || length(dpi) != 1) stop("dpi must be a numeric value")


  if (all(sample == "all")) {
    plotsample <- 1:length(x)
  } else {
    plotsample <- sort(sample)
  }

  centersplot <- lapply(plotsample, function(y) {

    p <- ggplot(as.data.frame(x[[y]]$data),
                aes_string(x = "Vic", y = "Fam")) +

      geom_point(size = 1) +

      geom_point(data = as.data.frame(x[[y]]$centers), size = 2,
                 aes(x = x[[y]]$centers[, 1],
                     y = x[[y]]$centers[, 2]), colour = "red") +

      labs(title = names(x)[y]) +

      theme(
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black"),
        plot.title = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 9)
      )

    if (isTRUE(save.plot)) {

      not.allowed1 <- paste("\\\\", "/", ":", "\\*", "\\?", "\"", ">", "<",
                            "\\|", sep = "|")

      if (grepl(not.allowed1, names(x)[y])) {
        new.name <- (gsub(not.allowed1, ".", names(x)[y]))
      } else {
        new.name <- names(x)[y]
      }

      ggsave(paste0(new.name, "_centersplot", ".", format), p, dpi = dpi)
    }

    graphics::plot(p)
  })
  names(centersplot) <- names(x)[plotsample]

  return(centersplot)
}

