#' Automated analysis of digital PCR data
#'
#' This function carries out the autometed clustering of digital PCR data.
#' @param  file character. The name or the path of csv file to be read. If it
#' does not contain an absolute path, the file name is relative to the current
#' working directory, (\code{\link[base]{getwd}}).
#' @param  system character. The name of digital PCR system used to generate
#'   the data. It must be either Thermo Fisher or Bio-Rad. Abbreviations are
#'   also accepted.
#' @param file.location character. Full path name to reference and sample
#' files location. The default corresponds to the working directory,
#' (\code{\link[base]{getwd}}). Tilde expansion (see
#' (\code{\link[base]{path.expand}})) is performed.
#' @param  reference.quality numeric. Between 0 and 1. Quality threshold to
#'   subset the data. If different thresholds have to be applied to various
#'   reference samples, a vectror of the same length of number of reference
#'   samples has to be provided. Used only when the \code{system} is Thermo
#'   Fisher.
#' @param  sample.quality numeric. Between 0 and 1. Quality threshold to subset
#'   data. If different thresholds have to be applied to various samples, a
#'   vectror of the same length of number of samples has to be provided. Used
#'   only when the \code{system} is Thermo Fisher.
#' @param  eps numeric. Input parameter for the DBSCAN algorithm.
#'   It represents the maximum distance between the elements within a cluster.
#'   See also \code{\link[dbscan]{dbscan}}. If different values have to be
#'   applied to various reference samples, a vectror of the same length of
#'   number of reference samples has to be provided.
#' @param  minPts numeric. Input parameter for the DBSCAN algorithm.
#'   It represents the number of minimum elements to assemble a cluster. See
#'   also \code{\link[dbscan]{dbscan}}. If different values have to be
#'   applied to various reference samples, a vectror of the same length of
#'   number of reference samples has to be provided.
#' @param  save.template logical. If TRUE a template of DBSCAN analysis of
#'   reference samples is saved. When \code{system} is Thermo Fisher,
#'   \code{save.template} can be also a character vector indicating the chipID.
#' @param rain logical. If TRUE the rain analysis is carried out.
#' @rdname plot.dPCP
#' @return An object of class \code{dPCP} containing the following components:
#'   \item{referenceDB}{an object of class \code{reference_dbscan}.}
#'   \item{samples}{a list of samples. Each sample sublist contains the
#'   information about the cluster analysis.}
#'   \item{results}{an object of class \code{replicates_quant}.}
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
#' #dPCP analysis
#' results <- dPCP(sampleTable, system = "bio-rad", file.location = fileLoc,
#'                 eps = 200, minPts = 50, save.template = FALSE, rain = TRUE)
#'
#' plot(results, sample = 1, type = "dPCP")
#' }
#' @export

dPCP <- function(file, system = NULL, file.location = ".",
                 reference.quality = 0.5, sample.quality = 0.5,
                 eps = 200, minPts = 50, save.template = FALSE,
                 rain = TRUE) {

  if (!is.logical(rain)) stop("rain must be logical")

  #Read samples table
  samTable <- read_sampleTable(file = file, system = system,
                               file.location = file.location)

  #Read reference samples
  refSample <- read_reference(sample.table = samTable, system = system,
                              file.location = file.location,
                              reference.quality = reference.quality,
                              eps = eps, minPts = minPts)

  #Read samples
  samples <- read_sample(sample.table = samTable, system = system,
                         file.location = file.location,
                         sample.quality = sample.quality)

  #Reference DBSCAN clustering
  refSampleDB <- reference_dbscan(reference.subquality = refSample,
                                  sample.table = samTable, eps = eps,
                                  minPts = minPts,
                                  save.template = save.template)


  #Predict clusters centers position
  centers <- centers_data(sample.subquality = samples, sample.table = samTable,
                          referenceDB = refSampleDB)


  #Fuzzy c-means clustering
  cmeansclus <- cmeans_clus(centers.data = centers)

  if (isTRUE(rain)) {
    clustering <- rain_reclus(cmeans.cluster = cmeansclus)
  } else {
    clustering <- cmeansclus
  }

  #Quantification
  target.quant <- target_quant(data.cluster = clustering,
                               sample.table = samTable)

  #Replicates pooling
  results <- replicates_quant(raw.results = target.quant,
                              sample.table = samTable)

  return.list <- list("referenceDB" = refSampleDB,
                      "samples" = lapply(seq_along(samples), function(x) {
                        list(
                          "quality" = samples[[x]]$quality,
                          "reference" = samTable$Reference[x],
                          "centers" = centers[[x]]$centers,
                          "data" = cbind.data.frame(
                            samples[[x]]$data,
                            "cmeans cluster" = cmeansclus[[x]]$data$cluster,
                            "final cluster" = clustering[[x]]$data$cluster),
                          "cmeans membership" = cmeansclus[[x]]$membership,
                          "raw results" = target.quant[[x]]$`raw results`)
                      }),
                      "results" = results
  )
  names(return.list$samples) <- names(samples)

  class(return.list) <- "dPCP"

  return(return.list)
}


#' Plot the results of dPCP analysis
#'
#' @param x an object of class \code{dPCP}
#' @param ... Arguments to be passed to methods
#' @param  sample 'all' to show all samples, or a numeric vector indicating
#'   the row number of samples in the sample table.
#' @param  reference 'all' to show all reference samples, or a character vector
#'   with chip ID (Thermo Fisher) or the file name (Bio-rad) of reference
#'   samples to be showed.
#' @param type string. Type of plot to be showed. Available plots:
#'   'reference dbscan', 'centers', 'cmeans', 'rain', 'dPCP'.
#' @export

plot.dPCP <- function(x, ..., sample = "all", reference = "all",
                      type = "dPCP") {

  if (all(type != c("reference dbscan", "centers", "cmeans", "rain", "dPCP")))
    stop("'type' must be one of following string:
         'reference dbscan', 'centers', 'cmeans', 'rain', 'dPCP'")

  if (type == "reference dbscan") {
    plot.reference_dbscan(x$referenceDB, reference = reference)
  } else {
    data <- x$samples
    if (type == "centers") {
      class(data) <- "centers_data"
    }
    if (type == "cmeans") {
      class(data) <- c("cmeans_clus", "dPCP")
    }
    if (any(type == c("rain", "dPCP"))) {
      class(data) <- c("rain_reclus", "dPCP")
    }

    graphics::plot(x = data, sample = sample)
  }
}

