#' Calculation of targets concentration.
#'
#' This function calculates the concentration of the targets according to the
#' Poisson distribution.
#' @param data.cluster an object of class \code{rain_reclus} or
#'   \code{cmeans_clus}.
#' @inheritParams read_sample
#' @return An object of class \code{target_quant} containing a sublist for
#'   each sample. Each sublist has the following components:
#'   \item{quality}{quality threshold used in \code{\link{read_sample}}.}
#'   \item{reference}{reference ID.}
#'   \item{raw results}{a data frame with the results of the quantification.}
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
#' #Fuzzy c-means clustering
#' cmclus <- cmeans_clus(cent)
#'
#' #Rain classification.
#' rainclus <- rain_reclus(cmclus)
#'
#' #Quantification
#' quantcm <- target_quant(cmclus, sample.table)
#' quant <- target_quant(rainclus, sample.table)
#' }
#' @export

target_quant <- function(data.cluster, sample.table) {

  if (all(class(data.cluster) != c("rain_reclus", "cmeans_clus")))
    stop("data.cluster must be an object of class rain_reclus or cmeans_clus")

  if (!inherits(sample.table,"sample_table"))
    stop("'sample.table' must be an object of class sample_table")

  targetQuant <- lapply(seq_along(data.cluster), function(x) {
    #Calculate the number of expected clusters
    mrks <- sample.table$No.of.targets[x]
    expectedClus <- rownames(data.cluster[[x]]$centers)

    #Count number of elements in each cluster
    sampleClus <- as.data.frame(table(data.cluster[[x]]$data[, 3]),
                                stringsAsFactors = FALSE)

    colnames(sampleClus) <- c("cluster", "size")

    #Find empty clusters and add 0 number of elements
    emptyClus <- expectedClus[!expectedClus %in% sampleClus[, 1]]
    emptyClusTab <- cbind("cluster" = emptyClus,
                          "size" = rep(0, length(emptyClus)))

    #Organize all clusters and their sizes in a table and count the total
    #number of elements
    ClusCount <- rbind.data.frame(sampleClus, emptyClusTab)
    ClusCount <- ClusCount[order(match(ClusCount$cluster, expectedClus)), ]
    ClusCount$size <- as.numeric(ClusCount$size)
    totalwells <- sum(ClusCount$size)

    #Check the number of targets and consequently count the number of positive
    #wells for each cluster.
    if (mrks == 1) {
      N.pos <- ClusCount$size[2]
      target <- ClusCount[2, 1]

      results <- cbind.data.frame(
        "Sample" = rep(names(data.cluster)[x], length(target)), target, N.pos,
        totalwells)
    } else if (mrks == 2) {
      target1 <- sum(ClusCount$size[c(2, 4)])
      target2 <- sum(ClusCount$size[c(3, 4)])
      N.pos <- c(target1, target2)

      target <- ClusCount[2:3, 1]

      results <- cbind.data.frame(
        "Sample" = rep(names(data.cluster)[x], length(target)), target, N.pos,
        totalwells)
    } else if (mrks == 3) {
      famtarget <- sum(ClusCount$size[c(2, 5, 6, 8)])
      mixtarget <- sum(ClusCount$size[c(3, 5, 7, 8)])
      victarget <- sum(ClusCount$size[c(4, 6, 7, 8)])
      N.pos <- c(famtarget, mixtarget, victarget)

      target <- ClusCount[2:4, 1]

      results <- cbind.data.frame(
        "Sample" = rep(names(data.cluster)[x], length(target)), target, N.pos,
        totalwells)
    } else if (mrks == 4) {
      famtarget <- sum(ClusCount$size[c(2, 6, 7, 9, 12, 13, 14, 16)])
      target3target <- sum(ClusCount$size[c(3, 6, 8, 10, 12, 13, 15, 16)])
      target4target <- sum(ClusCount$size[c(4, 7, 8, 11, 12, 14, 15, 16)])
      victarget <- sum(ClusCount$size[c(5, 9, 10, 11, 13, 14, 15, 16)])
      N.pos <- c(famtarget, target3target, target4target, victarget)

      target <- ClusCount[2:5, 1]

      results <- cbind.data.frame(
        "Sample" = rep(names(data.cluster)[x], length(target)), target, N.pos,
        totalwells)
    }

    #Calculate the results  assuming a Poisson distribution for the number of
    #copies in each partition.

    if (data.cluster[[x]]$quality == "Defined by Bio-Rad") {
      partition.volume <- 0.00085
    } else if (is.numeric(data.cluster[[x]]$quality)) {
      partition.volume <- 0.000755
    } else {
      partition.volume <- as.numeric(
        stringr::str_split(data.cluster[[x]]$quality, ": ")[[1]][2])
    }

    pos_rate <- results$N.pos / results$totalwells

    ci <- lapply(seq(length(pos_rate)), function(x) {

      if (results$N.pos[x] == 0) {
        pois <- exactci::poisson.exact(
          results$N.pos[x], results$totalwells[x], alternative = "two.sided",
          tsmethod = "central", midp = TRUE)

        c(pois$conf.int[1], pois$conf.int[2])

      } else if (results$N.pos[x] < 100) {
        pois <- exactci::poisson.exact(
          results$N.pos[x], results$totalwells[x], alternative = "two.sided",
          tsmethod = "central", midp = FALSE)

        c(pois$conf.int[1], pois$conf.int[2])

      } else {

        c(pos_rate[x] - 1.96 * sqrt(pos_rate[x] * (1 - pos_rate[x])/
                                      results$totalwells[x]),
          pos_rate[x] + 1.96 * sqrt(pos_rate[x] * (1 - pos_rate[x])/
                                      results$totalwells[x]))
      }
    })

    ci <- do.call(rbind, ci)

    lower_CI_pos_rate <- ci[, 1]

    upper_CI_pos_rate <- ci[, 2]

    lambda <- - log(1 - pos_rate)

    lower_CI_lambda <- - log(1 - lower_CI_pos_rate)

    upper_CI_lambda <- - log(1 - upper_CI_pos_rate)

    concentration <- lambda / partition.volume

    lower_CI_concentration <- lower_CI_lambda / partition.volume

    upper_CI_concentration <- upper_CI_lambda / partition.volume

    conc.dil <- concentration / sample.table$Dilution[x]

    lower_CI_conc.dil <- lower_CI_concentration / sample.table$Dilution[x]

    upper_CI_conc.dil <- upper_CI_concentration / sample.table$Dilution[x]

    spread <- cbind(abs(1 - (lower_CI_lambda)), abs(1 - (upper_CI_lambda)))

    precision <- apply(spread, 1, max) * 100

    results <- cbind.data.frame(
      results, lambda, lower_CI_lambda, upper_CI_lambda, concentration,
      lower_CI_concentration, upper_CI_concentration, conc.dil,
      lower_CI_conc.dil, upper_CI_conc.dil, round(precision, 3),
      rep(sample.table$Dilution[[x]], sample.table$No.of.targets[[x]]))

    names(results) <- c(
      "Sample", "Target", "Positive reactions", "Total reactions", "lambda",
      "Lower CI lambda", "Upper CI lambda", "Copies/ul", "Lower CI copies/ul",
      "Upper CI copies/ul", "Copies/ul at sample dilution",
      "Lower CI copies/ul at sample dilution",
      "Upper CI copies/ul at sample dilution", "Precision %", "Dilution")

    list("quality" = data.cluster[[x]]$quality,
         "reference" = data.cluster[[x]]$reference, "raw results" = results)
  })
  names(targetQuant) <- names(data.cluster)

  class(targetQuant) <- "target_quant"

  return(targetQuant)
}


#' Calculation of targets concentration, pooling the sample replicates
#'
#' This function calculates the concentration of the targets, combining the
#' results of the replicates of each sample.
#' @param raw.results an object of class \code{target_quant}, inherited
#'   from \code{\link{target_quant}}.
#' @inheritParams read_sample
#' @return An object of class \code{replicates_quant} containing a sublist for
#'   every sample. Each sublist has the following components:
#'   \item{quality}{quality threshold used in \code{\link{read_sample}}.}
#'   \item{reference}{reference ID.}
#'   \item{raw results}{a data frame with the results of quantification.}
#'   \item{replicates results}{a data frame with the results of quantification
#'   of pooled replicates.}
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
#'                    file.location = fileLoc)
#'
#' #Reference DBSCAN clustering
#' dbref <- reference_dbscan(ref, sample.table, save.template = FALSE)
#'
#' #Predict position of clusters centroid from reference DBSCAN results
#' cent <- centers_data(samp, sample.table,dbref)
#'
#' #Fuzzy c-means clustering
#' cmclus <- cmeans_clus(cent)
#'
#' #Rain classification.
#' rainclus <- rain_reclus(cmclus)
#'
#' #Quantification
#' quantcm <- target_quant(cmclus, sample.table)
#' quant <- target_quant(rainclus, sample.table)
#'
#' #Replicates pooling
#' rep.quant <- replicates_quant(quant, sample.table)
#' }
#' @export

replicates_quant <- function(raw.results, sample.table) {

  if (!inherits(raw.results, "target_quant"))
    stop("raw.results must be an object of class target_quant")

  if (!inherits(sample.table, "sample_table"))
    stop("'sample.table' must be an object of class sample_table")

  if (raw.results[[1]]$quality == "Defined by Bio-Rad") {
    partition.volume <- 0.00085
  } else if (is.numeric(raw.results[[1]]$quality)) {
    partition.volume <- 0.000755
  } else {
    partition.volume <- as.numeric(
      stringr::str_split(raw.results[[1]]$quality, ": ")[[1]][2])
  }

  if (any(duplicated(sample.table$Sample.name))) {

    replicates <- unique(
      sample.table$Sample.name[duplicated(sample.table$Sample.name)])

    repl.data <- lapply(replicates, function(x) {

      repl.index <- which(sample.table$Sample.name %in% x)
      n.target <- sample.table$No.of.targets[repl.index[1]]
      tot.rows <- length(repl.index) * n.target

      raw.res <- lapply(repl.index, function(y) {
        if (any(is.na(raw.results[[y]]$`raw results`))) {
          raw.results[[y]]$`raw results`[
            is.na(raw.results[[y]]$`raw results`)] <- 0
          raw.results[[y]]$`raw results`
        } else {
          raw.results[[y]]$`raw results`
        }
      })

      raw.res <- do.call(rbind.data.frame, raw.res)

      repl.calc <- lapply(seq(n.target), function(z) {

        sing.target <- raw.res[seq(z, tot.rows, n.target), c(3, 4)]

        pos_tot <- sum(sing.target[, 1])

        total_tot <- sum(sing.target[, 2])

        pos_rate <- pos_tot / total_tot

        if (pos_tot == 0) {
          pois <- exactci::poisson.exact(
            pos_tot, total_tot, alternative = "two.sided",
            tsmethod = "central", midp = TRUE)

          lower_CI_pos_rate <- pois$conf.int[1]

          upper_CI_pos_rate <- pois$conf.int[2]

        } else if (pos_tot < 100) {
          pois <- exactci::poisson.exact(
            pos_tot, total_tot, alternative = "two.sided",
            tsmethod = "central", midp = FALSE)

          lower_CI_pos_rate <- pois$conf.int[1]

          upper_CI_pos_rate <- pois$conf.int[2]

        } else {

          lower_CI_pos_rate <- pos_rate -
            1.96 * sqrt(pos_rate * (1 - pos_rate) / total_tot)

          upper_CI_pos_rate <- pos_rate +
            1.96 * sqrt(pos_rate * (1 - pos_rate) / total_tot)
        }

        lambda <- - log(1 - pos_rate)

        lower_CI_lambda <- - log(1 - lower_CI_pos_rate)

        upper_CI_lambda <- - log(1 - upper_CI_pos_rate)

        concentration <- lambda / partition.volume

        lower_CI_concentration <- lower_CI_lambda / partition.volume

        upper_CI_concentration <- upper_CI_lambda / partition.volume

        spread <- cbind(abs(1 - (lower_CI_lambda)), abs(1 - (upper_CI_lambda)))

        precision <- apply(spread, 1, max) * 100

        pool <- cbind.data.frame(
          "Copies/ul" = concentration,
          "Lower CI copies/ul" = lower_CI_concentration,
          "Upper CI copies/ul" = upper_CI_concentration,
          "Precision %" = precision, "No of replicates" = length(repl.index))
      })
      repl.calc <- do.call(rbind.data.frame, repl.calc)

      repl.calc <- cbind.data.frame(
        "Sample" = rep(sample.table$Sample.name[repl.index[1]],
                       sample.table$No.of.targets[repl.index[1]]),
        "Target" = raw.results[[repl.index[1]]]$`raw results`$Target,
        repl.calc)

      qual <- sapply(repl.index, function(k) {
        raw.results[[k]]$quality
      })

      ref <- sapply(repl.index, function(j) {
        raw.results[[j]]$reference
      })

      repl <- sapply(repl.index, function(i) {
        names(raw.results[i])
      })

      list("quality" = qual, "reference" = ref, "replicates" = repl,
           "raw results" = raw.res, "replicates results" = repl.calc)
    })
    names(repl.data) <- replicates

    repl.index <- match(replicates, sample.table$Sample.name)

    not.repl <- raw.results[which(!sample.table$Sample.name %in% replicates)]
    not.repl.index <- which(!sample.table$Sample.name %in% replicates)

    all.samples <- append(not.repl, repl.data)
    all.indexes <- append(not.repl.index, repl.index)

    all.samples <- all.samples[order(all.indexes)]
  } else {
    all.samples <- raw.results
  }
  class(all.samples) <- "replicates_quant"

  return(all.samples)
}


#' Export dPCP analysis results to a pdf report
#'
#' This function generates a pdf report of the dPCP analysis.
#' @inheritParams manual_correction
#' @inheritParams dPCP
#' @return A pdf file with the information and results of the dPCP analysis.
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
#'                 eps = 200, minPts = 50, save.template = FALSE,
#'                 rain = TRUE)
#'
#' report_dPCP(results, filename = "dPCRproject_1")
#' }
#' @export

report_dPCP <- function(data, filename, sample = "all") {

  if (all(class(data) != "dPCP"))
    stop("data must be an object of class dPCP")

  if (!is.character(filename) ||
      (is.character(filename) & length(filename) > 1))
    stop("filename must be a character string.")

  if (all(sample != "all") & !is.numeric(sample))
    stop("sample must be `all` or a numeric vector with row numbers of
         sample.table indicating samples to be showed")

  if (is.numeric(sample) &
      (any(sample < 1) | any(sample > length(data$samples))))
    stop("Undefined sample number selected")

  if (all(sample == "all")) {
    export <- seq_along(data$sample)
  } else {
    export <- sort(sample)
  }

  step1 <- lapply(seq_along(export), function(x) {

    tb1 <- ggpubr::ggtexttable(
      cbind(data$sample[[x]]$`raw results`[, 2:4],
            round(data$sample[[x]]$`raw results`[, 5], 5),
            paste0(round(data$sample[[x]]$`raw results`[, 8], 2),
                   " (", round(data$sample[[x]]$`raw results`[, 9], 2), " - ",
                   round(data$sample[[x]]$`raw results`[, 10], 2), ")"
            ),
            paste0(round(data$sample[[x]]$`raw results`[, 11], 2),
                   " (", round(data$sample[[x]]$`raw results`[, 12], 2), " - ",
                   round(data$sample[[x]]$`raw results`[, 13], 2), ")"
            ),
            round(data$sample[[x]]$`raw results`[, 14], 2)),
      rows = NULL,
      cols = c("Target", "Positive\nreactions", "Total\nreactions", "lambda",
               "Copies/ul\n(95% CI)", "Copies/ul at sample\ndilution (95% CI)",
               "Precision\n%"),
      theme = ggpubr::ttheme(base_size = 10))

    tx <- paste(
      "Quality threshold:", data$sample[[x]]$quality,
      "\nDilution:", data$sample[[x]]$`raw results`$Dilution,
      "\nReference:",
      paste0(data$sample[[x]]$reference, " (quality: ",
             data$referenceDB[[match(data$sample[[x]]$reference,
                                     names(data$referenceDB))]]$quality,
             ", ", "eps: ",
             data$referenceDB[[match(data$sample[[x]]$reference,
                                     names(data$referenceDB))]]$dbscan$eps,
             ", ", "minPts: ",
             data$referenceDB[[match(data$sample[[x]]$reference,
                                     names(data$referenceDB))]]$dbscan$minPts,
             ")"
      )
    )

    text.p <- ggpubr::ggparagraph(text = tx, color = "black", size = 12)


    cluscolors <- c(
      "gray70", "#004949", "#ff6db6", "#009292", "#ffb6db", "#490092",
      "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", "#920000", "#924900",
      "#db6d00", "#24ff24", "#ffff6d",
      "#000000")[1:nrow(data$sample[[x]]$centers)]
    names(cluscolors) <- row.names(data$sample[[x]]$centers)

    graph.p <- ggplot(as.data.frame(data$sample[[x]]$data),
                      aes_string(x = "Vic", y = "Fam")) +

      geom_point(
        size = 0.5,
        aes(color = data$sample[[x]]$data[, ncol(data$sample[[x]]$data)])) +

      scale_colour_manual(values = cluscolors) +

      labs(color = "Clusters", title = names(data$sample)[x]) +

      theme(
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black"),
        plot.title = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 10)
      ) +

      guides(colour = guide_legend(override.aes = list(size = 2)))

    tot <- ggpubr::ggarrange(graph.p, text.p, tb1, ncol = 1, nrow = 3,
                             heights = c(7, 2, 4))

  })

  ggpubr::ggexport(plotlist = step1, filename = paste0(filename, ".pdf"),
                   width = 8, height = 11, res = 300)
}


#' Export dPCP analysis results to a csv file
#'
#' This function exports dPCP analysis results to a csv file.
#' @param data an object of class \code{dPCP}, \code{target_quant} or
#'   \code{replicates_quant}.
#' @inheritParams manual_correction
#' @return A csv file with the information and results of dPCP analysis.
#' @examples
#' \donttest{
#' library(dPCP)
#'
#' #Find path of sample table and location of reference and input files
#' sampleTable <- system.file("extdata", "Template_sampleTable.csv",
#'                      package = "dPCP")
#'
#' fileLoc <- system.file("extdata",package = "dPCP")
#'
#' #dPCP analysis
#' results <- dPCP(sampleTable, system = "bio-rad", file.location = fileLoc,
#'                 eps = 200, minPts = 50, save.template = FALSE,
#'                 rain = TRUE)
#'
#' export_csv(results, filename = "dPCRproject_1")
#' }
#' @export

export_csv <- function(data, filename) {

  if (!is.character(filename) ||
      (is.character(filename) & length(filename) > 1))
    stop("filename must be a character string.")

  if (any(class(data) == c("target_quant", "replicates_quant"))) {

    data.results <- data
    exp.ref <- NULL

  } else if (any(class(data) == "dPCP")) {

    data.results <- data$results
    data.reference <- data$referenceDB

    c.data <- lapply(seq_along(data.reference), function(x) {

      paste(names(data.reference)[x], data.reference[[x]]$quality,
            data.reference[[x]]$dbscan$eps,
            data.reference[[x]]$dbscan$minPts, sep = ",")
    })

    c.name <- paste("Reference", "Quality", "DBSCAN eps", "DBSCAN minPts",
                    sep = ",")

    all.ref <- append(c.name, c.data)

    exp.ref <- paste(unlist(all.ref), sep = "\n")

  } else {
    stop("data must be an object of class target_quant, replicates_quant or
         dPCP")}

  exp.data <- lapply(seq_along(data.results), function(x) {

    n.repl <- length(data.results[[x]]$quality)
    n.target <- nrow(data.results[[x]]$`raw results`) /
      length(data.results[[x]]$quality)

    repetition <- rep(n.target, n.repl)

    raw.data <- cbind.data.frame(
      data.results[[x]]$`raw results`,

      if (is.character(data.results[[x]]$quality)) {
        "Quality" = rep(as.character(sapply(data.results[[x]]$quality, "[")),
                        repetition)
      } else {
        "Quality" = rep(as.numeric(sapply(data.results[[x]]$quality, "[")),
                        repetition)
      },

      "Reference" = rep(as.character(sapply(data.results[[x]]$reference, "[")),
                        repetition)
    )

    raw.data <- do.call(paste, c(lapply(raw.data, "["), sep = ","))


    if (length(data.results[[x]]) == 3) {

      exp.repl <- cbind.data.frame(
        data.results[[x]]$`raw results`[, c(1, 2, 11:14)],
        "No of replicates" = rep(1, nrow(data.results[[x]]$`raw results`))
      )
    } else {
      exp.repl <- data.results[[x]]$`replicates results`
    }

    exp.repl <- do.call(paste, c(lapply(exp.repl, "["), sep = ","))

    list("data_results" = raw.data, "results" = exp.repl)
  })
  names(exp.data) <- names(data.results)

  raw.name <- paste(
    "Sample", "Target", "Positive reactions", "Total reactions", "lambda",
    "Lower CI lambda", "Upper CI lambda", "Copies/microl",
    "Lower CI copies/microl", "Upper CI copies/microl",
    "Copies/microl at sample dilution",
    "Lower CI copies/microl at sample dilution",
    "Upper CI copies/microl at sample dilution",
    "Precision %", "Dilution", "Quality", "Reference", sep = ",")

  all.raw <- append(raw.name, lapply(exp.data, function(x) {
    x$data_results
  }))
  exp.raw <- paste(unlist(all.raw), sep = "\n")

  repl.name <- paste("Sample", "Target", "Copies/microl",
                     "Lower CI copies/microl", "Upper CI copies/microl",
                     "Precision %", "No of replicates", sep = ",")

  all.repl <- append(repl.name, lapply(exp.data, function(x) {
    x$results
  }))
  exp.repl <- paste(unlist(all.repl), sep = "\n")


  if (!is.null(exp.ref)) {

    all.exp <- rlist::list.append("Reference samples", exp.ref, "",
                                  "Results", exp.raw, "",
                                  "Replicates results", exp.repl)
  } else {

    all.exp <- rlist::list.append("Results", exp.raw, "", "Replicates results",
                                  exp.repl)
  }

  all.exp <- strsplit(all.exp, ",")

  all.exp <- sapply(all.exp, function(x) {
    if (length(x) > 10) {
      x[5:16] <- gsub("\\.", ",", x[5:16])
    } else if (length(x) > 6) {
      x[3:6] <- gsub("\\.", ",", x[3:6])
    }
    paste(x, collapse = "\t")
  })

  utils::write.table(all.exp, paste0(filename, ".csv"), row.names = FALSE,
                     col.names = FALSE, sep = "\t", na = "", quote = TRUE)
}
