#' Cluster analysis with fuzzy c-means algorithm
#'
#' This function carries out the c-means cluster analysis, using the centroids
#' position as initial values for cluster centers.
#' @param centers.data an object of class \code{centers_data}, inherited
#'   from \code{\link{centers_data}}.
#' @rdname plot.cmeans_clus
#' @return An object of class \code{cmeans_clus} containing a sublist for
#'   each sample. Each sublist has the following components:
#'   \item{quality}{quality threshold used in \code{\link{read_sample}}.}
#'   \item{reference}{reference ID.}
#'   \item{centers}{a data frame with the centroids coordinates.}
#'   \item{data}{a data frame with the fluorescence intensities and clusters
#'   name.}
#'   \item{membership}{a matrix with the membership values of the data elements
#'   to the clusters. See also \code{\link[e1071]{cmeans}}}
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
#' plot(cmclus, sample = "all")
#' }
#' @export

cmeans_clus <- function(centers.data) {

  if (class(centers.data) != "centers_data")
    stop("centers.data must be an object of class centers_data")

  cClus <- lapply(centers.data, function(x) {

    cm1 <- e1071::cmeans(x$data, centers = x$centers, iter.max = 3,
                         dist = "euclidean", method = "ufcl")

    #Use row names of centers table as cluster name
    clusname <- rownames(x$centers)[sort(unique(cm1$cluster))]
    cmclus <- factor(cm1$cluster)
    levels(cmclus) <- clusname

    cm.first <- cbind.data.frame(x$data, "cluster" = cmclus)

    sep.data <- split.data.frame(cm.first[,c(1,2)], cm.first[,3])

    clus.cent <- lapply(sep.data, colMeans)

    clus.cent <- do.call(rbind, clus.cent)

    noclus.cent <- subset(x$centers,
                          !(rownames(x$centers) %in% row.names(clus.cent)))

    new.cent <- rbind.data.frame(clus.cent, noclus.cent)

    new.cent <- new.cent[order(match(rownames(new.cent),
                                     rownames(x$centers))), ]

    cm <- e1071::cmeans(x$data, centers = new.cent, iter.max = 3,
                        dist = "euclidean", method = "ufcl")

    #Use row names of centers table as cluster name
    clusname <- rownames(x$centers)[sort(unique(cm$cluster))]
    cmclus <- factor(cm$cluster)
    levels(cmclus) <- clusname

    #Use row names of centers table as membership column names
    colnames(cm$membership) <- rownames(x$centers)

    clusters <- cbind.data.frame(x$data, "cluster" = cmclus)

    list("quality" = x$quality, "reference" = x$reference,
         "centers" = x$centers, "data" = clusters,
         "membership" = cm$membership)
  })

  class(cClus) <- "cmeans_clus"

  return(cClus)
}


#' Plot results of cmeans cluster analysis
#'
#' @param x an object of class \code{cmeans_clus}
#' @param ... Arguments to be passed to methods
#' @inheritParams dPCP
#' @export

plot.cmeans_clus <- function(x, ..., sample = "all") {

  if (all(sample != "all") & any(sample > length(x)))
    stop("sample must be `all` or a numeric vector indicating the row number of
         samples in sample.table")

  if (all(sample == "all")) {
    plotsample <- 1:length(x)
  } else {
    plotsample <- sort(sample)
  }

  cmeansplot <- lapply(plotsample, function(y) {

    if (all(class(x) == "cmeans_clus")) {
      cluster  <-  x[[y]]$data[, 3]
    } else if (all(class(x) == c("cmeans_clus", "dPCP"))) {
      cluster  <-  x[[y]]$data[, ncol(x[[y]]$data) - 1]
    }

    cluscolors <- c(
      "gray70", "#004949", "#ff6db6", "#009292", "#ffb6db", "#490092",
      "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", "#920000", "#924900",
      "#db6d00", "#24ff24", "#ffff6d",
      "#000000")[1:(nrow(x[[y]]$centers))]
    names(cluscolors) <- row.names(x[[y]]$centers)

    p <- ggplot(as.data.frame(x[[y]]$data),
                aes_string(x = "Vic", y = "Fam")) +

      geom_point(size = 1, aes(color = cluster)) +

      scale_colour_manual(values = cluscolors) +

      labs(color = "Clusters",
           title = names(x)[y]) +

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

    graphics::plot(p)
  })
  names(cmeansplot) <- names(x)[plotsample]

  return(cmeansplot)
}



#' Identification and clustering of "rain" data
#'
#' This function identifies the "rain" elements and re-clusters them using the
#' Mahalanobis distance. Each "rain" element is assigned to the cluster whose
#' Mahalanobis distance is the lowest.
#' @param cmeans.cluster an object of class \code{cmeans_clus}, inherited
#'   from \code{\link{cmeans_clus}}.
#' @rdname plot.rain_reclus
#' @return An object of class \code{rain_reclus} containing a sublist for
#'   each sample. Each sublist has the following components:
#'   \item{quality}{quality threshold used in \code{\link{read_sample}}.}
#'   \item{reference}{reference ID.}
#'   \item{centers}{a data frame with the centroids coordinates.}
#'   \item{data}{a data frame with the fluorescence intensities and clusters
#'   name.}
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
#' plot(rainclus, sample = "all")
#' }
#' @export

rain_reclus <- function(cmeans.cluster) {

  if (class(cmeans.cluster) != "cmeans_clus")
    stop("cmeans.cluster must be an object of class cmeans_clus")

  rainclus <- lapply(cmeans.cluster, function(x) {

    clus <- table(x$data$cluster)
    trueclus <- names(clus)[clus > 4]

    if (length(trueclus) <= 4) {
      probability <- 0.5
    } else if (length(trueclus) >= 5 & length(trueclus) <= 7) {
      probability <- 0.5
    } else {
      probability <- 0.5
    }


    max.mem <- apply(x$membership, 1, max)

    rain <- list("rain" = subset(x$data, max.mem < probability),
                 "norain" = subset(x$data, max.mem >= probability))

    mahala.dis <- if (nrow(rain$rain) < 2) {
      x$data
    }
    else {
      sapply(rownames(x$centers), function(y) {
        norain.sub <- subset(rain$norain[, c(1, 2)], rain$norain[, 3] == y)

        if (nrow(norain.sub) < 3 | rcond(stats::cov(norain.sub))  <= 1e-5) {
          eucl.dist <- raster::pointDistance(
            rain$rain[, c(1, 2)], x$centers[which(rownames(x$centers) == y), ],
            lonlat = FALSE)
        } else {
          stats::mahalanobis(rain$rain[, c(1, 2)], colMeans(norain.sub),
                             stats::cov(norain.sub), tol = 1e-5)
        }
      })
    }

    reClus <- if (!is.null(colnames(mahala.dis)) &
                  all(colnames(mahala.dis) %in%  c("Vic", "Fam", "cluster"))) {
      mahala.dis
    } else {
      rainClus <- cbind(
        rain$rain[, c(1, 2)],
        "cluster" = colnames(mahala.dis)[apply(mahala.dis, 1, which.min)])
      rbind(rain$norain, rainClus)
    }


    list("quality" = x$quality, "reference" = x$reference,
         "centers" = x$centers,
         "data" = reClus[order(as.numeric(rownames(reClus))), ])
  })

  class(rainclus) <- "rain_reclus"

  return(rainclus)
}


#' Plot results of "rain" analysis
#'
#' @param x an object of class \code{rain_reclus}
#' @param ... Arguments to be passed to methods
#' @inheritParams dPCP
#' @export

plot.rain_reclus <- function(x, ..., sample = "all") {

  if (all(sample != "all") & any(sample > length(x)))
    stop("sample must be `all` or a numeric vector indicating the row number of
         samples in sample.table")

  if (all(sample == "all")) {
    plotsample <- 1:length(x)
  } else {
    plotsample <- sort(sample)
  }


  dPCPplot <- lapply(plotsample, function(y) {

    if (all(class(x) == "rain_reclus")) {
      cluster  <-  x[[y]]$data[, 3]
    } else if (all(class(x) == c("rain_reclus", "dPCP"))) {
      cluster  <-  x[[y]]$data[, ncol(x[[y]]$data)]
    }

    cluscolors <- c(
      "gray70", "#004949", "#ff6db6", "#009292", "#ffb6db", "#490092",
      "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", "#920000", "#924900",
      "#db6d00", "#24ff24", "#ffff6d",
      "#000000")[1:nrow(x[[y]]$centers)]
    names(cluscolors) <- row.names(x[[y]]$centers)

    p <- ggplot(as.data.frame(x[[y]]$data),
                aes_string(x = "Vic", y = "Fam")) +

      geom_point(size = 1, aes(color = cluster)) +

      scale_colour_manual(values = cluscolors) +

      labs(color = "Clusters",
           title = names(x)[y]) +

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

    graphics::plot(p)
  })
  names(dPCPplot) <- names(x)[plotsample]

  return(dPCPplot)
}



#' Manual correction of dPCP cluster analysis
#'
#' This function builds an interactive app to manually correct the dPCP
#' cluster analysis.
#' @param data an object of class \code{dPCP}, inherited from
#'   \code{\link{dPCP}}.
#' @param filename character. File name (no extension) for csv and pdf files to
#'   create on disk.
#' @param  save.plot logical. If TRUE the plots are exported to a file.
#' @param  format a string indicating the file format for the export.
#'   Available formats: 'eps', 'ps', 'tex', 'pdf', 'jpeg', 'tiff', 'png',
#'   'bmp', 'svg', 'wmf'.
#' @param dpi numeric. Image resolution.
#' @return A Shiny session.
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
#' manual_correction(results, filename = "manual_dPCR", save.plot = FALSE)
#' }
#' @export
#' @import shiny

manual_correction <- function(data, filename, save.plot = FALSE,
                              format = "png", dpi = 300) {

  if (all(class(data) != c("dPCP")))
    stop("data must be an object of class dPCP")

  if (!is.character(filename) ||
      (is.character(filename) & length(filename) > 1))
    stop("filename must be a character string.")

  if (!is.logical(save.plot))
    stop("save.plot must be logical")

  if (all(format != c("eps", "ps", "tex", "pdf", "jpeg", "tiff", "png", "bmp",
                      "svg", "wmf")))
    stop("Uknown format")

  if (!is.numeric(dpi) || length(dpi) != 1) stop("dpi must be a numeric value")

  if ((is.character(data$samples[[1]]$quality)) &
      (data$samples[[1]]$quality == "Defined by Bio-Rad")) {
    partition.volume <- 0.00085
  } else if ((is.character(data$samples[[1]]$quality)) &
             (data$samples[[1]]$quality != "Defined by Bio-Rad")) {
    partition.volume <- as.numeric(
      str_split(data$samples[[1]]$quality, ": ")[[1]][2])
  } else {
    partition.volume <- 0.000755
  }

  clus.char <- lapply(data$samples, function(x) {
    x$data[, ncol(x$data)] <- as.character(x$data[, ncol(x$data)])
    return(x)
  })

  data$samples <- clus.char

  # Define UI
  ui <- fluidPage(

    h1("digital PCR Cluster Predictor"),


    fluidRow(
      #Set dropdown bar to select samples
      column(
        width = 10,
        selectInput(inputId = "sample", label = "Sample",
                    choices = names(data$samples))),
      #Set export action button
      column(
        width = 2,
        actionButton(inputId = "export", label = "Export all",
                     icon = icon("download"),
                     width = "100%",
                     style = "border-color: black; font-weight: bold"))
    ),

    fluidRow(
      #Set output graph
      column(
        width = 8,
        plotOutput("finalplot", dblclick = "finalplot_dblclick",
                   brush = brushOpts(id =  "finalplot_brush", fill = "white",
                                     stroke = "black", resetOnNew = TRUE))),
      #Set clusters action buttons
      column(width = 3,
             h5("Cluster", style = "text-align: center; font-weight: bold"),
             uiOutput("cluster_buttons")
      )
    ),
    shinyjs::useShinyjs(),
    #Set undo action button
    actionButton("undo", "Undo", width = "22%",
                 style = "border-color: black; font-weight: bold"),

    #Set reset action button
    actionButton("reset", "Reset", width = "22%",
                 style = "border-color: black; font-weight: bold"),

    #Set output table for results
    tableOutput("results")
  )

  server <- function(input, output) {

    #Create reactive variables for fluorence data, clusters and results
    cluster.data <- reactive({
      data$samples[[match(input$sample, names(data$samples))]]$data
    })

    centers.data <- reactive({
      data$samples[[match(input$sample, names(data$samples))]]$centers
    })

    results <- reactive({
      data$samples[[match(
        input$sample, names(data$samples))]]$`raw results`[, 2:15]
    })

    #Set reactive values for provisional changes (last modifications)
    values <- reactiveValues(
      new = NULL, #Store provisional cluster data
      results = NULL, #Store provisional results
      rows.name = NULL, #Store rows names of brushed points
      rows.index = NULL, #Store rows index of brushed points
      clus = NULL) #Store name of new cluster


    #Set reactive values to store definitive changes
    manual.mod <- reactiveValues(
      samples = c(), #Track sequence of visualized samples
      clus = lapply(seq_along(data$samples), function(x) {
        data$samples[[x]]$data
      }), #Store cluster data
      res = lapply(seq_along(data$samples), function(x) {
        data$samples[[x]]$`raw results`[, 2:15]
      })#Store results
    )

    #Set reactive values to store clusters action buttons and their actions
    in.clus  <- reactiveValues(buttons = list(), actions = 0)

    #Set action for dropdown bar
    observeEvent(input$sample, {
      #Track sample choise
      manual.mod$samples <- c(manual.mod$samples,
                              match(input$sample, names(data$samples)))

      #When sample changes and manual correction has been done in the previous
      #sample save changes in manual.mod
      if (length(manual.mod$samples) > 1 & !is.null(values$new)) {

        manual.mod$clus[[manual.mod$samples[length(manual.mod$samples)-1]]] <-
          values$new[, c(1:ncol(cluster.data()), ncol(values$new))]

        manual.mod$res[[manual.mod$samples[length(manual.mod$samples)-1]]] <-
          utils::tail(values$results, nrow(data$samples[[
            manual.mod$samples[length(manual.mod$samples)-1]]]$`raw results`))
      }

      values$new <- NULL #Reset values$new when a new sample is called
      values$results <- NULL #Reset values$results when a new sample is called
      #Reset cluster buttons counter when a new sample is called
      in.clus$count <- 0

      in.clus$buttons <- NULL #Reset action buttons when a new sample is called

      #Create an action button for each predited cluster
      lapply(seq(nrow(centers.data())), function(x) {

        buttons.colors <- c(
          "gray70", "#004949", "#ff6db6", "#009292", "#ffb6db", "#490092",
          "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", "#920000", "#924900",
          "#db6d00", "#24ff24", "#ffff6d",
          "#000000")[1:nrow(centers.data())]

        background.color <- paste0(
          "background-color:", buttons.colors[x], "; border-color: black", ";
          font-weight: bold", "; color: grey")


        in.clus$buttons[[x]] <- actionButton(
          inputId = paste0("button", x),
          label = rownames(centers.data())[x],
          style = background.color,
          width = "100%")
      })
      #Add the dynamic buttons into a single fluidRow
      output$cluster_buttons  <- renderUI({
        do.call(fluidRow, in.clus$buttons)
      })
    })

    #Create output graph
    #Set reactive values for zoom coordinates
    ranges <- reactiveValues(x = NULL, y = NULL)

    #if values$new is not NULL (new manual correction done), use values$new
    #for the plot.
    #If values$new is NULL check for saved manual correction data and use it
    #for the plot.
    dataplot <- reactive({
      if (!is.null(values$new)) {
        values$new
      } else {
        manual.mod$clus[[match(input$sample, names(data$samples))]]
      }
    })

    output$finalplot <- renderPlot( {

      cluscolors <- c(
        "gray70", "#004949", "#ff6db6", "#009292", "#ffb6db", "#490092",
        "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", "#920000", "#924900",
        "#db6d00", "#24ff24", "#ffff6d",
        "#000000")[1:nrow(centers.data())]
      names(cluscolors) <- rownames(centers.data())
      ggplot(dataplot(), aes_string(x = "Vic", y = "Fam")) +

        geom_point(size = 1.5,
                   aes(color = as.factor(dataplot()[, ncol(dataplot())]))) +

        scale_colour_manual(values = cluscolors, drop = FALSE) +

        theme(
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text = element_text(size = 14),
          legend.position = "none"
        ) +

        coord_cartesian(xlim = ranges$x, ylim = ranges$y) #Zoom coordinates
    })

    # When a double-click happens, check if there's a brush on the plot.
    # If so, zoom to the brush bounds; if not, reset the zoom.
    observeEvent(input$finalplot_dblclick, {

      brush <- input$finalplot_brush

      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })


    #Set action for each cluster action button
    lapply(seq(16), function(x) {

      #Change color according to manual correction
      observeEvent(input[[paste0("button", x)]], {

        #Store rows name of brushed points
        values$rows.name <- rownames(brushedPoints(dataplot(),
                                                   input$finalplot_brush))
        #Store the name of new selected cluster
        values$clus <- rep(in.clus$buttons[[x]]$children[[1]][[2]],
                           length(values$rows.name))

        #If values$new is NULL take saved manual correction data
        if (is.null(values$new)) {
          values$new <- manual.mod$clus[[match(input$sample,
                                               names(data$samples))]]
        }

        #If points are brushed recalculated results and graph according to
        #cluster button input
        if (length(values$rows.name) > 0) {

          #Track clusters action buttons actions
          in.clus$count <- in.clus$count + 1

          #Find indexes of brushed points
          values$rows.index <- match(values$rows.name, rownames(values$new))
          #Create a new column to store new clusters values
          values$new[, ncol(values$new) + 1] <- values$new[, ncol(values$new)]
          #Convert cluster of brushed points to new selected cluster
          values$new[values$rows.index, ncol(values$new)] <- values$clus


          #Recalculate results according to manual correction
          mrks <- log2(nrow(centers.data()))
          expectedClus <- rownames(centers.data())

          sampleClus <- as.data.frame(table(
            values$new[, ncol(values$new)]), stringsAsFactors = FALSE)

          colnames(sampleClus) <- c("cluster", "size")

          emptyClus <- expectedClus[!expectedClus %in% sampleClus[, 1]]
          emptyClusTab <- cbind("cluster" = emptyClus,
                                "size" = rep(0, length(emptyClus)))

          ClusCount <- rbind.data.frame(sampleClus, emptyClusTab)
          ClusCount <- ClusCount[order(match(ClusCount$cluster,
                                             expectedClus)), ]
          ClusCount$size <- as.numeric(ClusCount$size)
          totalwells <- sum(ClusCount$size)

          if (mrks == 1) {
            N.pos <- ClusCount$size[2]
            target <- ClusCount[2, 1]

            results <- cbind.data.frame(target, N.pos, totalwells)
          }else if (mrks == 2) {
            target1 <- sum(ClusCount$size[c(2, 4)])
            target2 <- sum(ClusCount$size[c(3, 4)])
            N.pos <- c(target1, target2)

            target <- ClusCount[2:3, 1]

            results <- cbind.data.frame(target, N.pos, totalwells)
          } else if (mrks == 3) {
            famtarget <- sum(ClusCount$size[c(2, 5, 6, 8)])
            mixtarget <- sum(ClusCount$size[c(3, 5, 7, 8)])
            victarget <- sum(ClusCount$size[c(4, 6, 7, 8)])
            N.pos <- c(famtarget, mixtarget, victarget)

            target <- ClusCount[2:4, 1]

            results <- cbind.data.frame(target, N.pos, totalwells)
          } else if (mrks == 4) {
            famtarget <- sum(ClusCount$size[c(2, 6, 7, 9, 12, 13, 14, 16)])
            target3target <- sum(
              ClusCount$size[c(3, 6, 8, 10, 12, 13, 15, 16)])
            target4target <- sum(
              ClusCount$size[c(4, 7, 8, 11, 12, 14, 15, 16)])
            victarget <- sum(
              ClusCount$size[c(5, 9, 10, 11, 13, 14, 15, 16)])
            N.pos <- c(famtarget, target3target, target4target, victarget)

            target <- ClusCount[2:5, 1]

            results <- cbind.data.frame(target, N.pos, totalwells)
          }

          pos_rate <- results$N.pos / results$totalwells

          ci <- lapply(seq(length(pos_rate)), function(x) {

            if (results$N.pos[x] == 0) {
              pois <- exactci::poisson.exact(
                results$N.pos[x], results$totalwells[x],
                alternative = "two.sided", tsmethod = "central", midp = TRUE)

              c(pois$conf.int[1], pois$conf.int[2])

            } else if (results$N.pos[x] < 100) {
              pois <- exactci::poisson.exact(
                results$N.pos[x], results$totalwells[x],
                alternative = "two.sided", tsmethod = "central", midp = FALSE)

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

          #lower_CI_pos_rate <- pos_rate - 1.96 *
          # sqrt(pos_rate * (1 - pos_rate) / results$totalwells)

          #upper_CI_pos_rate <- pos_rate + 1.96 *
          #sqrt(pos_rate * (1 - pos_rate) / results$totalwells)

          lambda <- - log(1 - pos_rate)

          lower_CI_lambda <- - log(1 - lower_CI_pos_rate)

          upper_CI_lambda <- - log(1 - upper_CI_pos_rate)

          concentration <- lambda / partition.volume

          lower_CI_concentration <- lower_CI_lambda / partition.volume

          upper_CI_concentration <- upper_CI_lambda / partition.volume

          conc.dil <- concentration / data$samples[[match(
            input$sample, names(data$samples))]]$`raw results`$Dilution[1]

          lower_CI_conc.dil <- lower_CI_concentration / data$samples[[match(
            input$sample, names(data$samples))]]$`raw results`$Dilution[1]

          upper_CI_conc.dil <- upper_CI_concentration / data$samples[[match(
            input$sample, names(data$samples))]]$`raw results`$Dilution[1]

          spread <- cbind(abs(1 - (lower_CI_lambda)),
                          abs(1 - (upper_CI_lambda)))

          precision <- apply(spread, 1, max) * 100

          if (is.null(values$results)) {
            values$results <- results()
          }

          values$results[nrow(values$results) + 1:nrow(results()), ] <-
            cbind.data.frame(
              results, lambda, lower_CI_lambda, upper_CI_lambda, concentration,
              lower_CI_concentration, upper_CI_concentration, conc.dil,
              lower_CI_conc.dil, upper_CI_conc.dil, precision,
              data$samples[[match(
                input$sample, names(data$samples))]]$`raw results`$Dilution)

          names(values$results) <- c(
            "Target", "Positive reactions", "Total reactions", "lambda",
            "Lower CI lambda", "Upper CI lambda", "Copies/ul",
            "Lower CI copies/ul", "Upper CI copies/ul",
            "Copies/ul at sample dilution",
            "Lower CI copies/ul at sample dilution",
            "Upper CI copies/ul at sample dilution", "Precision %", "Dilution")
        }
      })
    })

    #If clusters have changed activate undo action button
    observe({
      shinyjs::toggleState("undo", in.clus$count > 0)
    })

    #Set action for undo action button
    observeEvent(input$undo, {
      in.clus$count <- in.clus$count - 1
      values$new <- dataplot()[, 1:ncol(dataplot()) - 1]
      values$results <- values$results[1:(nrow(values$results) -
                                            nrow(results())), ]
    })

    #Set action for reset action button
    observeEvent(input$reset, {
      values$new <- NULL
      values$results <- NULL
      in.clus$count <- 0
      manual.mod$clus[[match(input$sample, names(data$samples))]] <-
        cluster.data()
      manual.mod$res[[match(input$sample, names(data$samples))]] <- results()
    })

    #Check if manual correction has been done.
    #If so show new results, if not show last saved results
    output$results <- renderTable({
      if (!is.null(values$results)) {
        utils::tail(values$results, nrow(results()))
      } else {
        manual.mod$res[[match(input$sample, names(data$samples))]]
      }
    }, digits = 3)

    observeEvent(input$export, {
      withProgress(message = "Exporting data", value = 0, {

        incProgress(0.5, detail = "Collecting data")

        #If manual correction has been done save the configuration of
        #current sample
        if (!is.null(values$new)) {
          manual.mod$clus[[manual.mod$samples[length(manual.mod$samples)]]] <-
            values$new
          manual.mod$res[[manual.mod$samples[length(manual.mod$samples)]]] <-
            utils::tail(values$results, nrow(results()))
        }


        c.data <- lapply(seq_along(data$referenceDB), function(x) {

          paste(names(
            data$referenceDB)[x], data$referenceDB[[x]]$quality,
            data$referenceDB[[x]]$dbscan$eps,
            data$referenceDB[[x]]$dbscan$minPts, sep = ",")
        })

        c.name <- paste("Reference", "Quality", "DBSCAN eps",
                        "DBSCAN minPts", sep = ",")

        all.ref <- append(c.name, c.data)

        exp.ref <- paste(unlist(all.ref), sep = "\n")

        sample.name <- stringr::str_extract(names(data$samples), ".+?(?=_)")

        exp.data <- lapply(seq_along(data$samples), function(x) {

          repl.index <- which(sample.name[x] == sample.name)

          n.repl <- sum(sample.name[x] == sample.name, na.rm = TRUE)
          n.target <- nrow(data$samples[[x]]$`raw results`)
          tot.rows <- n.repl * n.target


          all.col <- cbind.data.frame(
            rep(names(data$samples)[x], n.target),
            manual.mod$res[[x]],
            "Quality" = rep(data$samples[[x]]$quality, n.target),
            "Reference" = rep(data$samples[[x]]$reference, n.target)
          )

          raw.data <- do.call(paste, c(lapply(all.col, "["), sep = ","))

          tb1 <- ggpubr::ggtexttable(
            cbind(all.col[, 2:4],
                  round(all.col[, c(5:7)], 5),
                  round(all.col[, c(8:10)], 3)),
            rows = NULL,
            cols = c("Target", "Positive \nreactions", "Total \nreactions",
                     "lambda", "Lower CI \nlambda", "Upper CI \nlambda",
                     "Copies/ul", "Lower CI \ncopies/ul",
                     "Upper CI \ncopies/ul"),
            theme = ggpubr::ttheme(base_size = 10))

          tb1 <- ggpubr::ggtexttable(
            cbind(all.col[, 2:4],
                  round(all.col[, 5], 5),
                  paste0(round(all.col[, 8], 2),
                         " (", round(all.col[, 9], 2), " - ",
                         round(all.col[, 10], 2), ")"
                  ),
                  paste0(round(all.col[, 11], 2),
                         " (", round(all.col[, 12], 2), " - ",
                         round(all.col[, 13], 2), ")"
                  ),
                  round(all.col[, 14], 2)),
            rows = NULL,
            cols = c("Target", "Positive\nreactions", "Total\nreactions",
                     "lambda", "Copies/ul\n(95% CI)",
                     "Copies/ul at sample\ndilution (95% CI)", "Precision\n%"),
            theme = ggpubr::ttheme(base_size = 10))

          tx <- paste(
            "Quality threshold:", all.col$Quality,
            "\nDilution:", all.col$Dilution,
            "\nReference:",
            paste0(
              data$samples[[x]]$reference, " (quality: ",
              data$referenceDB[[match(data$samples[[x]]$reference,
                                      names(data$referenceDB))]]$quality,
              ", ", "eps: ",
              data$referenceDB[[match(
                data$samples[[x]]$reference,
                names(data$referenceDB))]]$dbscan$eps,
              ", ", "minPts: ",
              data$referenceDB[[match(
                data$samples[[x]]$reference,
                names(data$referenceDB))]]$dbscan$minPts, ")"
            )
          )

          text.p <- ggpubr::ggparagraph(text = tx, color = "black", size = 12)

          cluscolors <- c(
            "gray70", "#004949", "#ff6db6", "#009292", "#ffb6db", "#490092",
            "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", "#920000", "#924900",
            "#db6d00", "#24ff24", "#ffff6d",
            "#000000")[1:nrow(data$samples[[x]]$centers)]
          names(cluscolors) <- rownames(data$samples[[x]]$centers)

          graph.p <- ggplot(manual.mod$clus[[x]], aes_string(x = "Vic",
                                                             y = "Fam")) +

            geom_point(
              size = 0.5,
              aes(color = as.factor(
                manual.mod$clus[[x]][,ncol(manual.mod$clus[[x]])]))) +

            labs(color = "Clusters", title = names(data$samples)[x]) +

            scale_colour_manual(values = cluscolors, drop = FALSE) +

            theme(
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(fill = NA, color = "black"),
              plot.title = element_text(size = 16),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              axis.text = element_text(size = 10),
              legend.key = element_rect(fill = "white"),
              legend.text = element_text(size = 10),
              legend.title = element_text(face = "bold", size = 12)
            ) +

            guides(colour = guide_legend(override.aes = list(size = 2)))

          if (isTRUE(save.plot)) {
            ggsave(paste0(names(data$samples)[x], "_manual.", format), graph.p,
                   dpi = dpi, width = 9, height = 9)
          }

          tot <- ggpubr::ggarrange(graph.p, text.p, tb1, ncol = 1,
                                   nrow = 3, heights = c(7, 2, 4))


          if (n.repl == 1) {

            exp.repl <- cbind.data.frame(
              rep(names(data$samples)[x], tot.rows),
              manual.mod$res[[x]][, c(1, 10:13)],
              "No of replicates" = rep(n.repl, n.target)
            )
          } else {

            exp.repl <- lapply(seq(n.target), function(y) {

              sing.target <- lapply(repl.index, function(z) {
                if (any(is.na(manual.mod$res[[z]][y, c(4:6, 14)]))) {
                  NULL
                } else {
                  manual.mod$res[[z]][y, c(2, 3)]
                }
              })
              sing.target <- do.call(rbind.data.frame, sing.target)

              pos_tot <- sum(sing.target[, 1])

              total_tot <- sum(sing.target[, 2])

              pos_rate <- pos_tot / total_tot

              if (pos_tot == 0) {
                pois <- exactci::poisson.exact(pos_tot, total_tot,
                                      alternative = "two.sided",
                                      tsmethod = "central", midp = TRUE)

                lower_CI_pos_rate <- pois$conf.int[1]

                upper_CI_pos_rate <- pois$conf.int[2]

              } else if (pos_tot < 100) {
                pois <- exactci::poisson.exact(pos_tot, total_tot,
                                      alternative = "two.sided",
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

              spread <- cbind(abs(1 - (lower_CI_lambda)),
                              abs(1 - (upper_CI_lambda)))

              precision <- apply(spread, 1, max) * 100

              cbind.data.frame(
                "Sample" = sample.name[x],
                "Target" = data$samples[[x]]$`raw results`[y, 2],
                "Copies/ul" = concentration,
                "Lower CI copies/ul" = lower_CI_concentration,
                "Upper CI copies/ul" = upper_CI_concentration,
                "Precision %" = precision,
                "No of replicates" = n.repl)
            })

            exp.repl <- do.call(rbind.data.frame, exp.repl)

          }

          exp.repl <- do.call(paste, c(lapply(exp.repl, "["), sep = ","))

          list("data_results" = raw.data, "results" = exp.repl, "report" = tot)
        })
        names(exp.data) <- sample.name

        raw.name <- paste(
          "Sample", "Target", "Positive reactions", "Total reactions",
          "lambda", "Lower CI lambda", "Upper CI lambda", "Copies/ul",
          "Lower CI copies/ul", "Upper CI copies/ul",
          "Copies/ul at sample dilution",
          "Lower CI copies/ul at sample dilution",
          "Upper CI copies/ul at sample dilution",
          "Precision %", "Dilution", "Quality", "Reference", sep = ",")

        all.raw <- append(
          raw.name,
          lapply(exp.data, function(x) {
            x$data_results
          })
        )

        exp.raw <- paste(unlist(all.raw), sep = "\n")

        repl.name <- paste(
          "Sample", "Target", "Copies/ul", "Lower CI copies/ul",
          "Upper CI copies/ul", "Precision %", "No of replicates", sep = ",")

        all.repl <- append(
          repl.name,
          lapply(match(unique(names(exp.data)), names(exp.data)), function(x) {
            exp.data[[x]]$results
          })
        )

        exp.repl <- paste(unlist(all.repl), sep = "\n")

        all.exp <- rlist::list.append("Reference samples", exp.ref, "",
                                      "Results", exp.raw, "",
                                      "Replicates results", exp.repl)

        all.exp <- strsplit(all.exp, ",")

        all.exp <- sapply(all.exp, function(x) {
          if (length(x) > 10) {
            x[5:16] <- gsub("\\.", ",", x[5:16])
          } else if (length(x) > 6) {
            x[3:6] <- gsub("\\.", ",", x[3:6])
          }
          paste(x, collapse = "\t")
        })

        incProgress(1, detail = "Writing files")

        utils::write.table(all.exp, paste0(filename, ".csv"),
                           row.names = FALSE, col.names = FALSE, sep = "\t",
                           na = "", quote = TRUE)

        report <- lapply(exp.data, function(x) {
          x$report
        })

        ggpubr::ggexport(plotlist = report,
                         filename = paste0(filename, ".pdf"), width = 8,
                         height = 11, res = 300)
      })
    })
  }

  shinyApp(ui = ui, server = server)
}
