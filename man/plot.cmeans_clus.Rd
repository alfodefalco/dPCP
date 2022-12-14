% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3.clustering_part2.R
\name{cmeans_clus}
\alias{cmeans_clus}
\alias{plot.cmeans_clus}
\title{Cluster analysis with fuzzy c-means algorithm}
\usage{
cmeans_clus(centers.data)

\method{plot}{cmeans_clus}(x, ..., sample = "all", color.blind = FALSE)
}
\arguments{
\item{centers.data}{an object of class \code{centers_data}, inherited
from \code{\link{centers_data}}.}

\item{x}{an object of class \code{cmeans_clus}}

\item{...}{Arguments to be passed to methods}

\item{sample}{'all' to show all samples, or a numeric vector indicating
the row number of samples in the sample table.}

\item{color.blind}{logical. If TRUE colors optimized for colorblind readers
are used.}
}
\value{
An object of class \code{cmeans_clus} containing a sublist for
  each sample. Each sublist has the following components:
  \item{quality}{quality threshold used in \code{\link{read_sample}}.}
  \item{reference}{reference ID.}
  \item{centers}{a data frame with the centroids coordinates.}
  \item{data}{a data frame with the fluorescence intensities and clusters
  name.}
  \item{membership}{a matrix with the membership values of the data elements
  to the clusters. See also \code{\link[e1071]{cmeans}}}
}
\description{
This function carries out the c-means cluster analysis, using the centroids
position as initial values for cluster centers.
}
\examples{
\donttest{
library(dPCP)

#Find path of sample table and location of reference and input files
sampleTable <- system.file("extdata", "Template_sampleTable.csv",
                     package = "dPCP")

fileLoc <- system.file("extdata",package = "dPCP")

#Read sample table file
sample.table <- read_sampleTable(sampleTable, system = "bio-rad",
                                 file.location = fileLoc)

#Read reference files
ref <- read_reference(sample.table, system = "bio-rad",
                      file.location = fileLoc)

#Read samples files
samp <- read_sample(sample.table, system = "bio-rad",
                    file.location = fileLoc)

#Reference DBSCAN clustering
dbref <- reference_dbscan(ref, sample.table, save.template = FALSE)

#Predict position of clusters centroid from reference DBSCAN results
cent <- centers_data(samp, sample.table,dbref)

#Fuzzy c-means clustering
cmclus <- cmeans_clus(cent)

plot(cmclus, sample = "all")
}
}
