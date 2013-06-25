\name{samplePaths}
\alias{samplePaths}
\title{Creates a set of sample path p-values for each length given a weighted network.}
\description{Randomly traverses paths of increasing lengths within a set network to create an empirical pathway distribution for more accurate determination of path significance.}

\usage{
path.samples <- samplePaths(rnet,edge.weights,maxpathlength,numberofsamples,numberofwarmup)
}
\arguments{
\item{rnet}{The result from makeReactionNetwork.}
\item{edge.weights}{The result from computeMaxReactionCorrelations.}
\item{maxpathlength}{The maxmimum path length.}
\item{numberofwarmup}{The number of warm up paths to sample.}
}
\value{
A matrix where each row is a path length and each column is the number of paths sampled.
}
\details{
Can take a bit of time.
}