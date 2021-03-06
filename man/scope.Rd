\name{scope}
\alias{scope}
\title{Extracts the scope from a specified network}
\description{Searches all paths between the specified start and end vertices, and if a significant path is found it returns it as part of the scope.}

\usage{
cscopes <- scope(geneNetwork,edge.weights,sampledpaths,alpha = 0.01,echo = 0)
}
\arguments{
\item{geneNetwork}{The result from defineReactionNeighborhood.}
\item{edge.weights}{The result from computeMaxReactionCorrelations.}
\item{sampledpaths}{The result from samplePaths.}
\item{alpha}{The significance value used to test each extracted path.  This value is bonferoni corrected by the number of compounds within the neighborhood - default is 0.01.}
\item{echo}{Print out the scope identification process}
}
\value{
A scope object which is a list, where the first element contains the paths extracted and the second element contains a list of all paths found.
}
\details{
It doesn't search for the best path between the start and terminal vertices.  As there could be many paths which lead to the same terminal vertex, searching through all of them is time comsuming.  We just stop when the first significant path is found.
}
