\name{computeMaxReactionCorrelations}
\alias{computeMaxReactionCorrelations}
\title{Computes a correlation edge weight between two reactions.}
\description{Find the maximum correlation coefficient between all combination of the gene lists between two connected reactions}

\usage{
edge.weights <- computeMaxReactionCorrelations(reaction.network,microarray,y=NULL)
}
\arguments{
\item{reaction.network}{The result from makeReactionNetwork.}
\item{microarray}{The microarray data.}
\item{y}{The experimental conditions for the microarray experiment. (not essential).}
}
\value{
A data.frame where each column corresponds to to the maximum correlation between the gene lists of each connect reaction in the reaction.network edge list.
}
\details{
If the same gene occurs in two connected reactions it is not considered.

If the correlation cannot be computed because of NA's, some other reason it returns a correlation of -1.

NA's are skipped.
}