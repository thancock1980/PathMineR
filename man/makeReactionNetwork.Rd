\name{makeReactionNetwork}
\alias{makeReactionNetwork}
\title{Converts a compound-network into a reaction network}
\description{Takes the input from readSBML and returns a reaction network.}

\usage{
rnet <- makeReactionNetwork(kobj,removeSmallCompounds = NULL) 
}
\arguments{
\item{kobj}{The result from readSBML.}
\item{removeSmallCompounds}{A list of compounds to be removed before the reaction network is built.}
}
\value{
A list with the following elements:
\item{nodes}{The names of the network nodes.}
\item{edges}{The reaction network edge list. Each index corresponds to a node within the node list}
\item{reaction.info}{A list named by reaction name and each element corresponds to the list of genes required for that reaction to proceed.}
\item{sbml.info}{The complete SBML info.}
}
\details{
May take a little while.
}