\name{readSBML}
\alias{readSBML}
\title{Reads an SBML}
\description{Reads some elements of an SBML file using libsbml}

\usage{
sbmlinfo <- readSBML(filename) 
}
\arguments{
\item{filename}{The SBML file name.}
}
\value{
A list with the following elements:
\item{edge.list}{The compound network with edges labelled by reaction.}
\item{reaction.info}{A list named by reaction name and each element corresponds to the list of genes required for that reaction to proceed.}
\item{species.info}{A data.frame listing all compounds and their properties within the SBML file.}
\item{sbml.info}{The complete SBML info.}
}
\details{
If it can't find the file it will crash R.
}