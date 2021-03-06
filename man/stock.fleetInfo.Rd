% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions_FLFleets.R
\name{stock.fleetInfo}
\alias{stock.fleetInfo}
\title{stock.fleetInfo}
\usage{
stock.fleetInfo(fleets)
}
\arguments{
\item{fleets}{is an objetc of class FLFleetsExt.}
}
\value{
Return a matrix with rownames equal to the stocks names and colnames equal to names of fleet and metier.
If element (i,j) is equal to 0, the stock (i) is not caught by fleet/metier (j).
}
\description{
Indicates which stocks are caugth by each fleet-metier combination.
}
\examples{
\dontrun{
data(multi)
stock.fleetInfo(fl1)
}
}
