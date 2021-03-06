% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blandaltmanplots.R
\name{baplot}
\alias{baplot}
\title{Bland-Altman plots}
\usage{
baplot(x, y, regline = TRUE, la = c("log", "lin", "both"), main = "",
  xlab = "Average", ylab = "Difference")
}
\arguments{
\item{x}{a vector or matrix of the observations of one measurment type}

\item{y}{a vector or matrix of the observations of another measurment type}

\item{regline}{plot regression line default TRUE}

\item{la}{limits of agreement log, lin or both}

\item{main}{title of plot default ""}

\item{xlab}{x-axis title defaults to "Average"}

\item{ylab}{y-axis title defaults to "Difference"}
}
\value{
plot
}
\description{
Create Bland-Altman plots
}
\details{
Create Bland-Altman plots
Optionally with log-transformed limits of agreement on original scale
according to Euser et al. Clinical Chemistry 48, No. 5, 2002
}
\examples{
x <- runif(500, 0, 100)
y <- x + runif(500, -5, 5)
baplot(x, y)
plot(x, y)
}
\author{
mvaniterson
}
