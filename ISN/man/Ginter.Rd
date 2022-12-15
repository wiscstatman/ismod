\name{Ginter}
\alias{Ginter}
\title{Generate interactions from data}
\description{
  Take the data and generate single, pair, triple and quadruple
  interactions experienced by more than k1, k2, k3 and k4 tumors,
  respectively. 
}
\usage{
Ginter(data,k1,k2,k3,k4)
}
\arguments{
  \item{data}{ cCGH data}
  \item{k1}{threshold for single}
  \item{k2}{threshold for pair}
  \item{k3}{threshold for triple}
  \item{k4}{threshold for quadruple}
}
\value{
  interactions
}
