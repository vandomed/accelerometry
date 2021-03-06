% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inverse_rle2.R
\name{inverse_rle2}
\alias{inverse_rle2}
\title{Inverse Run Length Encoding (Alternate Implementation)}
\usage{
inverse_rle2(x)
}
\arguments{
\item{x}{Object returned by \code{\link{rle2}}.}
}
\value{
Integer or numeric vector.
}
\description{
Re-constructs vector compressed by \code{\link{rle2}}.
}
\examples{
# Create dummie vector x
x <- c(0, 0, 0, -1, -1, 10, 10, 4, 6, 6)

# Summarize x using rle2
x.summarized <- rle2(x)

# Reconstruct x
x.reconstructed <- inverse_rle2(x.summarized)


}
