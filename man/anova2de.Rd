% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anova2de.R
\name{anova2de}
\alias{anova2de}
\title{Use ANOVA to filter low epxressed genes}
\usage{
anova2de(data.exp, ncol.idx, pval.cut = 0.05, model = "expression~time")
}
\arguments{
\item{data.exp}{gene expression data matrix with variables in columns and samples in rows. The first columns are the indeces for comparison, such as
time points and conditions for comparision.}

\item{ncol.idx}{the number of index columns.}

\item{pval.cut}{a numeric value for significance cut-off. The default is \code{pval.cut=0.05}.}

\item{model}{the model used for comparison. The default is \code{model="expression~time"}.}
}
\value{
\code{anova2de} returns a list includes the ANOVA p-value results, the names of differentially expressed genes and the time-series data of differential expressed genes.
}
\description{
Analysis of variance is used to filter the low expressed genes. More details can be seen in \code{\link{aov}}.
}
\seealso{
\code{\link{aov}}, \code{\link{RLowPC}}
}

