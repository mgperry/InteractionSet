\name{Convert classes}

\alias{inflate}
\alias{inflate,GInteractions-method}
\alias{inflate,InteractionSet-method}
\alias{deflate}
\alias{deflate,ContactMatrix-method}

\title{Convert between classes}

\description{
Inflate a GInteractions or InteractionSet into a ContactMatrix, or deflate a ContactMatrix to an InteractionSet.
}

\usage{
\S4method{inflate}{GInteractions}(x, rows, columns, fill, ...)

\S4method{inflate}{InteractionSet}(x, rows, columns, assay=1L, sample, ...)

\S4method{deflate}{ContactMatrix}(x, unique=TRUE, ...)
}

\arguments{
\item{x}{
    A GInteractions or InteractionSet object for \code{inflate}, or a ContactMatrix object for \code{deflate}.
}

\item{rows, columns}{
    An integer vector, character vector or GRanges object indicating the regions of interest to be used as the rows or columns of the ContactMatrix.

    If it is an integer vector, it is assumed to refer to intervals in the \code{regions} slot of \code{x}.
    Values of the vector need not be sorted or unique, but must lie within \code{[1, regions(x)]}.

    If it is a character vector, it is assumed to contain the names of the reference sequences of interest (i.e., chromosome names).

    If it is a GRanges object, \code{\link[GenomicRanges]{overlapsAny}} will be called to identify the overlapping intervals of the \code{regions} slot.
}

\item{assay}{
    A numeric scalar indicating the assay of the InteractionSet object, from which values are extracted to fill the ContactMatrix.
}

\item{sample}{
    A numeric vector indicating the samples of the assay from which values are to be extracted, to fill the ContactMatrix.
}

\item{fill}{
    A vector of length equal to \code{nrow(x)}, containing values with which to fill the ContactMatrix.
    Alternatively, a matrix with number of rows equal to \code{nrow(x)}.
}

\item{unique}{
    A  logical scalar indicating whether duplicated interactions should be removed from \code{x} prior to deflation.
}

\item{...}{
    For \code{inflate}, additional arguments to pass to \code{\link[GenomicRanges]{overlapsAny}} when \code{rows} or \code{columns} is a GRanges.

    For \code{deflate}, additional arguments to pass to the InteractionSet constructor.
}
}

\details{
The \code{inflate} method will return a ContactMatrix containing one or more matrices, 
    where the rows and columns correspond to specified regions of interest in \code{rows} and \code{columns}.
For the GInteractions method with a vector \code{fill}, the ContactMatrix will contain only one matrix.
Each entry of \code{fill} corresponds to a pairwise interaction in \code{x} and is used to set the matrix value at the matching row/column.
Some entries of the matrix will correspond to interactions that are not present in the InteractionSet - these are filled with \code{NA} values.
If duplicated interactions are present in \code{x}, only the first interaction will be used.

Alternatively, if a \code{fill} matrix is supplied, the number of matrices in the returned ContactMatrix is defined from the number of columns in \code{fill}.
Each row contains values for a pairwise interaction, corresponding to an entry in \code{x}.
Each matrix is filled according to the values in the corresponding column of \code{fill}.
 
For the InteractionSet \code{inflate} method, the number of matrices in the ContactMatrix is equal to the number of samples.
Values in each matrix are filled in based on the assay values for the corresponding sample, by default.
The assay to be used and the samples to be inflated can be controlled by altering the \code{assay} and \code{sample} arguments.
By default, all samples are inflated from the first assay.

The \code{deflate} method will return an InteractionSet where each non-\code{NA} entry in the ContactMatrix is converted into a pairwise interaction.
The values of all matrix entries are stored as a one-sample assay.
Duplicate interactions are not included by default, but can be returned by setting \code{unique=FALSE}.
Additional arguments can be used to specify the \code{colData} and \code{metadata}, which are stored in the ContactMatrix itself.
}

\examples{
example(InteractionSet, echo=FALSE)

inflate(iset, 1:10, 1:10)
inflate(iset, 1:10, 1:5+10)
inflate(iset, "chrA", 1:5+10)
inflate(iset, "chrA", "chrB")
inflate(iset, "chrA", GRanges("chrB", IRanges(1, 10)))

y <- inflate(iset, 1:10, 1:10)
iset2 <- deflate(y)
iset2
assay(iset2)

# Testing with different fillings:
y <- inflate(iset, 1:10, 1:10, sample=2)
iset2 <- deflate(y)
assay(iset2)

y <- inflate(interactions(iset), 1:10, 1:10, fill=rowSums(assay(iset)))
iset2 <- deflate(y)
assay(iset2)

# Effect of 'unique'
y <- inflate(iset, c(8, 1:10), 1:10)
deflate(y)
deflate(y, unique=FALSE)
}

\author{
Aaron Lun
}

\seealso{
\code{\link{InteractionSet-class}},
\code{\link{GInteractions-class}},
\code{\link{ContactMatrix-class}}
}
