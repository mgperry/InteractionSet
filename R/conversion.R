# Inflate from InteractionSet to ContactMatrix.

setGeneric("inflate", function(x, ...) { standardGeneric("inflate") })

.make_to_indices <- function(regs, i, ...) {
    nregs <- length(regs)
    if (is.numeric(i)) { 
        i <- as.integer(i)
        if (any(!is.finite(i)) || any(i<=0L) || any(i > nregs)) { 
            stop("indices must be positive integers no greater than 'length(regions(x))'") 
        }
        return(i)
    } else if (is.character(i)) { 
        return(which(seqnames(regs) %in% i))
    } else if (is(i, "GRanges")) {
        return(which(overlapsAny(regs, i, ...)))
    } else {
        stop("invalid value for row/column selection")
    }
}

setMethod("inflate", "GInteractions", function(x, rows, columns, fill, ...) {
    row.chosen <- .make_to_indices(regions(x), rows, ...)
    col.chosen <- .make_to_indices(regions(x), columns, ...)
    if (is.vector(fill)) { 
        if (length(fill)!=nrow(x)) { 
            fill <- rep(fill, length.out=nrow(x))
        }
        fill <- as.matrix(fill)
    } else {
        if (nrow(fill)!=nrow(x)) { stop("nrow of 'fill' and 'x' should be the same") }
    }
     
    # Removing duplicated rows and resorting (we'll put them back in later)
    ro <- order(row.chosen)
    co <- order(col.chosen)
    row.chosen <- row.chosen[ro]
    col.chosen <- col.chosen[co]
    rnd <- !duplicated(row.chosen)
    cnd <- !duplicated(col.chosen)
    row.chosen <- row.chosen[rnd]   
    col.chosen <- col.chosen[cnd]

    # Duplicated interactions can't be handled.
    dx <- duplicated(x)
    if (any(dx)) { 
        warning("duplicated interactions in 'x' are removed")
        x <- x[!dx,]
        fill <- fill[!dx,]
    }

    # Matching.
    a1 <- anchors(x, type="first", id=TRUE)
    a2 <- anchors(x, type="second", id=TRUE)
    ar1 <- match(a1, row.chosen)
    ac1 <- match(a1, col.chosen)
    ar2 <- match(a2, row.chosen)
    ac2 <- match(a2, col.chosen)

    # Setting up the fill indices.
    nR <- length(row.chosen)
    nC <- length(col.chosen)
    relevantA <- !is.na(ar1) & !is.na(ac2)
    relevantB <- !is.na(ar2) & !is.na(ac1)
    indicesA <- (ac2[relevantA] - 1L) * nR + ar1[relevantA]
    indicesB <- (ac1[relevantB] - 1L) * nR + ar2[relevantB]

    # Also setting up the permutation to restore original order.
    original.rows <- cumsum(rnd)
    original.rows[ro] <- original.rows
    original.cols <- cumsum(cnd)
    original.cols[co] <- original.cols

    contacts <- Assays()
    for (lib in seq_len(ncol(fill))) { 
        out.mat <- matrix(NA, nR, nC)
        out.mat[indicesA] <- fill[relevantA,lib] 
        out.mat[indicesB] <- fill[relevantB,lib] 
        out.mat <- out.mat[original.rows,original.cols,drop=FALSE]
        contacts[[lib]] <- out.mat
    }

    return(ContactMatrix(contacts, row.chosen[original.rows], col.chosen[original.cols], 
                         regions(x), metadata=metadata(x)))
})
 
setMethod("inflate", "InteractionSet", function(x, rows, columns, assay=1, sample, ...) {
    fill <- assay(x, assay) 
    sample.data <- colData(x)
    if (!missing(sample)) { 
        fill <- fill[,sample]
        sample.data <- sample.data[sample,]
    }
    final <- inflate(interactions(x), rows, columns, fill=fill, ...)
    mcols(final) <- sample.data
    return(final)
})

setGeneric("deflate", function(x, ...) { standardGeneric("deflate") })

setMethod("deflate", "ContactMatrix", function(x, unique=TRUE, ...) {
    row.index <- rep(anchors(x, type="row", id=TRUE), ncol(x))
    col.index <- rep(anchors(x, type="column", id=TRUE), each=nrow(x))

    nlibs <- length(x)
    output.matrix <- NULL
    for (i in seq_len(nlibs)) {
        all.values <- as.vector(contacts(x, sample=i))

        if (i==1L) { 
            is.valid <- !is.na(all.values)
            output.matrix <- matrix(NA, sum(is.valid), nlibs)
            row.index <- row.index[is.valid]
            col.index <- col.index[is.valid]
        } else if (!identical(is.valid, !is.na(all.values))) { 
            stop("inconsistent 'NA' values between samples")
        }
        output.matrix[,i] <- all.values[is.valid]
    }

    out <- .enforce_order(row.index, col.index)
    dim(all.values) <- c(length(all.values), 1L)
    final <- InteractionSet(output.matrix, GInteractions(out$anchor1, out$anchor2, regions(x)), 
                    colData=mcols(x), metadata=metadata(x), ...)

    if (unique) { final <- unique(final) }
    return(final)
})

