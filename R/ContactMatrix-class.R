##############################################
# Defines the ContactMatrix class.

setClass("ContactMatrix",
    contains="Vector", 
    slots=list(
        contacts="Assays",
        anchor1="integer",
        anchor2="integer",
        regions="GRanges"
    ),
    prototype=prototype(
        contacts=Assays(),
        anchor1=integer(0),
        anchor2=integer(0),
        regions=GRanges()
    )
)

setValidity("ContactMatrix", function(object) {
    if (is.unsorted(object@regions)) {
        return("'regions' should be sorted")
    }
    msg <- .check_inputs(object@anchor1, object@anchor2, object@regions, same.length=FALSE)
    if (is.character(msg)) { return(msg) }

    if (nrow(object@contacts)!=length(object@anchor1)) { 
        return("'contacts' nrow must be equal to length of 'anchor1'")
    }
    if (ncol(object@contacts)!=length(object@anchor2)) {
        return("'contacts' ncol must be equal to length of 'anchor2'")
    }
    if (length(object@contacts)!=nrow(object@elementMetadata)) { 
        return("'contacts' length should be equal to 'elementMetadata' nrow")
    }
    return(TRUE)
}) 

setMethod("show", signature("ContactMatrix"), function(object) {
    cat("class:", class(object), "\n")
    cat("dim:", dim(object@contacts), "\n")

    expt <- names(metadata(object))
    if (is.null(expt))
        expt <- character(length(metadata(object)))
    scat("metadata(%d): %s\n", expt)

    nms <- names(contacts(object, sample=NA))
    if (is.null(nms))
        nms <- character(length(assays(object, withDimnames=FALSE)))
    scat("assays(%d): %s\n", nms)

    mcolnames <- names(mcols(object))
    fmt <- "metadata column names(%d): %s\n"
    scat(fmt, mcolnames)

    cat(sprintf("regions: %i\n", length(object@regions)))
})

##############################################
# Constructor:

.new_ContactMatrix <- function(contacts, anchor1, anchor2, regions, metadata) {
    anchor1 <- as.integer(anchor1)
    anchor2 <- as.integer(anchor2)
    
    msg <- .check_inputs(anchor1, anchor2, regions, same.length=FALSE)
    if (is.character(msg)) { stop(msg) }
    out <- .resort_regions(anchor1, anchor2, regions, enforce.order=FALSE)

    # Coercing it to be an assays object, if it isn't already.        
    if (!is(contacts, "Assays")) { 
        contacts <- Assays(contacts)
    }
    nlibs <- length(contacts)
    if (length(unique(names(contacts)))!=nlibs) {
        names(contacts) <- paste0("Sample", seq_len(nlibs))
    }
    elementMetadata <- new("DataFrame", nrows=nlibs) 

    new("ContactMatrix", contacts=contacts, anchor1=out$anchor1, anchor2=out$anchor2, 
        regions=out$regions, elementMetadata=elementMetadata, metadata=metadata)
}

setGeneric("ContactMatrix", function(contacts, anchor1, anchor2, ...) { standardGeneric("ContactMatrix") })
setMethod("ContactMatrix", c("ANY", "numeric", "numeric"), 
    function(contacts, anchor1, anchor2, regions, metadata=list()) { 
        .new_ContactMatrix(contacts, anchor1, anchor2, regions, metadata)
    }
)

setMethod("ContactMatrix", c("ANY", "GRanges", "GRanges"), 
    function(contacts, anchor1, anchor2, regions, metadata=list()) { 

        if (missing(regions)) { 
            collated <- .collate_GRanges(anchor1, anchor2)
            regions <- collated$ranges
            anchor1 <- collated$indices[[1]]
            anchor2 <- collated$indices[[2]]
        } else {
            anchor1 <- match(anchor1, regions)
            anchor2 <- match(anchor2, regions)
            if (any(is.na(anchor1)) || any(is.na(anchor2))) {
                 stop("anchor regions missing in specified 'regions'")
            }
        }
        
        .new_ContactMatrix(contacts, anchor1, anchor2, regions, metadata)
    }
)

##############################################
# Matrix dimensions

setMethod("dim", "ContactMatrix", function(x) { 
    dim(x@contacts)
})

setMethod("length", "ContactMatrix", function(x) { 
    length(x@contacts)
})

setGeneric("contacts", function(x, ...) { standardGeneric("contacts") })
setMethod("contacts", "ContactMatrix", function(x, sample=1) {
    if (is.na(sample)) { 
        return(x@contacts) 
    }
    return(x@contacts[[sample]])
}) 

setGeneric("contacts<-", function(x, ..., value) { standardGeneric("contacts<-") });
setReplaceMethod("contacts", "ContactMatrix", function(x, sample=1, ..., value) {
    if (is.na(sample)) { 
        if (!is(value, "Assays")) { 
            value <- Assays(value)
        }
        x@contacts <- value
    } else {
        x@contacts[[sample]][] <- value
    }
    return(x)
}) 

setMethod("$", "ContactMatrix", function(x, name) {
    return(mcols(x)[[name]])          
})

setReplaceMethod("$", "ContactMatrix", function(x, name, value) {
    mcols(x)[[name]] <- value
    return(x)    
})

##############################################
# Subsetting

setMethod("[", c("ContactMatrix", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    if (!missing(i) && !missing(j)) { 
        x@anchor1 <- x@anchor1[i]
        x@anchor2 <- x@anchor2[j]
        x@contacts <- x@contacts[i,j]
    } else if (!missing(i)) { 
        x@anchor1 <- x@anchor1[i]
        x@contacts <- x@contacts[i,]
    } else if (!missing(j)) {
        x@anchor2 <- x@anchor2[j]
        x@contacts <- x@contacts[,j]
    }
    return(x)
}) 

setMethod("[<-", c("ContactMatrix", "ANY", "ANY", "ContactMatrix"), function(x, i, j, ..., value) {
    if (!identical(regions(value), regions(x))) { 
        stop("replacement and original 'regions' must be identical")
    }
    if (!missing(i) && !missing(j)) { 
        x@anchor1[i] <- value@anchor1
        x@anchor2[j] <- value@anchor2
        x@contacts[i,j] <- value@contacts
    } else if (!missing(i)) { 
        x@anchor1[i] <- value@anchor1
        x@contacts[i,] <- value@contacts
    } else if (!missing(j)) { 
        x@anchor2[j] <- value@anchor2
        x@contacts[,j] <- value@contacts
    }
    return(x)
})

setMethod("subset", "ContactMatrix", function(x, i, j) {
    x[i,j]
})

##############################################
# Combining

setMethod("cbind", "ContactMatrix", function(..., deparse.level=1) {
    incoming <- list(...)
    ref <- incoming[[1]]
    for (x in incoming[-1]) {
        if (!identical(regions(ref), regions(x))) { 
            stop("'regions' must be identical for 'cbind'")
        }
        if (!identical(anchors(ref, type="row", id=TRUE),
                       anchors(x, type="row", id=TRUE))) {
            stop("row anchor indices must be identical for 'cbind'")
        }    
    }
    
    ref@contacts <- do.call(cbind, lapply(incoming, contacts, sample=NA))
    ref@anchor2 <- unlist(lapply(incoming, anchors, id=TRUE, type="column"))
    return(ref)
})

setMethod("rbind", "ContactMatrix", function(..., deparse.level=1) {
    incoming <- list(...)
    ref <- incoming[[1]]
    for (x in incoming[-1]) {
        if (!identical(regions(ref), regions(x))) { 
            stop("'regions' must be identical for 'rbind'")
        }
        if (!identical(anchors(ref, type="column", id=TRUE),
                       anchors(x, type="column", id=TRUE))) {
            stop("column anchor indices must be identical for 'rbind'")
        }    
    }
    
    ref@contacts <- do.call(rbind, lapply(incoming, contacts, sample=NA))
    ref@anchor1 <- unlist(lapply(incoming, anchors, id=TRUE, type="row"))
    return(ref)
})

setGeneric("sbind", function(..., deparse.level=1) { standardGeneric("sbind") }, signature="...")
setMethod("sbind", "ContactMatrix", function(..., deparse.level=1) {
    incoming <- list(...)
    ref <- incoming[[1]]
    nsamples <- length(ref)
    collected <- list(mcols(ref))

    for (x in incoming[-1]) {
        if (!identical(regions(ref), regions(x))) { 
            stop("'regions' must be identical for 'lbind'")
        }
        if (!identical(anchors(ref, type="row", id=TRUE),
                       anchors(x, type="row", id=TRUE))) {
            stop("row anchor indices must be identical for 'lbind'")
        }    
        if (!identical(anchors(ref, type="column", id=TRUE),
                       anchors(x, type="column", id=TRUE))) {
            stop("column anchor indices must be identical for 'lbind'")
        }    

        # Going through and adding assays together
        current <- seq_len(length(x))
        for (i in current) {
            ref@contacts[[nsamples+i]] <- x@contacts[[i]]
        }
        names(ref@contacts)[nsamples+current] <- names(x@contacts)
        nsamples <- nsamples + length(x)
        collected <- c(collected, mcols(x))
    }

    ref@elementMetadata <- do.call(rbind, collected)
    return(ref)        
})

setMethod("t", "ContactMatrix", function(x) {
    tmp <- Assays() 
    for (i in seq_along(x@contacts)) { tmp[[i]] <- t(x@contacts[[i]]) }
    x@contacts <- tmp

    tmp <- x@anchor1
    x@anchor1 <- x@anchor2
    x@anchor2 <- tmp
    return(x)
})

##############################################
# Sorting and ordering

setMethod("order", "ContactMatrix", function(..., na.last=TRUE, decreasing=FALSE) {
    incoming <- list(...)
    all.rows <- lapply(incoming, anchors, type="row", id=TRUE)
    all.columns <- lapply(incoming, anchors, type="column", id=TRUE)
    list(row=do.call(order, c(all.rows, na.last=na.last, decreasing=decreasing)),
         column=do.call(order, c(all.columns, na.last=na.last, decreasing=decreasing)))
})

setMethod("sort", "ContactMatrix", function(x, decreasing=FALSE, ...) {
    out <- order(x, decreasing=decreasing)
    x[out$row, out$column]
})

setMethod("duplicated", "ContactMatrix", function(x, incomparables=FALSE, ...) {
    r1 <- duplicated(x@anchor1, incomparables=incomparables, ...)
    r2 <- duplicated(x@anchor2, incomparables=incomparables, ...)
    return(list(row=r1, column=r2))
})

setMethod("unique", "ContactMatrix", function(x, incomparables=FALSE, ...) {
    is.dup <- duplicated(x, incomparables=incomparables, ...)
    return(x[!is.dup$row,!is.dup$column])
})

##############################################
# overlapsAny

setMethod("overlapsAny", c("ContactMatrix", "GRanges"), 
    function(query, subject, maxgap=0L, minoverlap=1L,
        type=c("any", "start", "end", "within", "equal"),
        algorithm=c("nclist", "intervaltree"), ignore.strand=TRUE) {
        a1 <- anchors(query, id=TRUE, type="row")
        a2 <- anchors(query, id=TRUE, type="column")
        
        is.used <- union(a1, a2)
        is.overlapped <- logical(length(regions(query)))
        is.overlapped[is.used] <- overlapsAny(regions(query)[is.used], subject, maxgap=maxgap,
                                        minoverlap=minoverlap, type=type, algorithm=algorithm, 
                                        ignore.strand=ignore.strand)
        return(list(row=is.overlapped[a1], column=is.overlapped[a2]))
})

# Use outer(output$row, output$column, "|" or "&") to get the logical area in the interaction space.
# Not sure it makes a great deal of sense to define 'findOverlaps' here.

##############################################
# End

