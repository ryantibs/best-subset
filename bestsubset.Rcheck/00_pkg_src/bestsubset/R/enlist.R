enlist <- function (...) 
{
    result <- list(...)
    if ((nargs() == 1) & is.character(n <- result[[1]])) {
        result <- as.list(seq(n))
        names(result) <- n
        for (i in n) result[[i]] <- get(i)
    }
    else {
        n <- sys.call()
        n <- as.character(n)[-1]
        if (!is.null(n2 <- names(result))) {
            which <- n2 != ""
            n[which] <- n2[which]
        }
        names(result) <- n
    }
    result
}
