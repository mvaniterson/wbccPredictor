.baplot <- function(x, y, regline, la, main, xlab, ylab){

    ##Log-transformed limits of agreement on original scale
    ## according to Euser et al. Clinical Chemistry 48, No. 5, 2002
    lafun <- function(x) 2*x*(10^(1.96*sd(a)) - 1)/(10^(1.96*sd(a)) + 1)

    b <- (x+y)/2
    a <- x - y

    ylim <- range(a, sd(a)*c(-1.96, 1.96))
    if(la == "log") {
        ylim <- range(a, -1*lafun(b), lafun(b))
    }
    else if(la == "both"){
        ylim <- range(ylim, a, -1*lafun(b), lafun(b))
    }
    ylim <- c(max(c(-100, ylim[1])), min(c(100, ylim[2]))) ##truncate if percentage

    ylim <- c(-1*max(abs(ylim)), max(abs(ylim))) ##make the ylim symmetric around zero

    plot(b, a, main=main, xlab=xlab, ylab=ylab,
         ylim=ylim, bty="n", yaxt="n")
    ticks <- axTicks(2)
    axis(2, at=ticks, labels=ticks, las=2)
    abline(h=0)
    abline(h=mean(a), col="blue")

    if(la == "log") {
        curve(-1*lafun(x), add=TRUE, lty=2, col="blue")
        curve(+1*lafun(x), add=TRUE, lty=2, col="blue")
    }
    else if(la == "both"){
        curve(-1*lafun(x), add=TRUE, lty=2, col="blue")
        curve(+1*lafun(x), add=TRUE, lty=2, col="blue")
        abline(h=mean(a) + sd(a)*c(-1.96, 1.96), lty=2, col="blue")
    }
    else
        abline(h=mean(a) + sd(a)*c(-1.96, 1.96), lty=2, col="blue")

    if(regline) {
        fit <- lm(a~b)
        abline(fit, col="red", lwd=2)
        print(summary(fit))
    }
}

##' Create Bland-Altman plots
##'
##' Create Bland-Altman plots
##' Optionally with log-transformed limits of agreement on original scale
##' according to Euser et al. Clinical Chemistry 48, No. 5, 2002
##' @title Bland-Altman plots
##' @param x a vector or matrix of the observations of one measurment type
##' @param y a vector or matrix of the observations of another measurment type
##' @param regline plot regression line default TRUE
##' @param la limits of agreement log, lin or both
##' @param main title of plot default ""
##' @param xlab x-axis title defaults to "Average"
##' @param ylab y-axis title defaults to "Difference"
##' @return plot
##' @export
##' @author mvaniterson
##' @examples
##' x <- runif(500, 0, 100)
##' y <- x + runif(500, -5, 5)
##' baplot(x, y)
##' plot(x, y)
baplot <- function(x, y, regline=TRUE, la=c("log", "lin", "both"), main="", xlab="Average", ylab="Difference"){

    la <- match.arg(la)

    if(is.null(dim(x)))
        x <- as.matrix(x, ncol=1)

    if(is.null(dim(y)))
        y <- as.matrix(y, ncol=1)

    if(ncol(x) > 4)
        op <- par(mfcol=c(3, 2), mar=c(4, 4, 3, 1))
    else
        op <- par(mfcol=c(3, 1), mar=c(4, 4, 3, 1))
    on.exit(op)

    for(k in 1:ncol(x))
        .baplot(x[,k], y[,k], regline=regline, la=la,
                main=paste0(colnames(x)[k], " (", round(icc(x[,k], y[,k]), 4), ")"),
                xlab=xlab,
                ylab=ylab)

}

##' intra-class correlation
##'
##' intra-class correlation
##' @title intra-class correlation
##' @param x numeric vector
##' @param y numeric vector should be same length as x
##' @return intraclass correlation
##' @export
##' @author mvaniterson
##' @examples
##' x <- runif(500, 0, 100)
##' y <- x + runif(500, -5, 5)
##' plot(x, y)
##' icc(x, y)
icc <- function(x, y){

    ##TODO: can be negative???

    if(length(x) != length(y))
        stop("x and y should contain same number of samples!")
    n <- length(x)
    sa2 <- var(x+y)
    sd2 <- var(x-y)
    dbar <- mean(x-y)
    (sa2 - sd2)/(sa2 + sd2 + 2*(n*dbar^2-sd2)/n)
}
