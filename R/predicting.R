##' prediction cell percentages based on trained model
##'
##' prediction cell percentages based on trained model
##' @title prediction
##' @param predictor results from train
##' @param data matrix of 450k beta values or RNAseq counts
##' @param covariates matrix of covariates to correct for in the model e.g. Age and Gender
##' @param transformation this will depend on the model e.g. 10^x -1 if model is log10(cellPercentages+1)~. Default is no transformation.
##' @param ncomp optimal number of components that will be used
##' @param impute TRUE
##' @param ... not used yet
##' @return predicted cell percentage
##' @author mvaniterson
##' @export
##' @importFrom stats predict
prediction <- function(predictor, data, covariates,
                       transformation=function(x) x, 
                       ncomp=NULL,
                       impute=TRUE,
                       ...) {

    ##some checking
    if(is.null(data) | is.null(covariates))
        stop("Data and covariates must be given!")
    if(ncol(data) != nrow(covariates))
        stop("The number of columns of data must match the number of rows of covariates!")
    if(!all.equal(colnames(data), rownames(covariates)))
        stop("Columns names of data must match rownames of covariates!")

    ##matching data covariates
    mId <- match(rownames(covariates), colnames(data))
    if(any(is.na(mId)))
        stop("colnames data do not match rownames covariates!")
    covariates <- covariates[mId,]

    if(class(predictor) == "mvr")
        names <- dimnames(coefficients(predictor))[[1]]
    else if(class(predictor) == "matrix")
        names <- rownames(predictor)
    else
        stop(paste("Unknown class", class(predictor), "of predictor!"))
    
    covarNames <- gsub("covariates", "", grep("covariates", names, value=TRUE))
    dataNames <- gsub("data", "", grep("data", names, value=TRUE))

    ##matching covariates and model
    mId <- match(covarNames, colnames(covariates))
    if(any(is.na(mId)))
        stop("colnames covariates do not match covariates of the predictor!")
    covariates <- covariates[,mId]

    if(any(is.na(covariates)) & !impute) {
        stop("covariates contains NA's not allowed!")
    } else if(any(is.na(covariates)) & impute) {
        print(paste("There are", sum(is.na(covariates)), "NA's in the covariate matrix.",
                    "These will be median imputed."))
        covariates <- apply(covariates, 2, function(x) {
            x[is.na(x)] = median(x, na.rm=TRUE)
            x})
    }

    ##matching data and model
    mId <- match(dataNames, rownames(data))
    if(any(is.na(mId)))
        warning("rownames data do not match data of the predictor!")
    data <- data[mId,]
    if(any(is.na(data)) & !impute) {
        stop("data contains NA's not allowed!")
    } else if(any(is.na(data)) & impute) {
        print(paste("There are", sum(is.na(data)), "NA's in the data matrix.",
                    "These will be median imputed."))
        nas <- apply(data, 1, function(x) any(is.na(x)))
        data[nas,] <- apply(data[nas, ], 1, function(x) median(x, na.rm=TRUE)) ##impute over rows
        data[is.na(data)] <- median(data, na.rm=TRUE) ##maybe some rows are completely NA
    }

    ##predict
    if(class(predictor) == "mvr") {
        predicted <- predict(predictor, newdata = list(covariates=covariates, data=t(data)), ncomp=ncomp, ...)
        predicted <- predicted[,,1]
    }
    else if(class(predictor) == "matrix") {
        predicted <-  cbind(1, covariates, t(data)) %*% predictor
    }
        
   invisible(transformation(predicted))
}
