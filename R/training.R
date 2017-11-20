##' train predictor
##'
##'
##' @title train
##' @param data matrix with gene expression(counts) or DNA methylation beta-values
##' @param covariates matrix of covariates to correct for in the model e.g. Age and Gender
##' @param cellPercentages matrix of cell percentages for which the predictor will be trained
##' @param model formula e.g. log10(cellPercentages+1)~covariates+data
##' @param ncomp total number of pls components that will be tested
##' @param keep.model logical default FALSE do not return full pls-model
##' @param ... additional parameters for plsr see ?plsr
##' @return prediction model plsr object
##' @author mvaniterson
##' @export
##' @importFrom pls plsr
train <- function(data, covariates, cellPercentages,
                  model=formula(cellPercentages~covariates+data),
                  ncomp = 50, keep.model = FALSE,
                  ...){

    ##some checking
    if(is.null(data) | is.null(covariates) | is.null(cellPercentages))
        stop("Data, covariates and cellPercentage must be given!")
    if(ncol(data) != nrow(covariates) | nrow(covariates) != nrow(cellPercentages))
        stop("The number of columns of data must match the number of rows of covariates and cellPercentages!")
    if(!all.equal(colnames(data), rownames(covariates)) | !all.equal(rownames(covariates), rownames(cellPercentages)))
        stop("Columns names of data must match rownames of covariates and cellPercentages!")

    ##matching
    mId <- match(rownames(covariates), colnames(data))
    if(any(is.na(mId)))
        stop("colnames data do not match rownames covariates!")
    covariates <- covariates[mId,]
    mId <- match(rownames(covariates), colnames(data))
    if(any(is.na(mId)))
        stop("colnames data do not match rownames cellPercentages!")
    cellPercentages <- cellPercentages[mId,]

    if(any(is.na(sum(data))) | any(is.na(sum(covariates))) | any(is.na(sum(cellPercentages))))
        stop("NA's are not allowed in data, covariates and cellPercentages!")

    ##train the model
    predictor <- plsr(model,
                      ncomp=ncomp,
                      data=list(cellPercentages = cellPercentages, covariates=covariates, data=t(data)),
                      ...)

    ##remove data model
    if(!keep.model)
        predictor$model <- NULL
    invisible(predictor)
}
