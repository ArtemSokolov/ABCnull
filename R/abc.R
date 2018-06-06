## Training of models and prediction on test data
##
## by Artem Sokolov

#' Model training
#'
#' Trains a GBM model on the provided dataset matrix
#'
#' @param XY One of the datasets in this packages, e.g., as loaded by data(MACCSbinary)
#' @return A GBM model constructed on the training examples (i.e., where Label is not NA)
#' @importFrom magrittr %>%
#' @export
ABCtrain <- function( XY )
{
    ## Use a fixed random seed to allow for reproducibility
    set.seed( 100 )

    ## Set up a cross-validation schema
    fc <- caret::trainControl( method="repeatedcv", repeats = 5 )

    ## Perform cross-validation on the training data, which is defined by non-NA Labels
    X1 <- dplyr::select( XY, -Drug, -pubchem_id ) %>% dplyr::filter( !is.na(Label) )
    caret::train( Label ~ ., data=X1, method="gbm", trControl=fc, verbose=FALSE )
}

#' Predictions on test data
#'
#' Applies a previously-trained ABC model to new data
#'
#' @param ABCmodel A model previously trained by ABCtrain()
#' @param XY A dataset where test data is taken to be those samples where Label is NA
#' @return Test data with predictions augmented as a new ABCpred column
#' @importFrom magrittr %>%
#' @export
ABCpredict <- function( ABCmodel, XY )
{
    ## Identify test data
    X1 <- dplyr::filter( XY, is.na(Label) )

    ## Compute predictions
    Ypred <- predict( ABCmodel, X1, type="prob" )[["Sensitive"]]

    ## Augment test data with predictions
    X1 %>% dplyr::mutate( ABCpred = Ypred ) %>% dplyr::select( -Label ) %>%
        dplyr::select( ABCpred, dplyr::everything() )
}
