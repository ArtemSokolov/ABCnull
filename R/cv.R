## Scores training data via cross-validation, for ROC curve generation
##
## by Artem Sokolov

#' Cross-validation
#'
#' Evaluates a single random cross-validation split for a set of models
#'
#' @param XY One of the datasets in this packages, e.g., as loaded by data(MACCSbinary)
#' @return A data frame containing cross-validation predictions
#' @importFrom magrittr %>%
#' @export
ABCcv <- function( XY )
{
    ## Set up a cross-validation schema
    fc <- caret::trainControl( method="cv", savePred = "all", classProbs = TRUE )

    ## Remove Drug and pubchem_id columns from consideration
    ## Isolate the training data, which is defined by non-NA labels
    X1 <- dplyr::select( XY, -Drug, -pubchem_id ) %>% dplyr::filter( !is.na(Label) )

    ## Each method has a slightly different flag for suppressing extraneous output
    ## Encapsulate the distinction into a functional module
    mytrain <- function( m )
    {
        cat( "Training", m, "...\n" )
        
        if( m == "gbm" )
            cv <- caret::train( Label ~ ., data=X1, method=m, trControl=fc, verbose = FALSE )
        else if( m == "nnet" )
            cv <- caret::train( Label ~ ., data=X1, method=m, trControl=fc,
                               trace = FALSE, MaxNWts=2000 )
        else
            cv <- caret::train( Label ~ ., data=X1, method=m, trControl=fc )
        cv$pred
    }

    ## Evaluate a set of models on the same cross-validation split
    RR <- c( "knn", "gbm", "glmnet", "svmLinear", "nnet" ) %>%
        rlang::set_names( c( "k-NN", "GBM", "Log.Reg.", "SVM", "NNet" ) ) %>%
        purrr::map( mytrain ) %>% dplyr::bind_rows( .id = "Method" )

    ## Average performance across parameter grid for each method
    R1 <- RR %>% dplyr::group_by( Method, rowIndex ) %>%
        dplyr::summarize( Pred = mean(Sensitive, na.rm=TRUE) ) %>% dplyr::ungroup()

    ## Match up predictions to the original pubchem_ids and labels
    XY %>% dplyr::mutate( rowIndex = 1:n() ) %>%
        dplyr::select( Label, Drug, pubchem_id, rowIndex ) %>%
            dplyr::inner_join( R1, by="rowIndex" ) %>% dplyr::select( -rowIndex )
}

