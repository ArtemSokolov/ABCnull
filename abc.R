## Data Analysis of the ABC / GM dataset
##
## by Artem Sokolov

library( tidyverse )
library( caret )

## "Quiet" version of read_csv
readq_csv <- function(...) {read_csv(..., col_types=cols())}

## Uses cross-validation to select the best model on training data
## Scores test data using this best model
pred1 <- function( XY, method = "gbm" )
{
    ## Use a fixed random seed to allow for reproducibility
    set.seed( 100 )
    
    ## Convert Label to a factor to have it treated as a binary classification task
    XY <- XY %>% mutate( Label = factor( Label, levels=c(0,1) ) )

    ## Set up a cross-validation schema
    fc <- trainControl( method="repeatedcv", repeats = 5 )
    
    ## Perform cross-validation on the training data, which is defined by non-NA Labels
    X1 <- select( XY, -Drug, -pubchem_id )
    if( method == "gbm" )
        cv <- train( Label ~ ., data = filter( X1, !is.na(Label) ), method=method, trControl = fc,
                    verbose=FALSE )
    else
        cv <- train( Label ~ ., data = filter( X1, !is.na(Label) ), method=method, trControl = fc )
    
    ## Predict labels on the test data
    y.pred <- predict( cv, filter( X1, is.na(Label) ), type="prob" )[["1"]]
    RR <- XY %>% filter( is.na( Label ) ) %>% select( Drug, pubchem_id ) %>% cbind( GMprob = y.pred )
    RR
}

## All predictions from models trained on final labels
pred.all <- function( method = "gbm" )
{
    ## Evaluate each feature space separately
    ## Combine the results into a single data frame afterwards
    R <- list()
    R[[1]] <- readq_csv("data/MACCSbinary.csv") %>% pred1(method) %>% rename( MACCSbin = GMprob )
    R[[2]] <- readq_csv("data/MACCScount.csv") %>% pred1(method) %>% rename( MACCScnt = GMprob )
    R[[3]] <- readq_csv("data/PChem.csv") %>% pred1(method) %>% rename( PChem = GMprob )
    R[[4]] <- readq_csv("data/Morgan.csv") %>% pred1(method) %>% rename( Morgan = GMprob )
    P <- plyr::join_all( R, type = "inner" )

    ## The following pubchem_id's appear in training data and should be excluded from final evaluation
    vEx <- c( 3559, 47576, 55283, 1549008 )
    P <- P %>% filter( !(pubchem_id %in% vEx ) )

    ## Write output to file
    fn <- str_c( method, "-pred.tsv" )
    write_tsv( P, fn )
}

## Retrieves feature importance scores from a GBM model for a given feature space
gbm.featImp <- function( XY )
{
    ## Use a fixed random seed to allow for reproducibility
    set.seed( 100 )

    ## Downsample to training data
    X1 <- XY %>% filter( !is.na( Label ) ) %>%
        mutate( Label = factor( Label, levels=c(0,1) ) ) %>%
        select( -Drug, -pubchem_id )
    levels( X1$Label ) <- c( "Resistant", "Sensitive" )
    
    ## Set up a cross-validation schema
    fc <- trainControl( method="cv", savePred = "all", classProbs = TRUE )

    ## Perform cross-validation
    cv <- train( Label ~ ., data=X1, method="gbm", trControl=fc, verbose=FALSE )
    summary( cv$finalModel, plot = FALSE ) %>% filter( rel.inf > 0 ) %>%
        rename( Feature = var, Score = rel.inf )
}

## Main driver of functions in this file
main <- function()
{
    ## Compute predictions for all feature spaces
    pred.all()

    ## Compute feature importance scores for a GBM model in MACCS-binary feature space
    read_csv( "data/MACCSbinary.csv", col_types = cols() ) %>% gbm.featImp %>%
        write_csv( "MACCSbin-gbm-fimp.csv" )
}
