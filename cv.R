## Scores training data via cross-validation, for ROC curve generation
##
## by Artem Sokolov

library( tidyverse )
library( caret )

## Runs a single cross-validation split for a given set of models
## All models are evaluated against the same split
cv1 <- function( XY, vModels = c( "knn", "gbm", "glmnet", "svmLinear", "nnet" ) )
{
    ## Set up a cross-validation schema
    fc <- trainControl( method="cv", savePred = "all", classProbs = TRUE )

    ## Remove Drug and pubchem_id columns from consideration
    ## Convert Label column to {"Resistant","Sensitive"} factor
    X1 <- select( XY, -Drug, -pubchem_id ) %>%
        mutate( Label = factor( Label, levels=c(0,1) ) )
    levels(X1$Label) <- c( "Resistant", "Sensitive" )

    ## Set up the results data frame
    RR <- data_frame( rowIndex = 1:nrow(X1) )
    
    ## Traverse the models
    for( m in vModels )
    {
        cat( "Working with", m, "\n" )

        ## Each method has a slightly different flag for suppressing extraneous output
        if( m == "gbm" )
            cv <- train( Label ~ ., data=X1, method=m, trControl=fc, verbose = FALSE )
        else if( m == "nnet" )
            cv <- train( Label ~ ., data=X1, method=m, trControl=fc, trace = FALSE, MaxNWts=2000 )
        else
            cv <- train( Label ~ ., data=X1, method=m, trControl=fc )
        
        ## Average the predicted probability across paramenter values
        RR <- cv$pred %>% group_by( rowIndex ) %>%
            summarize( !!m := mean(Sensitive, na.rm=TRUE) ) %>%
            inner_join( RR, ., by="rowIndex" )
    }

    ## Match up predictions to the original pubchem_ids
    XY %>% mutate( rowIndex = 1:n() ) %>% select( pubchem_id, rowIndex ) %>%
        inner_join( RR, by="rowIndex" ) %>% select( -rowIndex )
}

## Runs n iterations of cross-validation
## fnIn - filename of the data
## n - number of cross-validation iterations
cvn <- function( fnIn, n = 100 )
{
    ## Load the data
    XY <- read_csv( fnIn, col_types=cols() ) %>% filter( !is.na( Label ) )

    ## Identify and remove features with no variance
    S <- XY %>% select( -Label, -Drug, -pubchem_id ) %>% summarize_all( sd ) %>%
        select_if( ~(. > 0) )
    XY <- XY %>% select( Label, Drug, pubchem_id, one_of( colnames(S) ) )
    
    ## Compose an output filename
    fnOut <- tools::file_path_sans_ext( fnIn ) %>% basename %>% str_c( "-cv.csv" )
    cat( "Will write output to", fnOut, "\n" )

    ## Run n iterations of cross-validation
    ## Update output file after every iteration, in case something breaks partway
    RR <- data_frame()
    for( i in 1:n )
    {
        cat( "Iteration", i, "\n" )
        RR1 <- cv1( XY ) %>% mutate( Iteration = i )
        RR <- bind_rows( RR, RR1 )
        write_csv( RR, fnOut )
    }
}

## Main driver of functionality in this file
main <- function()
{
    cvn( "data/MACCSbinary.csv", 100 )
    cvn( "data/MACCScount.csv", 100 )
    cvn( "data/PChem.csv", 100 )
    cvn( "data/Morgan.csv", 100 )
}
