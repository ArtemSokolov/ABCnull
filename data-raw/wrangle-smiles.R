## Composes ABCSMILES.RData files for the data/ directory
##
## by Artem Sokolov

library( tidyverse )

## Retrieve the list of PubchemIDs to query with
MACCSbinary <- read_csv( "MACCSbinary.csv" ) %>% select( Drug, pubchem_id )

## Compose a PubChem query
ids <- str_flatten( MACCSbinary$pubchem_id, "," )
qry <- str_c( "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", ids, "/property/IsomericSmiles,CanonicalSMILES/CSV" )
SM <- read_csv( qry )

## Join everything into a common table
ABCSMILES <- rename( SM, pubchem_id = CID ) %>% inner_join( MACCSbinary, . ) %>% as.data.frame

## Write out the data to files
ABCSMILES %>% write_csv( "ABC-SMILES.csv" )
save( ABCSMILES, file="../data/ABCSMILES.RData" )
