## Composes .RData files for the data/ directory
##
## by Artem Sokolov

library( tidyverse )

## Load each raw dataset and store it as an .RData file to data/
MACCSbinary <- read_csv( "MACCSbinary.csv" ) %>% as.data.frame
save( MACCSbinary, file="../data/MACCSbinary.RData" )

##MACCScount <- read_csv( "MACCScount.csv" )
##Morgan <- read_csv( "Morgan.csv" )
##PChem <- read_csv( "PChem.csv" )

