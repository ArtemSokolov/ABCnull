## Composes .RData files for the data/ directory
##
## by Artem Sokolov

library( tidyverse )

## Load each raw dataset and store it as an .RData file to data/
MACCSbinary <- read_csv( "MACCSbinary.csv" ) %>% as.data.frame
save( MACCSbinary, file="../data/MACCSbinary.RData" )

MACCScount <- read_csv( "MACCScount.csv" ) %>% as.data.frame
save( MACCScount, file="../data/MACCScount.RData" )

Morgan <- read_csv( "Morgan.csv" ) %>% as.data.frame
save( Morgan, file="../data/Morgan.RData" )

PChem <- read_csv( "PChem.csv" ) %>% as.data.frame
save( PChem, file="../data/PChem.RData" )

ABCvaldata <- read_csv( "ABC-valdata.csv") %>% as.data.frame
save( ABCvaldata, file="../data/ABCvaldata.RData")
