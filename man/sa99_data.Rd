\name{sa99_data}
\alias{sa99_data}
\docType{data}
\title{
   Coordinates data table for the 20 U. S. NCI Seer Registry areas.
}
\description{
   Contains the required information and centroid coordinates for the 20 
   U. S. NCI Seer Registry areas.
}
\usage{data(sa99_data)}
\format{
  \preformatted{
  sa99_data structure:
     row.names - character - Seer Registry abbreviation.
     saName    - Character - Registry Name
     tracts    - integer   - Number of census tracts within state
     c_X_00    - integer   - centroid X in 2000 census 
     c_Y_00    - integer   - centroid Y in 2000 census
     c_X_10    - integer   - centroid X in 2010 census 
     c_Y_10    - integer   - centroid Y in 2010 census
  
  }
  /br
  The Seer Registry data set the information and centroid coordinates for
  the 20 NCI Seer Registries in the U. S.
  }
\source{  
    CENSUS.GOV web site (boundary data), HUD.gov Census Tract to ZipCode crossover table,
    USPS Zip code database (1999), and R MAPS package (state and county names) for 
    the U. S. 2000 Census.
}
\keyword{datasets}
