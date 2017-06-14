\name{tr99_data}
\alias{tr99_data}
\docType{data}
\title{
   Census Tract data table is used by the satscanMapper package to 
   process census tract level location IDs for the census years 2000 and 2010.
}
\description{
   Contains the census tract information for census years 2000 and 2010 to assist in 
   the validation of the census tract level location IDs in user data
   and the generation of the text reports by the package on the SatScan (TM)
   cluster results.
 }
\usage{data(tr99_data)}
\format{
  The information contained in the tr99_data data.frame is:
  \describe{
    \item{row.names}{- character - 11 digit State/County/Census Tract FIPS code.
       }
    \item{stcoID}{- character - 5 digit state/county FIPS code containing the census tract.
       }
    \item{plKey}{- character - a placename key associated with the census tract.  The 
       key contains the state/county 5 digit fips code and the placename.
       }
    \item{c_X_00 and c_Y_00}{- integers - the X/Y cartesian coordinates of the center 
       of the census tract based on an equal area projection in meters for the census
       year of 2000.
       }
    \item{c_X_10 and c_Y_10}{- integers - the X/Y cartesian coordinates of the center 
       of the census tract based on an equal area projection in meters for the census
       year of 2010.
       }
    \item{y}{- integer - a bit string identifying the census years the census tract existed.
       If the value is 1 or 3, the tract existed in the 2000 census.  If the value
       is 2 or 3, the tract existed in the 2010 census.
       }
  }
  /br
  The tr99 data contain the information for all census tracts.  
  /br
  The place names in the tr99_data structure are derived from the HUD.gov Census Tract to ZipCode
  crossover tables and the USPS.gov zip code database for the 2000 and 2010 census.  
  Additional place name information was provided by NCI as a supplement.
}
\source{  
    CENSUS.GOV web site (boundary data), HUD.gov Census Tract to ZipCode crossover table,
    USPS Zip code database (1999), and R MAPS package (state and county names) for 
    the U. S. 2000 Census.
}
\keyword{datasets}
\keyword{Census2000}
\keyword{Census2010}
