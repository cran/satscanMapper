\name{tr99_M_data}
\alias{tr99_M_data}
\docType{data}
\title{
   U.S. census tract data table for Alaska, Hawaii, and Puerto Rico census tracts 
   containing the projected coordinates of the state's census tract centroids.
}
\description{
   Contains the U. S. census tract information for Alaska, Hawaii, and Puerto Rico
   needed by the satscanMapper package to create the coordinates file to be 
   used by SaTScan (TM).  The main coordinates and boundary files  
   were modified to shift Alaska, Hawaii and Puerto Rico to better positions
   on the US map.  To be able to create a coordinates file based on
   the location IDs in the *.pop and *.cas input files to SaTScan (TM), 
   the coordinates must be based on the original Lat and Long and not the 
   moved locations.  This dataset contains the unmoved projected coordinates 
   of the of the centroids for all of the census tracts in the three
   affected states.  Once the analysis is completed, only the mapping coordinates
   are needed.
 }
\usage{data(tr99_M_data)}
\format{
  The information contained in the tr99_data data.frame is:
  \describe{
    \item{row.names}{- character - 11 digit State/County/Census Tract FIPS code.
       }
    \item{stID}{- character - 2 digit state FIPS code containing the census tract.
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
  The tr99 M data contain the coordinates for all of the Alaska, Hawaii and Puerto Rico census tracts.  
  /br
  The place names in the tr99_M_data structure are derived from the HUD.gov Census Tract to ZipCode
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
