\name{st99_M_data}
\alias{st99_M_data}
\docType{data}
\title{
   U. S. state centroid data for the states of Alaska, Hawaii, and Puerto Rico.
 }
\description{
   Contains the U. S. state information for Alaska, Hawaii, and Puerto Rico
   needed by the satscanMapper package to create the coordinates file to be 
   used by SaTScan (TM).  The main coordinates and boundary files  
   were modified to shift Alaska, Hawaii and Puerto Rico to better positions
   on the US map.  To be able to create a coordinates file based on
   the location IDs in the *.pop and *.cas input files to SaTScan (TM), 
   the coordinates must be based on the original Lat and Long and not the 
   moved locations.  This dataset contains the unmoved projected coordinates 
   of the of the centroids for the three moved states.
   Once the analysis is completed, only the mapping coordinates
   are needed.
 }
\usage{data(st99_M_data)}
\format{
  st99.ind00 structure:
  \describe{
    \item{row.names}{- character - 2 digit - State FIPS code.
      }
    \item{abbr}{- character - 2 character - State US PS Abbreviation.
      }
    \item{stName}{- character - State Name.
      }
    \item{rgID}{- integer - region ID number.
      }
    \item{rgName}{- character - region name. (NorthEast, South, MidWest, West, Offshore)
      }
    \item{dvID}{- integer - division ID number.
      }
    \item{dvName}{- character - division name 
      }
    \item{loc}{- character - postfix file string identifying the directory/package
      containing the census tract boundary datasets for this state. (Regs, East, West)
      }
    \item{county_00}{- integer - number of counties in state during census year 2000.
      }
    \item{county_10}{- integer - number of counties in state during census year 2010.
      }
    \item{tracts_00}{- integer - number of census tracts in state during census year 2000.
      }
    \item{tracts_10}{- integer - number of census tracts in state during census year 2010.
      }
    \item{change10}{- logical - indicates if the state had different counties and/or 
      boundaries in census year 2010.  If FALSE, use the 2000 boundary data.  If TRUE,
      use the supplemental 2010 boundary data.  Only 3 states required a 2010 set of 
      county boundary datasets (Alaska, Colorado, and Virginia)
      }
    \item{c_X and c_Y}{- integers - the X and Y cartesian coordinates of the center of 
      the state.
      }
  }
  This data.frame contains the centroid projected coordinates for Alaska, Hawaii,
  and Puerto Rico prior to their relocation on the map for better presentation.
  It has a duplicate structure and contains to the \var{st99_data} dataset to 
  allow the state data to be easily overlaid on the normal mapping data 
  when the coordinates file is being generated.
 }
\source{  
    CENSUS.GOV web site (boundary data), HUD.gov Census Tract to ZipCode crossover table,
    USPS Zip code database (1999), and R MAPS package (state and county names) for 
    the U. S. 2000 Census.
}
\keyword{datasets}
\keyword{Census2000}
\keyword{Census2010}

