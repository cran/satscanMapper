\name{co99_M_data}
\alias{co99_M_data}
\docType{data}
\title{
   U. S. County data table for Alaska, Hawaii, and Puerto Rico counties containing
   the projected coordinates of the state's county centroids.
}
\description{
   Contains the U. S. county information for Alaska, Hawaii, and Puerto Rico
   needed by the satscanMapper package to create the coordinates file to be 
   used by SaTScan (TM).  The main coordinates and boundary files  
   were modified to shift Alaska, Hawaii and Puerto Rico to better positions
   on the US map.  To be able to create a coordinates file based on
   the location IDs in the *.pop and *.cas input files to SaTScan (TM), 
   the coordinates must be based on the original Lat and Long and not the 
   moved locations.  This dataset contains the unmoved projected coordinates 
   of the of the centroids for all of the counties in the three
   affected states.  Once the analysis is completed, only the mapping coordinates
   are needed.
   }
\usage{data(co99_M_data)}
\format{
  The co99_M_data structure is:
  \describe{
     \item{row.names}{- character - 5 digit State/County FIPS code.
       }
     \item{stID}{- character - 2 digit - State FIPS code.
       }
     \item{stName}{- character - state full name.
       }
     \item{coName}{- character - county full name.
       }
     \item{saID}{- character - Seer Registry Area identifier (abbreviation).
       }
     \item{c_X_00}{- integer - the X value of the cartesian coordinate for the 
       center of the county in an equal area projection for the census year 2000.
       }
     \item{c_Y_00}{- integer - the Y value of the cartesian coordinate for the 
       center of the county in an equal area projection for the census year 2000.
       }
     \item{c_X_10}{- integer - the X value of the cartesian coordinate for the 
       center of the county in an equal area projection for the census year 2010.
       }
     \item{c_Y_10}{- integer - the Y value of the cartesian coordinate for the 
       center of the county in an equal area projection for the census year 2010.
       }
     \item{tracts_00}{- integer - the number of census tracts in the county
       in census year 2000.
       }
     \item{tracts_10}{- integer - the number of census tracts in the county
       in census year 2010.
       }
    \item{y}{- integer - a bit string identifying the census years the census tract existed.
       If the value is 1 or 3, the tract existed in the 2000 census.  If the value
       is 2 or 3, the tract existed in the 2010 census.
       }
  }
  /br
  The county data information contains the information for the states of 
  Alaska, Hawaii, and Puerto Rico.}
\source{  
    CENSUS.GOV web site (boundary data), HUD.gov Census Tract to ZipCode crossover table,
    USPS Zip code database (1999), and R MAPS package (state and county names) for 
    the U. S. 2000 Census.
}
\keyword{datasets}
\keyword{Census2000}
\keyword{Census2010}
