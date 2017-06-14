\name{st99_data}
\alias{st99_data}
\docType{data}
\title{
   Data and index table for the state boundary information.
}
\description{
   Contains the required state information for satscanMapper to 
   be able to validate state location IDs in the SatScan (TM) data
   and produce the text reports.
 }
\usage{data(st99_data)}
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
  The state data table contains the information for all of the U. S. states, 
  DC and Puerto Rico.  
  This dataset is load for state, state/county, and state/county/census tract
  level data.
 }
\source{  
    CENSUS.GOV web site (boundary data), HUD.gov Census Tract to ZipCode crossover table,
    USPS Zip code database (1999), and R MAPS package (state and county names) for 
    the U. S. 2000 Census.
}
\keyword{datasets}
\keyword{Census2000}
\keyword{Census2010}

