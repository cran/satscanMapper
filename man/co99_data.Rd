\name{co99_data}
\alias{co99_data}
\docType{data}
\title{
   U.S. County data table for all U. S. counties containing information to 
   support the county fips code validation and report generation 
   in the satscanMapper package.
}
\description{
   Contains all of the U. S. county information for all 50 states, the 
   District of Columbia and Puerto Rico needed by the satscanMapper package to 
   validate the user provided county level location IDs and support the 
   generation of the text reports on the SatScan (TM) results.
   }
\usage{data(co99_data)}
\format{
  The co99_data structure is:
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
  }
  /br
  The county data information contains the information for all U. S. counties.  }
\source{  
    CENSUS.GOV web site (boundary data), HUD.gov Census Tract to ZipCode crossover table,
    USPS Zip code database (1999), and R MAPS package (state and county names) for 
    the U. S. 2000 Census.
}
\keyword{datasets}
\keyword{Census2000}
\keyword{Census2010}
