\name{pl99_data}
\alias{pl99_data}
\docType{data}
\title{
   A data.frame containing place name information used by the county and census tract
   data levels when creating the detailed run report in the satscanMapper package.
}
\description{
   Contains the place name information table to support the text 
   report on county and census tract place names.
}
\usage{data(pl99_data)}
\format{
  The pl99_data data.frame contains:
  \describe{
    \item{row.names}{- character - a search key containing the county 5 digit fips code
      and the place name within the county.
      }
    \item{stcoID}{- character - 5 digit state/county FIPS code containing the place name.
      }
    \item{plName}{- character - the place name.
      }
    \item{tracts_00}{- integer - the number of census tracts in the county
      in census year 2000.
      }
    \item{tracts_10}{- integer - the number of census tracts in the county
      in census year 2010.
     }
  }
  /br
  The place name index contains the information for all counties and census tracts.
  /br
  The data used to provide place names for the census tracts are derived 
  the Missouri Census Data Center correlation tables, the HUD.gov Census Tract to ZipCode
  crossover tables for 2000 and 2010, and the USPS.gov zip code database 
  for the 2000 and 2010 census years.  
  Additional place name information was provided by NCI and via location research in 
  the Internet.
}
\source{  
    CENSUS.GOV web site (boundary data), HUD.gov Census Tract to ZipCode crossover table,
    USPS Zip code database (1999), and R MAPS package (state and county names) for 
    the U. S. 2000 Census.
}
\keyword{datasets}
\keyword{Census2000}
\keyword{Census2010}
