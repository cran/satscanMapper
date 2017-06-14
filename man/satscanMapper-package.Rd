\name{satscanMapper-package}
\alias{satscanMapper-package}
\docType{package}
\title{A graphics package to create U. S. maps from SaTScan (TM) analysis 
    results data at the state, state/county or state/census tract levels.}
\description{
The satscanMapper package provides a means of creating U. S. maps at the state, 
county within a state, or census tracts within a state level (locations) using the 
results data from a SaTScan (TM) Spatial and Spatial-Temporal analysis.  
The package categorizes and color the observed over expected rations 
for locations used in the analysis.  The locations can be states, or 
counties or census tracts within a state.  
The package requires the location ID be the U. S. FIPS codes to be 
able to create a geo file and be able to map the locations (states, 
counties or census tracts within a state or DC.)  
This allows the package to link up 
the data to the US Census 2000 and 2010 boundary data.

The package uses the observed/expected ratios calculated by SaTScan (TM) 
to categorize and color the locations within each cluster.  
The package can use the ratios of each location or the ratio of the cluster being mapped.

The SaTScanMapper package uses the SeerMapper and Seer2010Mapper packages 
for the state, county and census tract boundaries and mapping functions.
SaTScanMapper then adds the cluster labels and outlines as requested.
}

\details{
    \tabular{ll}{
      Package:  \tab satscanMapper\cr
      Type:     \tab Package\cr
      Version:  \tab 1.0.0\cr
      Date:     \tab 2017-06-14\cr
      License:  \tab GPL-2\cr
      LazyLoad: \tab no\cr
    } 

The package contains a function to validate the population and case file location IDs 
against the contained boundary data and generate an appropriate geographic file for use
with SaTScan (TM). The minimum SaTScan results data files are:
\enumerate{
  \item Saved Parameter File
  \item Cluster Information File
  \item Location Information File
  \item Relative Risk Estimates File
 }
and must be in DBF format.

The boundary data for the states, counties and census tracts and mapping
functions from the \var{SeerMapper} package is used.

A state, county and place name (for census tracts) information is 
provided to help make the cluster reports more readable to the user.
Only the needed boundary data are loaded into memory.
}

\author{
Jim Pearson \email{jpearson@statnetconsulting.com}\cr
Maintainer: Jim Pearson \email{jpearson@statnetconsulting.com}\cr 
Package compiled by Jim Pearson \email{jpearson@statnetconsulting.com}
}

