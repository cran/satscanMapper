\name{CreateGeo4SS}
\alias{CreateGeo4SS}
\title{Function reads the population and case SaTScan (TM) input files and 
   creates a coordinates file based on the U. S. state, county, or census 
   tract FIPS code as the location IDs.  
 }
\description{
The satscanMapper creates maps based on the information and statistics 
data provide by in the SaTScan (TM) Results Files.
The package can map the 50 U.S. states, District of Columbia and the Puerto Rico 
territory at the state, county and census tract levels.

The satscanMapper requires the location IDs to be the U. S. FIPS codes 
for the States (and DC/Territories), State/Counties, or 
State/County/Census Tracts depending on the mapping level required for the data.  
The FIPS codes and boundary data used by this package is based on 
the 2000 and 2010 U. S. Census and were pulled from the Census Bureau web site.  
The \var{SeerMapper} package is used for the mapping and the source of 
all boundary data.

The \var{CreateGeo4SS} function creates a coordinates file based on the 
location IDs in the input population and case files to the SaTScan (TM) 
package.  The \var{satscanMapper} function takes the results and generates
maps.  To make sure the analysis and the clustering coordinates matches
the state, county and census tract boundary data used by \var{SeerMapper} to 
map the results, it is CRITICAL the initial input coordinates file matches 
the projection and coordinates used by \var{SeerMapper}.  

\var{CreateGeo4SS} reads the population and case files, collects the locations
IDs used, validates the Location IDs against the U. S. FIPS standards, and 
creates a coordinates (.geo) file for use as input to the SaTScan (TM) package.

}
\usage{CreateGeo4SS( path      = NULL,   # optional path to the case and population files.
                     pop       = NULL,   # filename of the population file
                     cas       = NULL,   # filename of the case file
                     
                     geo       = NULL,   # output filename for the coordinates file
                     overwrite = FALSE,  # overwrite permission on coordinates file
                                         #   def = FALSE
                     header    = TRUE,   # headers present on the pop and case files
                                         #   def = TRUE
                     censusYear= NULL    # def = 2000.
                 )
}

\arguments{
  \item{path}{is the name of the directory containing the input population and case
    files.  If present, it will also be used for the coordinates output file.
    The path= value will be prepended to the pop, case and geo filenames if needed 
    to access the file.
   }
  \item{pop}{is the file name of the input population file.  The value
    can be a full filename with a path or just a filename to be used with
    the path= call parameter.  
    The format of the population file must be the Location ID followed 
    by a white-space separator and any number of fields required. Only the 
    location ID is evaluation.  A valid pop or case file must be 
    specified.
   }  
  \item{cas}{is the file name of the input case file.  The value can be a 
    full filename with a path or just a filename to be used with the path=
    call parameter.
    The format of the case file must be the Location ID followed 
    by a white-space separator and any number of fields required. Only the 
    location ID is evaluation.  A valid pop or case file must be specified.
   }  
  \item{geo}{is the file name of the output coordinates file.  If no extension is
    provided, ".geo" will be used.  It the geo= call parameter is not 
    specified, the base filename and path of the case and population 
    will be used to create a coordinates file.  The value can be a full filename
    with path or a filename to be used with the path= call parameter.
    The file specified is check to see if it already exists.  If the overwrite=
    call parameter is set to TRUE, the existing file will be overwritten.
    If overwriting is not permitted and the file exists, no coordinates
    file is created and the caller informed .
    Note: A data.frame representing the coordinates file, is returned to the 
    \var{CreateGeo4SS} caller and the functions value.
   }  
  \item{overwrite}{is a logical value.  The \var{overwrite} call parameter
    determines what will happen if the output coordinates filename already
    exists.  If TRUE, the fill will be overwritten.  If FALSE, the file will 
    not be overwritten and no output file will be generated.  The default value
    is FALSE.
   }
  \item{header}{is a logical value.  The \var{header} call parameter is 
    passed through to the read.table functions to indicate if the input files have a
    header line as the first line of the file.  If set to TRUE, 
    the first line is processed as a field name line instead of data.
    If set to FALSE, the first line of the file is treated as data.
    The default value is TRUE.  The \var{header} call parameter is also 
    used to control whether or not to include a header line in the coordinate output
    file.  Both the population and case files must have headers or not have headers.
   }
  \item{censusYear}{This is a character value of the census year of the location IDs
    in the population and case data records or the SatScan (TM) results data files.
    The acceptable values for the censusYear parameter are "2000" and "2010".
    The default value is "2000" for the 2000 census year.
   }
 }
\details{
  The \var{CreateGeo4SS} function creates a coordinates file for input to 
  the SaTScan (TM) package.  By using the \var{CreateGeo4SS} function,
  the caller an ensure the location IDs used in the population and 
  case files match the FIPS codes for the specified census years and that 
  the coordinates used in the analysis are equal-area projection coordinates 
  that match the boundary data stored in the \var{SeerMapper} package (which draws the 
  maps.)  Once the coordinate file is created, the population, case and coordinates
  file can be used by the SaTScan (TM) package and the results mapped by the 
  \var{satscanMapper} function in this package.

  The locations IDs used in the population and case files must be the 
  U. S. FIPS states (2-digit), state/counties (5-digit), 
  or state/county/census tract (11-digit) codes based on the geographical units
  being mapped.  Since the FIPS codes are different in the 2000 and 2010 
  census years, the caller needs to notify \var{CreateGeo4SS} and \var{satscanMapper}
  via the \var{censusYear} call parameter if the location identifiers are from 
  the 2010 census years.  The packages default to using the 2000 census FIPS codes 
  and boundary data.  If the location ID based on the 2010 census, set \var{censusYear}
  call parameter to "2010" and the packages adjust.

  This package will only map the boundaries of states that are linked to data or have 
  counties or tracts linked to data within their boundaries.
  
}
\value{None}
\author{ 
Jim Pearson and Linda Pickle of StatNet Consulting, LLC, Gaithersburg, MD 
}
\examples{

######
##
##  Example # 1 - creating a coordinates file for a SaTScan (TM) analysis
##
##  This example access the population and case data files used in a 
##  SaTScan (TM) analysis:
##     a) read the files to build a list of all of the location IDs 
##          involved in the analysis, 
##     b) validates the location IDs against the U. S. Census web site's 
##          list of Location IDs for that year and geographic level 
##          (State, State/County or State/County/Census Tract), and 
##     c) returns a list of centroid X, Y coordinates for each location ID.
##
##  The coordinates file can be written to disk or returned as a data.frame 
##  to the caller as the return value of the CreateGeo4SS function.
##
##  The coordinates values returned are based on an equal-area projection 
##  of the lat and long values of the area centroids based on the 
##  area boundary data used to map the areas.
##
##  The example population and case data files distributed with the package 
##  are located in the installation directory for the satscanMapper package.
##
#
##  return a DF of the coordination data
#
\dontrun{

SSInstDir <- system.file(package="satscanMapper", "extdata")

PopFile   <- "inc_noadj_cir.pop"

CasFile   <- "inc_noadj_cir.cas"

GeoDF     <- CreateGeo4SS(path=SSInstDir, 
                          pop=PopFile, 

                          cas=CasFile, 
                          censusYear="2010",
                          header=FALSE)

str(GeoDF)

}
#  Example removed due to time constraints
#
######


#####
#
#  This example #2 creates a coordinates file.
#

\dontrun{

SSInstDir <- system.file(package="satscanMapper", "extdata")

PopFile   <- "inc_noadj_cir.pop"

CasFile   <- "inc_noadj_cir.cas"

GeoDF     <- CreateGeo4SS(path=SSInstDir, 
                          pop=PopFile, 
                          cas=CasFile, 
                          geo="inc_noadj_cir.geo",
                          censusYear="2010",
                          header=FALSE)

# output file written to the extdata directory of the package installation.

}

####

print(" End of CreateGeo4SS examples.")

####



}
