\name{satscanMapper}
\alias{satscanMapper}
\title{Function to map Results data produced by a SaTScan (TM) cluster analysis}
\description{
 The satscanMapper creates maps based on the information and statistics 
 data provide by in the SaTScan (TM) Results Files.
 The package can map the 50 U.S. states, District of Columbia and the Puerto Rico Territory at 
 the state, county and census tract levels.
 
 The locations are identified using the U. S. FIPS codes for the State (and DC/Territories), 
 Counties with in a state, and Census Tracts with a state.  
 The boundary data used by this package is based on the 2000 and 2010 U. S. 
 Census and were pulled from the Census Bureau web site.  The \var{SeerMapper}
 package is used for the mapping and the source of all boundary data.
}
\usage{satscanMapper( resDir      = NULL,           # no default, must be provided
                      prmFile     = NULL,           # no default, must be provided
                      outDir      = NULL,           # default = NULL, uses value in resDir
                      
                      censusYear  = NULL,           # default = census year 2000.
                      categ       = NULL,           # default = 5
                      title       = NULL,           # no default, optional
                      pValue      = NULL,           # default = 0.05
                      report      = NULL,           # default = TRUE
                      runId       = NULL,           # default = ""
                      bndyCol     = NULL,           # default = grey(0.9)
                      label       = NULL,           # default = TRUE
                      outline     = NULL,           # default = TRUE
                      locO_EMap   = NULL,           # default = TRUE
                      clusO_EMap  = NULL            # default = TRUE
                 )
       }

\arguments{
  \item{resDir}{is the name of the directory containing the SaTScan (TM) saved session parameter 
    file (.PRM).  This directory will be used to hold the output graphics (PDF) and 
    text output report (.txt) files. 
   }
  \item{prmFile}{is the file name of the SaTScan (TM) saved session parameter file (.PRM).
    The saved session parameter file contains all of the information about 
    the SaTScan analysis to be mapped, including the location of the 
    results cluster, location and RR files.
    The base of the .PRM filename (without extension) is used to create the 
    filenames for the graphic and text output files.  If the output filename already exists, 
    a new name based on the original name is created and used. or you can specify
    the \option{runId} to prepend a few characters to the filenames to make them 
    unique for each \var{satscanMapper} run.
   }  
  \item{outDir}{is the name of the directory to write the PDF graphic map file and 
    the TXT summary report file.  Normally, the output files are written to the 
    directory specified in the \option{resDir} parameter to keep the input and output information
    in the same location.  However, if the output files need to be written to a different
    directory, the \option{outDir} parameter can be used to specify the different
    directory. If the \option{outDir} parameter is missing or \var{NULL}, the value of 
    the \option{resDir} parameter is used as the output directory.  
    If parameter is set to \var{""}, the current working directory will be used when the 
    output files are written.  If the parameter is set to a
    character string, it must be a valid/existing directory and will be used as the destination when 
    the output files are written.    The default value is \var{NULL}.
   }
  \item{censusYear}{This is a character value of the census year of the location IDs
    in the population and case data records or the SatScan (TM) results data files.
    The acceptable values for the censusYear parameter are "2000" and "2010".
    The default value is "2000" for the 2000 census year.
   }
  \item{categ}{This is an integer value specifying the number of categories to be 
    used when coloring the locations in the cluster maps.  A preset of Observed/Expected 
    value break points is used based on the number of categories specified 
    in the \option{categ} argument.  \option{categ} can range from \var{3} to var{10}.
    The default \option{categ} value is \var{5}.
    The breakpoints for each categ value from \var{3} to \var{10} are 
    hard coded into the package as follows:
    \preformatted{
       categ   breakpoint list
         3     -Inf,                       0.9, 1.1,                Inf
         4     -Inf,                       0.9, 1.1,           3.0, Inf
         5     -Inf,            0.5,       0.9, 1.1,           3.0, Inf
         6     -Inf,            0.5,       0.9, 1.1,      2.0, 3.0, Inf
         7     -Inf, 0.33,      0.5,       0.9, 1.1,      2.0, 3.0, Inf
         8     -Inf, 0.33,      0.5,       0.9, 1.1, 1.5, 2.0, 3.0, Inf
         9     -Inf, 0.33,      0.5, 0.67, 0.9, 1.1, 1.5, 2.0, 3.0, Inf
        10     -Inf, 0.33, 0.4, 0.5, 0.67, 0.9, 1.1, 1.5, 2.0, 3.0, Inf
      }
    Future plans include allowing the caller to specify their own set of break points.
   }
  \item{title}{A character string(vector) used as the main title for the cluster
    maps.  A \var{title} may be one or two lines. If two lines, code 
    \code{title=c("line 1", "line 2")}. 
    for the \var{title} call parameter.  If no title is provided (\var{NULL} or \var{""}), 
    the prmFile argument is used for 
    the main title.  The default is \var{NULL}, the package uses the prmFile argument
    as the title.
   }
  \item{pValue}{The package uses a default P Value of \var{0.05}. Only cluster with a 
    P Value < \var{0.05} are colored on the resulting map by default.  
    If the P Value needs to be a different value, then the user
    can specify the override P Value using the PValue argument.  
    The \code{PValue} argument must be between var{0.001} and var{0.5} to be valid.
   }
  \item{report}{This is logical value.  If \var{TRUE} a text report of the statistics from each 
    cluster and location in a tier format by State, County, Place, and Location (Census Tract) 
    will be created.  If var{FALSE}, no report is generated.  The default is \var{TRUE}.
   }
  \item{runId}{If a user is creating multiple mappings of the same SaTScan (TM) 
    Results Data Files, but with different arguments and options, they will want 
    to save the PDF and report text files for each run.  Normally the PDF and report text
    files would be overlaid. The \option{runId} argument provides a way to modify the 
    output PDF and text filenames to make a unique filename for the output.  
    The \option{runId} string is inserted near the end of the output filename.  This allows
    multiple mapping runs against the same data and a collection of resulting maps.
    If \option{runId} is not specified, the package will append to the '.prm' filename the 
    string of "-Rnn", where 'nn' is a number from 01 to 99.  If an output set already exists, 
    'nn' is increased by 1 until an unused output file name is discovered.  
   }
  \item{label}{This is a logical argument that indicates whether or not the outlines of 
    each cluster on a map is labeled with the cluster number from the SaTScan (TM) data.
    If a cluster spans multiple years, the number of years is added to the label. 
    For example: "1-2", this is cluster number 1 and the cluster existing over a two
    year period including the current year being mapped.  
    The default value for \option{label} is var{TRUE}.
   }
  \item{outline}{This is a logical value specifying whether or not draw the outline of the 
    cluster on the cluster map. If \code{TRUE}, an outline is draw.  If \var{FALSE}, no outlines
    are drawn.
    The default value is \var{TRUE}.
   }
  \item{bndyCol}{is a character string that defines the color to use for the 
    area boundaries.  The default color value is \var{grey50}.
    The value provided is verified against the values returned from 
    the \var{colors()} function.
   }
  \item{locO_EMap}{a logical variable. This parameter enables the use of the location
    Obs/Exp ratio to categorize and color the data areas.  If \var{TRUE}. the location 
    ratios are used.  If \var{FALSE}, the location Obs/Exp ratio is not used.
    The default value is \var{TRUE}.  
   }
  \item{clusO_EMap}{a logical variable. This parameter enables the use of the cluster
    Obs/Exp ratio to categorize and color the data areas.  If \var{TRUE}. the cluster 
    ratios are used.  If \var{FALSE}, the cluster Obs/Exp ratio is not used.
    The default value is \var{FALSE}. 
   }
 }
\details{
The satscanMapper function creates one or more maps for the SaTScan (TM) cluster analysis 
results data.   The following must be done to facilitate the use of SaTScan (TM) cluster 
analysis and the satscanMapper package:
\enumerate{
  \item{The coordinates file must be generated by the \var{CreateGeo4SS} function
    in the package to ensure the area centroids used in the analysis match the 
    boundary projection used in the area mapping.
   }
  \item{The Coordinates type must be set to "Cartesian".  The centroid values 
    returned by the \var{CreateGeo4SS} function are based on an equal-area projection
    and not Lat/Long coordinates.
   }
  \item{Time Precision must be set to "Years".  The mapping package only supports 
    time units of years.
   }
  \item{Since mapping is spatial, only the Retrospective Purely Spatial and Space-Time 
    and the Prospective Space-Time analysis are supported.
   }
  \item{There restrictions on the Probability Model used are imposed by the type of 
    analysis.
   }
  \item{Time Aggregation is support for units of "years" for any number of years.
    The default is 1 year.  This parameter is generally not used.
   }
  \item{Both the circular and elliptic Spatial Window Shapes are supported.
   }
  \item{The Minimum Temporal Cluster Size is supported for any number of years.  When 
    specified, maps are generated on a sliding time line for group of years.
   }
  \item{The main results file and all Column Output files use the same "Main Results File"
    name. This information is stored in the saved session *.PRM file.  Any value
    may be specified.
   }
  \item{The Cluster Information dbase box must be checked.
   }
  \item{The Location Information dbase box must be checked.
   }
  \item{The Risk Estimate for Each Location dbase box must be checked.
   }
  \item{IMPORTANT: Just before executing the SaTScan (TM) analysis, save the 
    session *.prm file in an appropriate directory. Menu -> File -> Save Session or Save Session as.
    If you make changes to the analysis parameters make sure to re-save the 
    session file.
   }
 }
    
The saved session *.PRM file provides the main mapping function with 
all of the information to generate the maps of the analysis results.
The saved session parameters are checked to ensure the mapping can
be done and any errors reported to the user.

The \var{CreateGeo4SS} function should be used to evaluate the location IDs
in the population and/or case data files and generate the coordinates(geo) file
for use in the SaTScan (TM) analysis.  

The Location IDs used in the original population, case and coordinates(geo) input files 
and the output results files must be in one of three formats based on 
the mapping level desired: the 2 digit U. S. fips code for state data, 
the 5 digit U. S. fips code for county data, 
or 11 digit U. S. fips code for census tract data.  
If the Location IDs for the census year being mapped, do not match the 
boundary data, warnings will be generated to inform the user of the issue.
Check to make sure the correct census year was specified.
The boundary data and FIPS codes were pulled from the CENSUS.GOV web site
for the 2000 and 2010 census years and should be complete.

If 2 digit fips code IDs are used, the SaTScan (TM) results information 
is mapped at the state, district and territory level. When states are mapped,
the boundaries for all of the states in the US will be drawn.

If 5 digit fips code IDs are used, the SaTScan (TM) results information 
is mapped at the county level for the states containing the counties with data.  
State boundaries are drawn around the counties as a reference. Only the county 
boundaries in states with data will be drawn.  Therefore, is the analysis results
contain counties in multiple states, all of the boundaries for the counties in these
states are drawn.  

If 11 digit fips code IDs are used, the SaTScan (TM) location information is 
mapped at the census tract level for the states containing the census tracts with data.
State and county boundaries are drawn for any state contain census tracts referenced
by an analysis.

The area categorization and coloring are based on either the location and cluster ODE ratios
reported in the results files.  The ratio used is controlled by the \var{locO_EMap} and 
\var{clusO_EMap} parameters.

The \option{bndyCol} parameter specifies the color used for the data area boundaries 
(state, county or census tract).  It's value must be a colors that is either from 
the \var{colors()} name list or is a
\var{#hhhhhh} or \var{#hhhhhhhh} format as a character vector.

Only areas with data are mapped.  To drawn all of the data level boundaries up 
to the next level, specify \option{fill} = \var{TRUE}.   
For example: All tract boundaries will be drawn up to 
their county level for tract data. For county data, all counties will be drawn up 
to the state boundary.

At the current time, the boundary data is based on the boundary data in the SeerMapper
and Seer2010Mapper packages.  An additional parameter is required to specify which
census year is to be used:  2000 or 2010.

If a text report is required (\option{report} = \var{TRUE}), the package provides 
a listing of the Cluster Information, Location Information and Risk data 
in a tabular format for the user's review.  The location information is listed 
for each cluster in a tier format with the following levels:
\preformatted{
  State
     County
        PlaceName
           Census Tract
}
If the locations are counties, county names are provided, but not the 
Place names or the Census Tract IDs levels.

The critical SaTScan (TM) run options the must be set are:
\preformatted{
   StartDate = yyyy/mm/dd   
   EndDate = yyyy/mm/dd

   AnalysisType = 4   (default)
       Valid Types: 1 = Purely Spatial, 
                    3 = Retrospective Space-time, 
                    4 = Prospective Space-Time. 
       Not-Supported Types:
                    2 = Purely Temporal, 
                    5 = Spatial Variation in Temporal Trends,
                    6 = Prospective Purely Temporal 
       

   IncludeRelativeRisksCensusAreasDBase       = y     # *.rr.dbf
   CensusAreasReportedClustersDBase           = y     # *.col.dbf
   MostLikelyClusterEachCentroidDBase         = y     # *.gis.dbf
   MostLikelyClusterCaseInfoEachCentroidDBase = y     # *.sci.dbf
}
The \var{satscanMapper} package uses the \var{SeerMapper} package to generate 
all of the maps.  Other R packages used by \var{satscanMapper} are: base, utils, 
graphics, grDevices, stats, RColorBrewer, stringr, sp, and foreign.
}
\value{None}
\author{ 
Jim Pearson and Linda Pickle of StatNet Consulting, LLC, Gaithersburg, MD 
}
\examples{
######
#
#  These examples focus on mapping the SaTScan (TM) data after the analysis
#  has been completed.
#  See the section on the CreateGeo4SS function for an example of how 
#  to build a coordinates file (.geo) that matches the satscanMapper 
#  boundary data.
#
######

######
#
#  Example # 1 - Mapping existing SaTScan (TM) results - "USStateLung" result files
#     Make sure SaTScan (TM) results files are still located in the directories
#     documented in the session saved (.prm) file.
#
#     This example maps data for the US States containing circular clusters.
#
#   Get location of .prm file and location to write mapping output files
#
#\dontrun{
#tempdirOut <- tempdir()
#cat("tempDirOut:",tempdirOut,"\n")

SSMInstDir <- system.file(package="satscanMapper", "extdata")
cat("SSMInstDir:",SSMInstDir,"\n")

TT =c("Contiguous US States All Lung Cancer Mortality, 2004")

satscanMapper(resDir   = SSMInstDir,                # path to .prm file location for output files. 
              prmFile  = "stateLung.prm", 
              outDir   = "",
              title    = TT    
             )
#}
#
##  end of example # 1  Removed due to time constraints.
##
######

######
#
#  Example # 2 - Mapping existing SaTScan (TM) results - "USCountyLung" result files
#     named "inc_noadj_cir" files
#     Make sure SaTScan (TM) results files are still located in the directories
#     documented in the session saved (.prm) file.
#
\dontrun{
#tempdirOut <- tempdir()
#cat("tempDirOut:",tempdirOut,"\n")

SSMInstDir <- system.file(package="satscanMapper", "extdata")

TT = "Contiguous US County Female Breast Cancer Incidence, 2009-2013"

satscanMapper(resDir=SSMInstDir,              # path to .prm file location for output files. 
              prmFile    = "inc_noadj_ellip_hilo_10_nosp.prm", 
              outDir     = "",
              categ      = 7,
              title      = TT,
              censusYear = "2010"
             )

}
#
#  end of example # 2
#
#####

print("end of examples")
}
