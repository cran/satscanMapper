#     
#     
#  R-Code: Generate cluster circles or elipses on US map based on 
#     clusters identifed by SatScan.
#     
#     Program uses SatScan's output files to identify the clusters and 
#     size of the cluster.  This is then ploted on an US map.
#
#  Updating package from NCI SaTScan R Code mapper to general package.
#     Based on release  0.94 to NCI 
#     Prepared by Jim Pearson - December 12, 2013
#        Re-designed for release by Jim Pearson - January and Feb, 2017 
#        to use the SeerMapper boundary data and mapping logic.
#        Uses new data.frame support structure that has been expanded to 
#        support census years 2000 and 2010.
#        The original NCI.Ind00 dataset containing several data.frames has 
#        been replaced with a dataset for each data.frame:  
#          st99_data -> state data and information
#          co99_data -> state/county data and information
#          tr99_data -> state/county/tract data and information
#          pl99_data -> place name data to support text reports involving census tracts.
#
#  Load supporting libraries and routines. Make sure all library/packages have 
#    been installed.   The only key private package is SeerMapper and it's 
#    census tract boundary files.
#
#
debugFlag <- FALSE
#
#library(graphics)      - yes
#library(foreign)       - yes - reading dbf and csv files.
#library(stringr)       - yes
#library(RColorBrewer)  - yes
#library(stats)         - yes - categorization
#library(SeerMapper)    - yes

##library(sp)           - no  - deferred to SeerMapper
##library(data.table)   - no
##library(maptools)     - no
##library(gplots)       - no
##library(maps)         - no
##library(mapproj)      - no
##library(rgdal)        - no
##library(matrixStats)  - no 
#
# Change Log
#
#   17/05/09 - completed conversion to use boundary projections and centroids
#            - calculated adjustment to user cartesian coordinates to convert
#              radius and ellipical values to correct magnitude.
#            - looked at how to handle Lat/Long coordinates when used.
#            - looked at how to handle overlapping clusters.  Reorder 
#              color fill to keep lead cluster's color first on clusO_E.
#   17/05/22 - Need to be able to match up categories and colors legend between SeerMapper 
#              and satscanMapper's selection.  Initial solution, change satscanMapper
#              to select colors based on number of categories not 11, pell off top and bottom, and 
#              select from remainder for the colors.  The cluster outline colors are still picked from
#              the extremes of the 11 color group. 
#            - Corrected calculation of label Y value.
#            - Increased size of labels from 0.3 to 0.5 cex.
#   17/06/10 - Fix problems with census tract reports place index
#            - change off ODE to obs/exp in outputs.
#            - change call parameters locODEMap and clusODEMap to locO_EMap and clusO_EMap.
#   18/06/29 - add outDir call parameter to allow output graphic and 
#              text report files to be written a directory other than
#              the distribution/test directory containing the SatScan(TM) 
#              results files.  CRAN is their release testing no longer 
#              allows a package to write to it's own release/test data
#              directory, like the user do.  An alternate directory 
#              must be provided to allow the package to pass
#              release testing.  outDir allows no output files to be 
#              written, write to a dedicated directory, write to 
#              the default directory.
# 
iRes <- require(SeerMapper)
if (!iRes) {
    stop("SeerMapper package could be found and loaded.")

}
#
#
#######
#
# Instructions:
#
#  To run the program correctly --> use the function call.
#
#  Provide the Dir and ResultsFile parameters to permit the package to find the SaTScan result files.
#
#  The package scans the *.prm file and verifies the correct results files should be present
#  and gathers most of the run parameters:
#      SSForm, SSShape, Run Period (start and end dates)
#
#  Different from previous versions, the *.pop, *.cas, and *.geo files are not used.
#
#  A new function has been included to read the location ids from the *.pop and *.cas
#  files and generate the *.geo file with cartisan coordinates for the centroids that
#  match the boundary data contained in the mapping portion of the package.
#
#  The *.rr.dbf file is used to get a full list of locations involved in the analysis, if available.
#  otherwise the GEO file is read.
#
#  The Rate data files and mapping are options and may not be provided.  The rate data file must
#   be located in the same directory as the result files.
#
#  The package relies on the SeerMapper package to contain the boundry data 
#   and mapping functions for use by this package.  Several tables are keep 
#   in this table complement the SeerMapper data and must be sync'd between 
#   the package whenever boundary data is changed.
#  satscanMapper uses the SM_GlobInit, SM_Build,  and SM_Mapper functions from 
#   SeerMapper. The call parameter validation, hatching and categorizing 
#   functions are not use.
#  Once the boundary data is loaded (once per call to satscanMapper), it is 
#   reused for each required map.  satscanMapper will make changes to the 
#   SeerMapper structures (rPM and MV) to effect the changes and then call 
#   SM_Mapper to do the mapping on a single page.  SM_Mapper returns the 
#   x-limit and y-limit plot values used to permit satscanMapper to overlay 
#   additional labels, circles, ellipses, etc. accurately on the map.
#
########################################################
#
#   Functions --------------
#
#  draw the ellipse on the plot.
#

ellipsePoints <- function(a,b, alpha = 0, loc = c(0,0), n = 201)
{
    ##
    ## Purpose: ellipse points,radially equispaced, given geometric par.s
    ## -------------------------------------------------------------------------
    ## Arguments: a, b : length of half axes in (x,y) direction
    ##            alpha: angle (in degrees) for rotation
    ##            loc  : center of ellipse
    ##            n    : number of points
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 19 Mar 2002, 16:26
    ##
    
    B <- min(a,b)
    A <- max(a,b)
    ## B <= A  - given
    d2 <- (A-B)*(A+B)                   #= A^2 - B^2
    phi <- 2*pi*seq(0,1, len = n)
    sp <- sin(phi)
    cp <- cos(phi)
    r <- a*b / sqrt(B^2 + d2 * sp^2)
    xy <- r * cbind(cp, sp)
    ## xy are the ellipse points for alpha = 0 and loc = (0,0)
    al <- alpha * pi/180   # convert to radiants
    ca <- cos(al)
    sa <- sin(al)
    xy %*% rbind(c(ca, sa), c(-sa, ca)) + cbind(rep(loc[1],n),
                                                rep(loc[2],n))
}

####
#
#  routines to handle the problem with CUT and non-unique break point lists.
#
####

RateQuan <- function(brkpt, data) 
  {
    #  generate the quantile list for the break points for Rate data.
    
    wQ <- quantile(data,probs=brkpt,na.rm=T)
    wQ <- c(wQ,Inf)    # may not be needed if last break point is the maximum.
    wQ[1] <- -Inf      # change "0%" value to -Inf
    return(wQ)
  }
  
RateCutAdj <- function(brkpt) 
   {
      #
      # review break point list and adjust when duplicate values are found.
      # also create a category label list for use later with "No Values" if needed.
      #
      wC = brkpt

      #
      #  build literal category list for legend
      #    and make adjustments to the RateCut break points for equal pointsl.
      #

      wCat <- rep("",length(brkpt))
      for (i in 1:(length(brkpt)-1))
         {
           if (brkpt[i] >= brkpt[i+1]) 
             {
              vv = "No Values"
              wC[i+1] = wC[i]+0.000001
             } else {
              vv = gsub(" ","",paste("(",brkpt[i],",",brkpt[i+1],"]",seq=""))
             }
           wCat[i] = vv
         }
      #  return data.frame with adjusted break point list and category labels.
      wDF = data.frame(cuts = wC, cats = wCat)
   }

#
#  Common code as function to plot cluster outlines 
#

  PlotClusterOutline <- function ( ESetLen, ESet, CLab=TRUE, CLabCol="black")
        {
          if (ESetLen > 0 )
            {
              # overlay rate map with cluster outlines.
              for (ind in c(1:ESetLen))
                {
                  if (sum(ESet$E_MAJOR[ind], ESet$E_MINOR[ind]) > 0)
                    {
                      aAngle <- ESet$E_ANGLE[ind]
                      aLoc <- c(ESet$t_X[ind],ESet$t_Y[ind])   # center of Ellipse or Circle
          
                      # Generate set of points for Ellipse or Circle
                      epset <- ellipsePoints(
                             ESet$E_MAJOR[ind],ESet$E_MINOR[ind]
                            ,alpha=aAngle
                            ,loc=aLoc)
                           
                      #  Plot ellipse or circles for Cluster with P_value <= 0.05
                      #cat("plotting E/C: Ang:",aAngle,"  Loc:",aLoc,"  col:",ESet$OutCol[ind],"\n")
                      #cat("points:\n")
                      #print(epset)
                      
                      polygon(epset
                             ,border=ESet$OutCol[ind]
                             ,lwd=0.75   # changed from 0.5   17/06/12
                           )
                    }
                }
              if (CLab)
                {
                  par(new=T) 
                  xuser <- par("usr")   # get x and y widths.
          
                  labX  <- ESet$Lab_X
                  rX    <- xuser[2]
                
                  OffR       <- (labX >= rX) 
                  labX[OffR] <- rX - strwidth(ESet$Lab[OffR],cex=0.6)
                                   
                  #cat("Labels:\n")
                  #print(ESet[,c("Lab","Lab_X","t_X","t_Y")])
                  text(labX,ESet$t_Y,ESet$Lab,cex=0.6,col=CLabCol)       
                }
            }
        }
 
 #
 #  Get Year from date
 #
 getYear <- function(x) {
   as.POSIXlt(x)$year + 1900L
 }
 
 #
 #  print callVarList
 #
 printCallVarList <- function(x) {
    if (!is.list(x)) {
       xmsg <- paste0("callVarList is not a list structure. Printed in raw format.")
       warning(xmsg, call.=FALSE)
       print(x)
    } else {
       xN <- sort(names(x))
       #cat("xN:",xN,"\n")
       if (is.null(xN)) {
          print(x)
       } else {
          xNMax <- max(nchar(xN))
          xSp   <- paste0(rep(" ",xNMax),collapse="")
          for (N in xN) {
              wN <- str_sub(paste0(N,xSp),1,xNMax)
              wM <- paste0(wN,":",paste0(x[N],collapse=", "))
              cat(wM,"\n")
          }
          cat("\n")     
       }
    }
 }
 
####
#
# Set Version String for build...
#

SSMVersion <- "SaTScanMapper V1.0.1 - built 1806030 at 06:08am"
   
#
# function to get version information
#

satscanMapper.Version <- function() {
    return(SSMVersion)
  }
#
#

#### As function or direct R code...

#####
#
# Call Variable Defaults
#
#####
#
#  default values
#
# resDir    <- NULL
# prmFile   <- NULL  
# outDir    <- NULL
# censusYear<- NULL
#
# categ     <- 5
# title     <- NULL
# report    <- TRUE
# runId     <- NULL
# pValue    <- 0.05
# bndyCol   <- "grey50"  
# label     <- TRUE
# outline   <- TRUE
# proj4     <- NULL
# locO_EMap <- TRUE 
# clusO_EMap<- TRUE
# palColors <- NULL
#
#
#####
#
#  Test parameters when not running as a function
#
#####  Two Circular
#
# resDir    <- "c://projects//statnet//SatScan-R//NCIData//" 
# prmFile   <- "Two circular 1year"
# outDir    <- NULL
# categ     <- 5
# title     <- ""

# report    <- TRUE
# runId     <- "01"
# pValue    <- 0.05
# bndyCol   <- "grey50"
# label     <- TRUE
# outline   <- TRUE
# proj4     <- NULL
# locO_EMap <- TRUE
# clusO_EMap<- TRUE
# palColors <- NULL
#

#####  ellipse
#
# resDir    <- "c://projects//statnet//SatScan-R//NCIData//" 
# prmFile   <- "ellipse" 
# outDir    <- NULL
# report    <- TRUE
# runId     <- "01"
# pValue    <- 0.05
#
# categ     <- 5
# bndyCol   <- "grey50"
# label     <- TRUE
# outline   <- TRUE
# proj4     <- NULL
# locO_EMap <- TRUE
# clusO_EMap<- TRUE
# palColors <- NULL
#

######


satscanMapper <- function (resDir       = NULL,           # Required
                           prmFile      = NULL,           # Required
                           
                           outDir       = NULL,           # Optional (default = NULL)
                          
                           censusYear   = NULL,           # Optional (default = "2000")
                           categ        = NULL,           # Optional (default = 7) 
                           title        = NULL,           # Optional (default = .prm filename)
                           pValue       = NULL,           # Optional (default = 0.05()
                           report       = NULL,           # Optional (default = TRUE) 
                           runId        = NULL,           # Optional (default = "")
                           bndyCol      = NULL,           # Optional (default = "grey50")
                           label        = NULL,           # Optional (default = TRUE)
                           outline      = NULL,           # Optional (default = TRUE)
                           locO_EMap    = NULL,           # Optional (default = TRUE)
                           clusO_EMap   = NULL            # Optional (default = TRUE)
                        )
{

   ####
   #
   #  Call arguments:
   #
   #  File Directories and Names
   #
   #    resDir = Name of directory ending in "//" the contains the .prm saved 
   #           session file.
   #           Any reports will also be written to this directory. 
   #           Required and must be a character string
   #
   #    prmFile = the filename of the .prm file.  If .prm is ommited, 
   #           it is appended to the filename.
   #           Required and must be a character string.
   #
   #    outDir = the name of the directory to write the PDF graphics map 
   #           file and the TXT summary file. Parameter is optional.  If not 
   #           specified, the resDir value is used.  If outDir="", the current 
   #           work directory is used. Otherwise the outDir value must 
   #           be an valid/existing directory.
   #
   #    pValue = is a numerical value to be used for the significant test.  
   #           It must range from 0.001 to 1.  The default value is 0.05.  
   #           This value is used to determine which clusters identified by 
   #           SaTScan (TM) should be mapped.  Only cluster with pValues <
   #           this parameter are mapped.
   #
   #    report = is a logical value.  If TRUE, the package generates a text 
   #           report of the calculations and numerical information from 
   #           SaTScan (TM) related to the clusters and locations.
   #           If FALSE, no report is created.  The default value is TRUE.
   #
   #    runId = is a character string that is appended to the result files basename when creating 
   #           the filenames for the output PDF and report test files.  This allows the same SaTScan (TM) 
   #           data to be used in multiple mapping runs of the package without overlaying 
   #           the output PDF and report text files.  The value must be NULL, NA or a character string.
   #           The default value is NULL.  The character string should be no more than 2-4 characters.
   #
   #  Mapping Options
   #
   #    categ = Number of categories to present for the Cluster Low/High levels.  
   #           The valid range is 3 to 9.  
   #           The default is 5.
   #
   #    locO_EMap = is a logical variable.  If TRUE, the package will create 
   #           the mapping of the results using the Location Obs/Exp ratio data 
   #           for categorization and coloring of the areas. If FALSE, the 
   #           mapping of the Location Obs/Exp ratio data will not be done.
   #           Note: One of the locO_EMAP or clusO_EMap call parameters must 
   #           be set TRUE.
   #           Default is TRUE.
   #
   #    clusO_EMap = is a logical variable. If TRUE, the package will create 
   #           the mapping of the results using the Cluster Obs/Exp ratio data 
   #           for categorization and coloring of the areas.
   #           If FALSE, the mapping of the Location Obs/Exp ratio data will 
   #           not be done.
   #           Note: One of the locO_EMAP or clusO_EMap call parameters must 
   #           be set TRUE.  Default is TRUE.
   #
   #    title = is a character vector of 1 or 2 element to be used as the title 
   #           for the Cluster Maps.  If value is "NULL" or empty (""), 
   #           the prmFile file name will be will be used as the title.   
   #           At this time, the second title position is used to identify 
   #           the legend or specific map charateristics.
   #
   #    label = is a logical variable. If set to TRUE, the package will display 
   #           cluster number next to the cluster outline (circle or ellipse) 
   #           on the each map.  The default value is TRUE.
   #
   #    outline = is a logical variable. If set to TRUE, the package display 
   #           an circle or ellipse around the center of each cluster on the maps.
   #           The default value is TRUE.
   #
   #    proj4  = specifies the proj4 projection parameters to be used when 
   #           the output maps are drawn.  
   #           If set to NULL or "", no transformation is made to the boundary 
   #           data when the maps are drawn.  If proj4 parameter is a "string", 
   #           it must be a set of valid proj4 parameters.  The proj4 parameter 
   #           is passed to SeerMapper to transform the projection of the drawn maps. 
   #           By default the map projection is:
   #           CRS("+proj=aea +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=96w +units=m")
   #
   #    bndyCol = is a character vector reprecenting a color.  The border 
   #           color to use when mapping area boundaries in Cluster Maps.  
   #           This color is not used when boundaries are drawn for higher 
   #           level areas. The value must be a valid  color as listed by 
   #           the R "colors()" function.  
   #           The default is "grey50", the global value is set via the 
   #           internal variable CTBorder.  To disable the mapping of the 
   #           borders, set bndyCol option to "NA".
   #
   #    censusYear = is a character vector representing the census year of 
   #           the location IDs.  This is important to make sure the 
   #           boundary data matches the location FIPS codes used in the 
   #           cluster analysis and the mapping.  This package supports
   #           the census years of "2000" and "2010".  
   #	       The default value is "2000".
   #
   #    palColors = is the name from RColorBrewer of the color palette 
   #           for the mapping.  This replaces the default value of "-RdYlBu".  
   #           The use of the "-" at the beginning of the name requests 
   #           the color palette be reversed.  See RColorBrewer for the 
   #           valid list of palette names.
   #    
   #
   #  Need to Add:  ReportTier - 1 to 4, Options for passing to SeerMapper.
   #
   ####

   #####
   #
   #   Local Functions
   #

   ####
   #
   #   Validate location IDs and set the Loc type.  (2, 5, 11 and now 3 HSA)
   #
   #   Upgraded to handle HSA location ID determination and conflict with state FIPS IDs.
   #
   CleanLocID <- function(xTable) {   
      # assume it's FIPS - 2, 5, or 11 digits (goal)

      ErrFnd         <- FALSE
      xTable$LOC_ID  <- str_trim(as.character(xTable$LOC_ID))   # convert from factor to character
      
      #   State should be 1 or 2 digits
      #   County should be 4 or 5 digits
      #   Tract can be 1 or 2 + 3 + 4 or 6 -> 8 to 11 (pad front(1) and back(2))  
      #   Valid patterns are  1, 2, 4, 5, 8, 9, 10, 11.   Bad patterns are  3, 6, 7
      xIDType        <- 0
      Len_ID         <- nchar(xTable$LOC_ID)
      idList         <- xTable$LOC_ID
      
      Len_ID         <- nchar(idList)
      xM             <- (Len_ID == 1 | Len_ID == 4 | Len_ID == 8 | Len_ID == 10)
      idList[xM]     <- paste0("0",idList[xM])
      Len_ID         <- nchar(idList)
      xM             <- (Len_ID == 9)
      idList[xM]     <- paste0(idList[xM],"00")
      Len_ID         <- nchar(idList)
      
      #
      #  State:     S           (1)     or SS          (2) 
      #  County:    SCCC        (4)     or SSCCC       (5)
      #  Tract:     SCCCTTTT,   (8)     or SSCCCTTTT   (9)  
      #             SCCCTTTTTT, (10)    or SSCCCTTTTTT (11)
      #
      inVal_ID       <-(Len_ID <= 0 | Len_ID == 3 | Len_ID == 6 | Len_ID == 7 | Len_ID > 11)
      
      if (any(inVal_ID)) {
         # range of LOC_ID lengths must match the IDType (2, 5, 11) or one less.
         xmsg <- paste0("LOC_ID values in results are not in FIPS code format.\n",
                      "They must be numeric and have lengths of 1-2 digits for states FIPS, ",
                      "4-5 digits for county FIPS,\n",
                      "or 8-11 digits for census tract FIPS.  The lengths found range from ",locMin," to ",locMax,".\n",
                      "Correct and rerun package.\n")
         stop(xmsg,call.=FALSE)
         ErrFnd      <- TRUE
      } 
      xIDRange       <- range(Len_ID)
      xIDType        <- xIDRange[2]
      locMin         <- xIDRange[1]
      locMax         <- xIDRange[2]
            
      #cat("xIDRange:",xIDRange,"  xIDType:",xIDType,"\n")
  
      xTable$LOC_ID  <- idList
      #cat("xTable$LOC_ID:\n")
      #print(head(xTable$LOC_ID,10))
      
      xTable$stID    <- str_sub(idList,1,2)
        
      xLocList       <- sort(unique(idList))
      xStateList     <- sort(unique(xTable$stID))
      
      return (list(Table=xTable,IDType=xIDType,LocList=xLocList,StateList=xStateList,Err=ErrFnd))
   }
   
   
   
   #   End of CleanLocID
   #
   #####
   
   #####
   #
   #   fixDir function looks at a directory path string and adds a slash at the end
   #    if needed.
   #
   fixDir <- function(xDir) {   
      if (str_sub(xDir,-1) != "\\" && str_sub(xDir,-1) != "/")  {
          # xDir does not end with "\\" or "/", add one.
          xDir <- paste0(xDir, "/")
      }
      return(xDir)
   }
   #
   #   end of fixDir
   #
   #####
   
   #
   #   End of functions
   #
   ########

   ####
   #
   #
   #  A Few Global Variables:
   
   callVarList   <- NULL
   
   #cat("Default path:",getwd(),"\n")
   
   #
   #  Temp - Variables
   #
   proj4     <- NULL
   CRSproj4  <- NULL
   palColors <- NULL
      
   #
   #  bndyCol Colors
   #
   
   ColorB_State  <- "black"
   ColorB_County <- "grey40"
   ColorB_Tract  <- "grey50"
   LocbndyCol    <- ColorB_Tract
   
   #  Cluster Label Colors
   ClusLabCol    <- "black"
   
   ######
   #
   #  Verification Step 1 - Check ranges of subroutine call variables
   #
   FNDError     <- FALSE      # indicator a fault error was found and execution can't continue
   FNDWarn      <- FALSE      # indicator a warning message issues, but execution can continue
   NumErrors    <- 0          # number of errors
   NumWarns     <- 0          # number of warnings
   
   RateMapping  <- FALSE      # enable Rate Mapping map generation
   NoRRFile     <- FALSE      # no RR file, work around it.
   
   DoOutput     <- TRUE       # Output control flag - TRUE do output.
   
   PrmFileName  <- NULL          # PRM file name (dir, results file, extension)
   ColFileName  <- NULL          # Cluster Info file name (dir, results file, extension)
   GisFileName  <- NULL          # Location Info file name
   RRFileName   <- NULL          # RR file name
   RateFileName <- NULL          # Rate data file name
   outDirName   <- NULL          # Output Report/File Directory
   
   L2Off        <- FALSE
   
   categMode    <- 1             # 1 = single value, 2 = breakpoint list.
   
   #####
   # 
   #  Set up project 4 strings  -  Original and Projected
   #
   
   #OrigCRS <- CRS("+proj=longlat +datum=NAD83")
   
   #
   #  Transform the State, State/County, State/County/Census Tract
   #  boundary polygons from long/lat to Equidistance Conic projection.
   #
   #   Projection = Equidistance-Conical => simpleconic
   #   Lat Parallel 1   = 33
   #   Lat Parallel 2   = 45
   #   Origin of Lat    = 39
   #   central Meridian = -96   (96W)
   #   
   #  ERROR below - 
   #ProjCRS <- CRS("+proj=eqdc +lat_1=33 +lat_2=49 +lat_0=39 +lon_0=96w +units=m")
   #  Corrected string
   #ProjCRS <- CRS("+proj=eqdc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=96w +units=m")
   
   #  projections no longer used in this package - all transformations done prior
   #  when dataset are created and stored.
   #
   #####
   
   #####
   #
   #  Special Feature in satscanMapper to be able to support examples included in package.
   #  It's very difficult to set up result files from the SaTScan (TM) package that can 
   #  be distributed as examples in this package.   Since the package relies on the filenames and paths
   #  in the saved session file (.prm) and it's file references, the file paths may be difference
   #  on each installed platform.  The data supporting the example would therefore, be wrong.
   #  To compensate, the package will provide generic path substitution based on the installed
   #  R version and platform.  The following symbolic strings must be at the beginning of 
   #  a path/filename and be followed by a "/":  
   #                %HOME%/    ->   getwd() path at the start of the package.
   #                %EXTDATA%/ ->   system.file(package="satscanMapper",extdata)
   #                               This is the path to the installed library directory for satscanMapper
   #                               and it's extpath subdirectory contining the example data files.
   #                %RESDIR%/  ->   The path value provided in the resDir= call parameter.
   #  
   #  If an output file does not have a pathname specified, the %HOME%/ path will be used.
   #
   #  Writing files to the %EXTDATA%/ directory or any subdirectory is strongly not recommended,
   #
   #  These can be used in the resDir string
   #
   #
   #####
   #
   # Set up the values for the symbolics.
   #
   
   vHOMEDir    <- paste0(getwd(),"/")
   vEXTDATADir <- paste0(system.file(package="satscanMapper","extdata"),"/")
   vRESDIR     <- paste0(getwd(),"/")
   vOUTDIR     <- vRESDIR
   
   #####
   #
   #  Process and validate call parameters
   #
   
   #
   #  resDir = 
   #
   if (is.null(resDir) || is.na(resDir) || !is.character(resDir) || resDir == "") {
       FNDError   <- TRUE
       NumErrors  <- NumErrors + 1
       xmsg       <- paste0("The required resDir= call argument has not been provided.")
       warning(xmsg,call.=FALSE)
   } else {
       resDir <- as.character(resDir[[1]][1])
   }
   #
   #  prmFile =
   #
   if (is.null(prmFile) || is.na(prmFile) || !is.character(prmFile) || prmFile == "") {
       FNDError   <- TRUE
       NumErrors  <- NumErrors + 1
       xmsg       <- paste0("The required prmFile= parameter has not been provided.")
       warning(xmsg,call.=FALSE)
   } else {
       prmFile <- as.character(prmFile[[1]][1])
   }
   
   #
   #  outDir = 
   #
   if (is.null(outDir) ) {  # outDir is NULL 
      outDir   <- resDir    # use resDir as the default value.
   }
   if (is.na(outDir) )   {  # outDir is NA - no output - special case
      outDir   <- NA        # special case.
      DoOutput <- FALSE
   } else {
      outDir   <- as.character(outDir[[1]][1])
      # must be a character string
      outDir   <- str_trim(outDir)    # trim leading and trailing blanks.
      if ( outDir == "" ) {     # outDir is an empty string "" - set to the current getwd() value.
         # if value is ""   use current working directory
         outDir <- getwd()
      }
   }
   #  All values are now character strings.  
   #     Other structures and factors have been eliminated.
   
   #  Clean up path and file name strings - replace "\\" with "/"  
   #             (Windows to Unix style format)
   
   # determine if resDir ends with "/" or "\\".  If it doesn't add "/".
   #   Note "\\" turns into a single "\", but shows up as one character "\\"
 
   resDir  <- gsub("[\\]","/",resDir)
   if (str_sub(resDir,-1) != "\\" && str_sub(resDir,-1) != "/")  {
       # resDir does not end with "\\" or "/", add one.
       resDir <- paste0(resDir, "/")
   }
   # process path substitutions variables (%HOME%/ and %EXTDATA%/)
   resDir  <- sub("^%HOME%/",   vRESDIR,    resDir)
   resDir  <- sub("^%EXTDATA%/",vEXTDATADir,resDir)

   vRESDIR <- resDir
    
   # prmFile checks 
   # determine if .prm is at the end of the prmFile parameter
   orgPrmFile    <- prmFile
   baseDataFN    <-  sub("[.][^.]*$", "", prmFile,perl=TRUE)
  
   prmFile <- gsub("[\\]","/",prmFile)
   
   if (toupper(str_sub(prmFile,-4,-1)) != ".PRM") {
       # not ending in .prm - add it.
       prmFile <- paste0(prmFile,".prm")
       #xmsg <- paste0("The prmFile parameter does not end with the extension of .prm.  The extension of .prm has been added.")
       #warning(xmsg,call.=FALSE)
   }
   
   # Setup to read and validate .PRM file
   wPrmFileName <- paste0(resDir,prmFile)    # append resDir
   
   if (!dir.exists(resDir)) {
       # resDir option directory does not exist
       FNDError    <- TRUE
       NumErrors   <- NumErrors + 1
       xmsg        <- paste0("The directory specified by the resDir argument does not exist :",resDir," ")
       warning(xmsg,call.=FALSE)
   } else {
       # directory exists - does the .prm file (fullname)?
       if (!file.exists(wPrmFileName)) {
          # file does not exist in directory provided
          FNDError    <- TRUE
          NumErrors   <- NumErrors + 1
          xmsg        <- paste0("The .prm file specified by the resDir and prmFile call argument does not exist: ",prmFile)
          warning(xmsg,call.=FALSE)
          xmsg        <- paste0("Make sure the session parameters were SAVED before executing the SaTScan (TM) analysis.")
          warning(xmsg,call.=FALSE)
       }
   }

   if (DoOutput) {
      outDir  <- gsub("[\\]","/",outDir)
      if (str_sub(outDir,-1) != "\\" && str_sub(outDir,-1) != "/")  {
          # outDir does not end with "\\" or "/", add one.
          outDir <- paste0(outDir, "/")
      }
      # process path substitutions variables (%HOME%/ and %EXTDATA%/)
      outDir  <- sub("^%HOME%/",   vRESDIR,    outDir)
      outDir  <- sub("^%EXTDATA%/",vEXTDATADir,outDir)
      if (!dir.exists(outDir)) {
          # resDir option directory does not exist
          FNDError    <- TRUE
          NumErrors   <- NumErrors + 1
          xmsg        <- paste0("The directory specified by the outDir argument does not exist :",resDir," ")
         warning(xmsg,call.=FALSE)
      }
   }
   vOUTDIR <- outDir
   
   # if the resDir and prmFile and outDir arguments are not provide or are wrong, 
   #   we have to stop
   
   if (FNDError) {
      xmsg <- paste0("Call argument errors found - execution stopped.")
      stop(xmsg,call.=FALSE)
   }
   
   # no problems continue   
   
   #  Get full path to .prm file.
   
   #
   #  xRF  <- paste(prmFile,"*.*",sep="")  # form search wildcard
   #   xresDir <- str_sub(resDir,1,nchar(resDir)-1)
   #   xx <- list.files(resDir, pattern=xRF, ignore.case=TRUE)  # check master directory and Results File argument for files
   #   if (length(xx) == 0)
   #     {
   #       xmsg = paste("There are no files with Results file base name of ",ResultsFile," in \n", 
   #               "the directory - ",resDir,
   #               ".  Verify the resDir= and ResultsFile= arguments are correct.\n",
   #               sep="")
   #       message(xmsg)
   #       FNDError  <- TRUE
   #       NumErrors <- NumErrors + 1
   #       # if the resDir and ResultsFile arguments is not valid - can't continue to check the existance of  
   #       #   individual files
   #       
   #     } else {
   #       #  resDir and ResultsFile are good and reference several files. See if the right ones are there.
   #     }
   
   # resDir is the base directory for the results
   # baseDataFN is the base filename from the .prm file
   # prmFilel is the full filename and ext for the .prm
   # outDir is the base directory for writing the output.

   callVarList$resDir       <- resDir
   callVarList$prmFile      <- prmFile
   callVarList$outDir       <- outDir
   callVarList$wPrmFileName <- wPrmFileName
   callVarList$baseDataFN   <- baseDataFN
   
   vRESDIR                  <- resDir
   vOUTDIR                  <- outDir
   
   #cat("results Dir   :",resDir,"\n")
   #cat(".prm File     :",prmFile,"\n")
   #cat(".prm path&name:",wPrmFileName,"\n")
   #cat("base filename :",baseDataFN,"\n")
   #cat("output Dir    :",outDir,"\n")
   
   #
   #  end of resDir and prmFile validation checks
   #
   
   #### General Call Arguments
   
   #
   #  censusYear - "2000" or "2010"
   #
   censusYear_def  <- "2000"
   
   if (is.null(censusYear) || is.na(censusYear)) {
      # caller did not provide the censusYear call parameter
      censusYear <- censusYear_def
   
   } else {
      censusYear <- as.character(censusYear)   # make character vector (can take numerics.
      censusYear <- censusYear[[1]][1]
      
      if (censusYear != "2000" && censusYear != "2010") {
         # invalid census year value.
         xmsg <- paste0("The censusYear call parameter value is invalid: ",censusYear,". It must be 2000 or 2010. The default of 2000 will be used.")
         warning(xmsg,call.=FALSE)
         censusYear <- "2000"
      }
   }
   callVarList$censusYear <- censusYear
   #cat("censusYear:",censusYear,"\n")
   
   
   #
   #  runId
   #
   runId_Def <- ""
   
   if (is.null(runId)) {
      # runId not present - set to default
      runId <- runId_Def
   } else {
      if ( is.na(runId) || !is.character(runId) || length(runId) <= 0) {
         # runId argument is NA or not correctly formated. Set to default value of "".
         xmsg        <- paste0('The runId call argument is set to NA, not a character string or empty.  runId set to the default of "".')
         NumErrors   <- NumErrors + 1
         warning(xmsg,call.=FALSE)
         runId       <- runId_Def
      } else {
         if (nchar(runId) > 8) {
            xmsg     <- paste0("The length of the runId string is greater than 8.  Please reduce the number of characters.")
            warning(xmsg,call.=FALSE)
         }
      }
   }
   
   callVarList$runId <- runId
   
   #
   # pValue
   #
   pValue_Def <- 0.05
   
   if (is.null(pValue) || is.na(pValue)) {
      # pValue is not provided - set to default
      pValue <- pValue_Def
   
   } else {
      if (!is.numeric(pValue) || pValue < 0.005 || pValue > 0.5) {
         # the pValue argument is set to NA, a non-numeric or out of range.
         xmsg     <- paste0("pValue call argument specified is not a numeric value or not within range of (0.005 to 0.5).\n",
                 "The default value of 0.05 will be used.\n")
         warning(xmsg,call.=FALSE)
         NumWarns <- NumWarns + 1
         pValue   <- pValue_Def
     }
   }
   callVarList$pValue <- pValue
   
   #
   #### Mapping Arguments
   #

   #
   # title   
   #
   vTitle <- ""
   if (is.null(title)) {
      # no title argument provided
      vTitle      <- orgPrmFile
   
   } else {
      if (is.na(title) || !is.character(title) || nchar(title) <= 0 || length(title) <= 0) {
         # title is set to NA, not a character vector or empty.
         xmsg     <- paste0("The title call argument is NA, a non-character string, a non-vector, or empty.\n",
                "The title is set to the prmFile name of ",prmFile,"\n")
         warning(xmsg,call.=FALSE)
         NumWarns <- NumWarns + 1
         vTitle   <- orgPrmFile
      } else {
         if (!is.vector(title)) {
            xmsg <- paste0("The title call argument is must be a vector of strings for multiple line titles.\n",
                  "The title is set to the prmFile name of ",prmFile,"\n")
                
            warning(xmsg,call.=FALSE)
            NumWarns <- NumWarns + 1
            vTitle   <- orgPrmFile
         } else {
           title <- as.character(title)
           if (length(title) > 2) {
              title <- title[1:2]   # only take first two
           } else {
              title <- title[[1]][1]
           }
         }
         vTitle <- title
      }
   }
   callVarList$title <- vTitle

   #
   # label   
   #
   label_def <- TRUE
   vLabel    <- label_def
   
   if (is.null(label) || is.na(label) ) {
      # no label argument provided or set to NA
      vLabel <- label_def
   
   } else {
      if ( !is.logical(label) ) {
         # label must be TRUE or FALSE logical
         xmsg <- paste0("The label call argument is a non-logical variable.\n",
                "The label option is set to the default value of TRUE.\n")
         warning(xmsg,call.=FALSE)
         NumWarns = NumWarns + 1
         vLabel    <- label_def
     } else {
         vLabel    <- label
     }   
   }
   callVarList$label <- vLabel
   label             <- vLabel

   #
   # bndyCol 
   #
   bndyCol_def <- ColorB_Tract
   
   if (is.null(bndyCol)) {
      # no argument provided.
      bndyCol <- bndyCol_def
   
   } else {
      # bndyCol can be either TRUE/FALSE (on or off) or a color.
   
      if (is.na(bndyCol)) {
         # bndyCol set to NA - turn it off.
         bndyCol <- ""
      } else {
         if (is.logical(bndyCol)) {
            # if logical = TRUE or FALSE
            if (bndyCol) {
               # TRUE - enable
               bndyCol  <- bndyCol_def
            } else {
               # FALSE - disable
               bndyCol  <- ""
            }
         } else {
         
            # not logical, must be character => color name
            if (is.character(bndyCol)) {
               # if character = bndyCol color.     
               if (!any(bndyCol == colors())) {
                  # bndyCol color does not match a valid name
                  xmsg      <- paste0("The bndyCol call argument of ",bndyCol," is not a valid color name or empty. The default value of ",ColorB_Tract," is used. \n")
                  warning(xmsg,call.=FALSE)
                  FNDWarn   <- TRUE
                  NumWarns  <- NumWarns + 1
                  bndyCol   <- bndyCol_def 
               }
            } else {
               # bndyCol is not character vector
               xmsg     <- paste0("The bndyCol call argument is not valid data type.  Boundaries will be drawn using the default color of ",ColorB_Tract,"\n")
               warning(xmsg,call.=FALSE)
               FNDWarn  <- TRUE
               NumWarns <- NumWarns + 1
               bndyCol  <- bndyCol_def
            }
         }
      }
   }
   #cat("bndyCol:",bndyCol,"\n")
   callVarList$bndyCol <- bndyCol
   
   #
   # categ
   #
   categ_Def <- 5 
   
   if (is.null(categ)) {
      # no call argument provides
      categ   <- categ_Def   # set to the default.
   } else {
      categ <- as.integer(categ[[1]][1])
      if ( is.na(categ) || categ < 3 || categ > 9) {
        # categ call argument is set to NA or out of range.
        xmsg     <- paste0("The categ call argument is set to NA or is out of range ( < 3 or > 9). \n",
                           "The default value of ",categ_Def," will be used.\n")
        warning(xmsg,call.=FALSE)
        FNDWarn  <- TRUE
        NumWarns <- NumWarns + 1
        categ    <- categ_Def
     }
   }
   callVarList$categ <- categ
   
   #
   #   ####>>>>  Update to handle break point list specification.
   #
   
   #
   # locO_EMap
   #
   locO_EMap_def <- TRUE
   
   if (is.null(locO_EMap) || is.na(locO_EMap)) {
      # no argument provides or it is set to NA
      locO_EMap <- locO_EMap_def
   } else {
      if (!is.logical(locO_EMap)) {
         # locO_EMap  is not a logical variable
         xmsg        <- paste0("The locO_EMap call argument is not a logical variable.  The default value of ",locO_EMap_def," will be used.\n")
         warning(xmsg,call.=FALSE)
         FNDWarn     <- TRUE
         NumWarns    <- NumWarns + 1
         locO_EMap   <- locO_EMap_def
      } else {
         locO_EMap   <-  locO_EMap[[1]][1]  # get first value
      }
   }
   
   #
   # clusO_EMap
   #
   clusO_EMap_Def <- TRUE
   
   if (is.null(clusO_EMap) || is.na(clusO_EMap)) {
      # no argument provides or it is set to NA
      clusO_EMap <- clusO_EMap_Def
   } else {
      if (!is.logical(clusO_EMap)) {
         # clusO_EMap  is not a logical variable
         xmsg        <- paste0("The clusO_EMap call argument is not a logical variable.  The default value of ",clusO_EMap_Def," will be used.\n")
         warning(xmsg,call.=FALSE)
         FNDWarn     <- TRUE
         NumWarns    <- NumWarns + 1
         clusO_EMap  <- clusO_EMap_Def
      } else {
         clusO_EMap  <-  clusO_EMap[[1]][1]  # get first value
      }
   }
   
   # Check to make sure one was set.
   
   #if (!locO_EMap && !clusO_EMap) {
   #   xmsg    <-  paste0("Both locO_EMap and clusO_EMap have been set to FALSE.  No map will be drawn in this situation.  Both parameters have been se to TRUE.")
   #   warning(xmsg, call.=FALSE)
   #   locO_EMap  <- TRUE
   #   clusO_EMap <- TRUE
   #}
   
   #
   #
   #  ????  Do we let loc and clus both be SET ????
   #
   #
   
   callVarList$locO_EMap  <- locO_EMap
   callVarList$clusO_EMap <- clusO_EMap
   
   # 
   # outline
   #
   outline_Def <- TRUE
   
   if (is.null(outline)) {
      # no argument provided
      outline <- outline_Def
      
   } else {
      if (is.na(outline) || !is.logical(outline)) {
         xmsg     <- paste0("The outline call argument is not a logical variable. The default of ",outline_Def," will be used.\n")
         message(xmsg,call.=FALSE)
         FNDWarn  <- TRUE
         NumWarns <- NumWarns + 1
      }
   }
   callVarList$outline <- outline
   
   #
   #  Future call parameters:
   #
   #     palColors
   #     stateB
   #     countyB
   
   #
   # End of call argument checks
   #
   #####
 
   #cat("Call parameters - checked.\n") 
   #printCallVarList(callVarList)
   
   
   ###
   #
   #  At the same time as checking the arguments, their defaults are set, if there 
   #    is an error or warning called by the user.
   #
   #  Arguments are check as best we can to this point.  Now to try and use the information.
   #
   ###

   if (FNDWarn) {  
      # warnings found
      xmsg <- paste0("\n\nThe following number of warnings were found checking the call arguments:",NumWarns,".\n",
             "     The package may not run as expected. \n",
             )
      warning(xmsg,call.=FALSE)
   }
   
   if (FNDError) {  
      # terminal errors found
      xmsg <- paste0("\n\nThe following number of terminal errors were found checking the call arguments:",
                     NumErrors,".\n",
                     "  Run will be terminated.  Please fix errors and re-run.\n")
      stop(xmsg,call.=FALSE)
   }

   #
   # End resDir, prmFile, outDir and call argument processing and verification.
   #
   ###########################################################

   ########################
   #
   #   Parameter File Processing and Verification
   #
   #
   #  Variable/Tables used to analyze PRM file.
   #
   
   
   PrmChkN <- c(               # PRM file parameter (field) names
      "CaseFile"                    # nu Case file path/filename
     ,"PopulationFile"              # nu Population file path/filename
     ,"CoordinatesFile"             # nu Coordinates file path/filename
     
     ,"PrecisionCaseTimes"          # 0=None, 1=Year, 2=Month, 3=Day, 4=Generic (1)
     ,"CoordinatesType"             # 0=Cartesian 1=Lat/Long  (0, 1)
     ,"StartDate"                   # a date (-1)              
     ,"EndDate"                     # a date (-1)
     ,"AnalysisType"                # 1=Purely Spatial, 
                                    # nu 2=Purely Temporal(No) 
                                    # 3=Retrospective ST 
                                    # 4=Prospective ST, 
                                    # nu 5=Spatial Var. in Temporal Trends (No), 
                                    # nu 6=Prospective Purely Temporal (No)
                                    #   1, 3, and 4 are acceptable  
  
     ,"ModelType"                   # 0=Discrete Poisson, 
                                    # 1=Bernoulli, 
                                    # 2=Space-Time Permutation, 
                                    # nu 3=Ordinal, 
                                    # 4=Exponential, 
                                    # nu 5=Normal, 
                                    # 6=Continuous Poisson, 
                                    # nu 7=Multinomial
                                    
     ,"TimeAggregationUnits"        # 0=None, 1=Year, 2=Month, 3=Day, 4=Generic  (1)
     ,"TimeAggregationLength"       # (Positive Integer) in above units. (-2)  (0 to n)
     ,"ResultsFile"                 # Base Filename - string (-3)  (text file)
     ,"IncludeRelativeRisksCensusAreasDBase" # Must be "y" (RR)
     ,"CensusAreasReportedClustersDBase"     # Must be "y" (Cluster)
     ,"MostLikelyClusterEachCentroidDBase"   # Must be "y" (location)
     ,"MostLikelyClusterCaseInfoEachCentroidDBase"   # nu

     ,"SpatialWindowShapeType"      # 0=Circular, 1=Elliptic   - either are acceptable
     ,"MaxTemporalSize"             # Max temporal size of cluster over time, positive integer (-2) (???) 
     ,"Version"                     # Versions - match on first character in case of sub-releases.                  
   )
   
   #
   # PrmChk values => 0   value
   #               = -1   don't care - data
   #               = -2   positive integer
   #               = -3   string value
   #               = -4   SaTScan Version number
   #               = -5   string path/filename with symbolic substitutions
   #               = -6   string path/filename with symbolic substitutions (OPTIONAL)
   #
   
   PrmChkV <- list(   # Acceptable Values
       c(-6)               # Case File (path/filename)
      ,c(-6)               # Population File (path/filename)
      ,c(-6)               # Coordnates File (path/filename)
      ,c(1)                # PrecisionCaseTime 0=None, 1=Year, 2=Month, 3=Day, 4=Generic
      ,c(0,1)              # CoordinatesType   0=Cartesian 1=Lat/Long 
      ,c(-1)               # StartDate         a date (-1)              
      ,c(-1)               # EndDate           a date (-1)
      ,c(1,3,4)            # AnalysisType      1=Purely Spatial, 2=Purely Temporal(No) 3=Retrospective ST 4=Prospective ST, 5=Spatial Var. in Temporal Trends (No), 6=Prospective Purely Temporal (No)
      ,c(0,1,2,4,6)        # ModelType         0=Discrete Poisson, 1=Bernoulli, 2=Space-Time Perm, 3=Ordinal, 4=Exponential, 5=Normal, 6=Continuous Poisson, 7=Multinomial
      ,c(0,1,4)            # TimeAggregationUnits  0=None, 1=Year, 2=Month, 3=Day, 4=Generic
      ,c(-2)               # TimeAggregationLength Positive Integer) in above units. (-2)
      ,c(-5)               # ResultsFile       Base Filename - path/filename (-5)
      ,c("y")              # IncludeRRDBase    Must be "y" 
      ,c("y")              # CensusAreaRDBase  Must be "y"
      ,c("y")              # MostLikeClusDBase Must be "y"
      ,c("y")              # MostLikeClusCentroidDBase Must be "y"
      ,c(0,1)              # SpatialWindowShapeType 0=Circular, 1=Elliptic) 
      ,c(-2)               # MaxTemporalSize   Max temporal size of cluster over time, positive integer (-2) 
      ,c(-4)               # Version           Versions - match on first character in case of sub-releases.                 
   )
   
   ###
   #
   # Aceptable values of -1, -2, -3 indicate a type of value is acceptable.
   #    -1 = any date
   #    -2 = any positive integer
   #    -3 = any character string
   #    -4 = SaTScan Version number
   #    -5 = character string - path/filename with substitutes.
   #
   ###
    
   PrmChkT <- c(         # Target variable name within package
      "CaseFile"                 # CaseFile
     ,"PopFile"                  # PopFile
     ,"GeoFile"                  # GeoFile
     ,"PrecisionCaseTimes"       # 0=None, 1=Year, 2=Month, 3=Day, 4=Generic
     ,"CoordinatesType"          # 0=Cartesian 1=Lat/Long 
     ,"Prm.StartDate"            # a date (-1)              
     ,"Prm.EndDate"              # a date (-1)
     ,"AnalysisType"             # 1=Purely Spatial, 
                                 # 2=Purely Temporal(No), 
                                 # 3=Retrospective Space-Time,
                                 # 4=Prospective Space-Time, 
                                 # 5=Spatial Var. in Temporal Trends (No),
                                 # 6=Prospective Purely Temporal (No)
     ,"ModelType"                #  Acceptable: 
                                 #   0=Discrete Poisson
                                 #   1=Bernoulli
                                 #   2=Space-Time Permutation, 
                                 #   3=Ordinal (NO), 
                                 #   4=Exponential, 
                                 #   5=Normal (NO), 
                                 #   6=Continuous Poisson, 
                                 #   7=Multinomial (NO)
     ,"TimeAggUnits"             # 0=None, 1=Year, 2=Month, 3=Day, 4=Generic
     ,"TimeAggLen"               # (Positive Integer) in above units. (-2)
     ,"baseFile"                 # Base Filename - string (-3)
     ,"RRFileDBase"              # Must be "y" 
     ,"ColFileDBase"             # Must be "y"
     ,"GisFileDBase"             # Must be "y"
     ,"ColCenFileDBase"          # Must be "y"
     ,"Prm.ShapeType"            # 0=Circular, 1=Elliptic) 
     ,"MaxTemporalSize"          # Max temporal size of cluster over time, positive integer (-2) 
     ,"SatScanVer"               # Versions - match on first character in case of sub-releases.                 
     )
   
   PrmChkDef <- c(   # Variable Default Values
       # default values
        NA                      # CaseFile => NA
       ,NA                      # PopFile  => NA
       ,NA                      # GeoFile  => NA
       ,1                       # Precision Case Time   # Year
       ,0                       # Coordinates Type      # Cartisian
       ,""                      # Start Date
       ,""                      # End Date
       ,NA                      # Analysis Type - don't care  # 1 to 6
       ,NA                      # Model Type - don't care but some don't work
       ,1                       # Time Agg Unit         # Year
       ,1                       # Time Agg Len          # 1
       ,NA                      # ResultsFile           # user provided
       ,"y"                     # RR DBF File           # DBF = YES
       ,"y"                     # Col DBF File          # DBF = YES
       ,"y"                     # Gis DBF File          # DBF = YES
       ,"y"                     # ColCent DBF File      # DBF = YES
       ,NA                      # Shape Type - don't care  # 0 or 1
       ,1                       # Max Temporal Size (???)  # 1 
       ,NA                      # SatScanVersion - don't care  # version number for checking.
       )
    
   #
   #
   #####
   
   PrmDF <- data.frame(Name=PrmChkN,Target=PrmChkT,Expected=I(PrmChkV),Def=PrmChkDef)
   
   #####
   #
   #  Keep RCMD Check Happy
   #
   ColFileDBase    <- "n"        # cluster info
   GisFileDBase    <- "n"        # location info
   RRFileDBase     <- "n"        # RR
   #
   #  place holders 
   CoordinatesType <- -1
   Prm.ShapeType   <- -1
   MaxTemporalSize <- -1
   
   #######
   #####
   ###
   #
   #   Open, Read and decode prm file.
   #
   
   # Create connector and open
   PrmCon    <- file(description=wPrmFileName,open="",blocking=TRUE)
   open(con=PrmCon,open="r")
   
   #  Read PRM file lines.
   TxLns     <- readLines(con=PrmCon,n=-1,ok=TRUE,warn=TRUE)
   
   #  Read Done - close connection
   close(con=PrmCon)
   
   #  How many lines do we have?
   TxLnsLen  <- length(TxLns)
  
       
   if (TxLnsLen > 0) {     # continue processing is the file contained lines.
     
      # Delete extra blanks around text.
      TxLns     <- str_trim(TxLns, side="both")
   
      # Identify and delete items of no interest  - blank lines and lines starting with ";" and "["
      TxLnsGood <- (TxLns != ""  & str_sub(TxLns,1,1) != ";" & str_sub(TxLns,1,1) != "[")
      TxLns2    <- TxLns[TxLnsGood]
   
      TxLnsLen  <- length(TxLns2)   # get current length
       
      TxLns3    <- str_split(TxLns2,"[=,]")  # split lines on "=" and "," characters
       
      PrmLen    <- sapply(TxLns3, function(x)  length(x))   # get list of number of elements on each line.
      PrmMax    <- max(PrmLen)               # save the maximum number of elements on a line.
   
      # build column names for split results.
      PrmLab    <- paste("P",formatC(seq(from=1,to=(PrmMax-1)),format="d",flag="0",digits=0,width=2),seq="")
      PrmLab    <- str_replace_all(PrmLab," ","")
   
      # build PrmList - field: name, length, value(s) 
      PrmList   <- t(sapply(TxLns3, function(x) c(x[1], length(x)-1, x[2:PrmMax]) ))
      # at to far left the original line text.
      PrmList   <- cbind(PrmList,TxLns2)
       
      # Set the column names for processing.
      colnames(PrmList) <- c("Name","Len",PrmLab,"Orig")
   
      # 
      # At this point the PRM file parameters are parsed into a data frame "PrmList".
      # Columns:
      #    parameter name
      #    number of elements in assignment
      #    V1 to Vn for the elements in the assignment
      #
   
   } else {
       # The file is empty - no lines to scan - probable error.  How far can be go?
       xmsg  <- paste0("There is no information or data in the .PRM file.  Run Terminated.")
       stop(xmsg,call.=FALSE)
       # for now - if no file or data, we may be making a wrong default setting.
       #   How much can the data we have tell us?
   }
   
   #
   # Set up variables to process the PrmList and set the R code variables appropriately.
   #
   # Set Default Values in variables - in case they are not in the PRM file.
   #
   # may be best to set these to NULL - then check afterword for issues or omissions.
   #
   
   #  Assign the NA value to each variable
   
   for (j in seq(from=1, to=length(PrmChkN))) {
     
       assign(PrmChkT[j], NA)    # set all to NA.
   }
   
   # Get a list of matched items = PRM file to the table of interesting parameters
   
   PrmMatch  <- match(PrmList[,1],PrmChkN)
   
   PrmLen    <- length(PrmMatch)
   FNDError  <- FALSE
   NumErrors <- 0
   NumFnd    <- 0
  
   
   # Process each matched parameter.
   
   needPrmHdr    <- TRUE

    
   for (inx in c(1:PrmLen)) {               # loop through all *.Prm parameters.
      if (!is.na(PrmMatch[inx])) {        
         # match - process parameter
   
         NumFnd <- NumFnd + 1
         jx     <- PrmMatch[inx]
         
         yNAME  <- PrmList[inx,"Name"]
         yP01   <- PrmList[inx,"P01"]      # Get first value (the only one we want to see.)
         
         if (yP01 == "") yP01 <- NA      # if value = "", then NA
          
         yChkVal <- PrmChkV[[jx]][1]     # get first value to see if neg and needs different processing
  
         #cat("mode:",yChkVal,"  inx:",inx,"  jx:",jx,"  yP01:",yP01," valid:",PrmChkV[[jx]],"\n")
          
         if (!(yChkVal < 0))  {          # No "special" processing - number or string
            # Verify the value is supported.
            
            assign(PrmChkT[jx],yP01)
              
         } else {
            
            # if minus special handling
            if (yChkVal == -1) {
               # -1   Date      yyyy/mm/dd
                   
               xDate <- try( as.Date(yP01,"%Y/%m/%d" ) )
               if (is.null(xDate) || is.na(xDate) || class(xDate) == "try-error"  ) {
                  # error in validating/converting date
                  #  Error message
                  xmsg      <- paste0("The following date in the *.PRM file is not a valid DATE format: ",PrmList[inx,"Orig"])
                  warning(xmsg,call.=FALSE)
                  FNDError  <- TRUE
                  NumErrors <- NumErrors + 1
               
               } else {
                  # String to Date completed OK
                  assign(PrmChkT[jx],xDate)   # assign date to variable
               }
                   
            }
            if (yChkVal == -2) {
               # -2   Positive Integer
                  
               x <- strtoi(yP01)
               if (is.null(x) || is.na(x) || x < 0 ) {
                   #  Error message
                   xmsg      <- paste0("The following value in the *.PRM file was not a positive integer: ",PrmList[inx,"Orig"],"\n")
                   warning(xmsg,call.=FALSE)
                   FNDError  <- TRUE
                   NumErrors <- NumErrors + 1
                       
               } else {
                   assign(PrmChkT[jx],as.integer(yP01))  # assign integer to variable
               }
            }
            if (yChkVal == -3 ) {
               # -3   String (name, file, etc.)
                  
               if (is.null(yP01) || is.na(yP01) || !is.character(yP01) || length(yP01)<=0)  {
                  # Error message
                  xmsg      <- paste0("The following SaTScan (TM) parameter is empty and must be set to a value in the *.PRM file: ",PrmList[inx,"Orig"],"\n")
                  warning(xmsg,call.=FALSE)
                  FNDError  <- TRUE
                  NumErrors <- NumErrors + 1
               } else {
                  assign(PrmChkT[jx],yP01)   # assign string to variable
               }
            }
            if (yChkVal == -4) {
               # -4   Version Match - wildcard.
               if (is.null(yP01) || is.na(yP01) || !is.character(yP01) || length(yP01) <= 0 ) {
                  xmsg       <- paste0("The version parameter is empty or missing in the .PRM file.\n")
                  warning(xmsg,call.=FALSE)
                     
                  SatScanVer <- "9.1.1"
                
               } else {
                  # have good field?
                       
                  vPat <- "^9\\.[1234].*"  # 9.1, 9.2, 9.3 or 9.4 leading versions
               
                  if (!grepl(vPat,yP01))  {
                  
                     # version is not supported.
                     xmsg     <- paste0("SatScan Version ",yP01," has not been tested with this package.\nRun will be attempted, please report any problem.\n")
                     warning(xmsg,call.=FALSE)
                     NumWarns <- NumWarns + 1
                     FNDWarn  <- TRUE
                  }
               }
            }
            if (yChkVal == -5 ) {
               # -5   String (path / filename with symbolic substituions)
                  
               if (is.null(yP01) || is.na(yP01) || !is.character(yP01) || length(yP01)<=0)  {
                  # Error message
                  xmsg      <- paste0("The following SaTScan (TM) parameter is empty and must be set to a value in the *.PRM file: ",PrmList[inx,"Orig"],"\n")
                  warning(xmsg,call.=FALSE)
                  FNDError  <- TRUE
                  NumErrors <- NumErrors + 1
               } else {
                  yP01 <- gsub("[\\]","/",yP01)
                  yP01 <- sub("^%HOME%/",vHOMEDir,yP01)
                  yP01 <- sub("^%EXTDATA%/",vEXTDATADir,yP01)
                  yP01 <- sub("^%RESDIR%/",vRESDIR,yP01)
                  
                  assign(PrmChkT[jx],yP01)   # assign string to variable
               }
            }
            if (yChkVal == -6 ) {
               # -6   String (path / filename with symbolic substituions)  (Optional)
                    
               if (is.null(yP01) || is.na(yP01) || !is.character(yP01) || length(yP01)<=0)  {
                  #  If missing - no problem - ignore.
                  assign(PrmChkT[jx],"")
                  
                  ## Error message
                  #xmsg      <- paste0("The following SaTScan (TM) parameter is empty and must be set to a value in the *.PRM file: ",PrmList[inx,"Orig"],"\n")
                  #warning(xmsg,call.=FALSE)
                  #FNDError  <- TRUE
                  #NumErrors <- NumErrors + 1
               } else {
                  yP01 <- gsub("[\\]","/",yP01)
                  yP01 <- sub("^%HOME%/",vHOMEDir,yP01)
                  yP01 <- sub("^%EXTDATA%/",vEXTDATADir,yP01)
                  yP01 <- sub("^%RESDIR%/",vRESDIR,yP01)
                  
                  assign(PrmChkT[jx],yP01)   # assign string to variable
               }
            }
         }
      }  # end of check for match
   }  # loop to next entry in list.
   
   ####
   #
   #  Have all important parameters from *.PRM decoded and in memory.
   #
   #  Y and N entries are not check here.  
   #
   ####


   ######
   #
   #  Clean up some variables from .prm
   #
   #   baseFile  -> results text report   path/filename.txt
   #

   ResultsText <- baseFile
   callVarList$ResultsText <- ResultsText
   
   x <- str_sub(baseFile,-4,-1)  # get possible extension
   if (x == ".txt")   baseFile <- str_sub(baseFile,1,-5)
   
   callVarList$baseFile    <- baseFile
   
   #cat("baseFile:\n")
   #print(baseFile)

   #
   ######
   
   ######
   ##
   ##   Check PRM parameters on files - Read and setup Tables
   ##
   ##    Done after .prm file read to determine if flags set saying DBF file should exist.
   ##
   #cat("wColFileName:\n")
   wColFileName  <- paste0(baseFile,".col.dbf")    # cluster info filename
   #print(wColFileName)

   #cat("wColFile:\n")
   wColFile      <- basename(wColFileName)
   #print(wColFile)
   
   #cat("dirname(wColFileName):\n")
   #print(dirname(wColFileName))
   
   if (ColFileDBase != "y") {
      # Cluster Information File (*.col.dbf) not requested.  
      xmsg      <- paste0("The PRM file indicates the results Cluster Information DBF file (",wColFile,") was not requested.\n",
                 "Please check the Output Cluster Info. DBF box and rerun SaTScan (TM).\n")
      warning(xmsg,call.=FALSE)
      FNDError  <- TRUE
      NumErrors <- NumErrors + 1
   } else {
      # should have the Required "*.col.dbf"  cluster information
      if (!file.exists(wColFileName)) {
         # cluster info file missing
         xmsg       <- paste0("The cluster information results DBF file (",wColFile,
                        ") does not exist.\n",
                        "If you moved the results files, please copy them back to the original results directory.\n",
                        "Results Directory:",dirname(wColFileName),"\n")
         warning(xmsg,call.=FALSE)
         FNDError  <- TRUE
         NumErrors <- NumErrors + 1
      }
   }
   callVarList$ColFileName <- wColFileName
   
   
   wGisFileName   <- paste0(baseFile,".gis.dbf")   # location filename
   wGisFile       <- basename(wGisFileName)		   # path and filename
   #print(wGisFileName)
   
   if (GisFileDBase != "y") {
      # location Information File (*.gis.dbf) not requested.
      xmsg      <- paste0("The PRM file indicates the results Location Information DBF file (",wGisFile,") was not requested.\n",
                 "Please check the Output Location Info. DBF box and rerun SaTScan (TM).\n")
      warning(xmsg,call.=FALSE)
      FNDError  <- TRUE
      NumErrors <- NumErrors + 1
   } else {
      # should have the Required "*.GIS.dbf" file
      if (!file.exists(wGisFileName)) {
         # location info file missing
         xmsg       <- paste0("The location information results DBF file (",wGisFile,
                        ") does not exist.\n",
                        "If you moved the results files, please copy them back to the original results directory.\n",
                        "Results Directory:",dirname(wGisFileName),"\n")
         warning(xmsg,call.=FALSE)
         FNDError  <- TRUE
         NumErrors <- NumErrors + 1
      }
   }
   callVarList$GisFileName <- wGisFileName
   
   
   wRRFileName   <- paste0(baseFile,".rr.dbf")   # RR path and filename
   wRRFile       <- basename(wRRFileName)	 # RR filename
   #print(wRRFileName)
   
   if (RRFileDBase != "y") { 
      # relative risk data for location (*.rr.dbf) not requested.
      xmsg <- paste0("The PRM file indicates the results Relative Risk DBF file (",wRRFile,") was not requested.\n",
                 "Please check the output RR DBF box and rerun SaTScan (TM).\n")
      warning(xmsg,call.=FALSE)
      FNDError = TRUE
      NumErrors = NumErrors + 1
   } else {
      # should have the Required "*.RR.dbf" file
      if (!file.exists(wRRFileName)) {
         #  Relative Risk file missing 
         xmsg       <- paste0("The Relative Risk results DBF file (",wRRFile,
                        ") does not exist.\n",
                        "If you moved the results files, please copy them back to the original results directory.\n",
                        "Results Directory:",dirname(wGisFileName),"\n")
         warning(xmsg,call.=FALSE)         
         xmsg     <- paste0("If your type of analysis has the relative risk estimates grayed out,\n", 
                       "this package cannot identify and outline the all of the areas involved in the analysis.\n",
                       "The package will map all of the area boundaries in the states containing\n",
                       "location IDs appearing in the results.\n" )
         warning(xmsg,call.=FALSE)
         FNDWarn  <- TRUE
         NumWarns <- NumWarns + 1
         NoRRFile <- TRUE
      } else {
         NoRRFile <- FALSE
      }
   }
   callVarList$ResultsText <- baseFile
     
   #
   #   End of results file verification and locating.
   #
   ####
   
   #
   #  CoordinatesType
   #
   CoordinatesTypeLit <- ""
   
   if (CoordinatesType != 0) {
      # cartisian coordinates were not used.
      #xmsg <- paste0("Saved sesson file (.prm) indicates the coordinates used was Latitude and Longitude.\n",
      #           "Please select the cartesian coordinates and rerun SaTScan (TM) with the geo file generated by this package.\n")
      #warning(xmsg,call.=FALSE)
      #FNDError = TRUE
      #NumErrors = NumErrors + 1
      CoordinatesTypeLit <- "Lat/Long (Converted)"
   } else {
      CoordinatesTypeLit <- "Castesian"
   }
   
   #cat("CoordinatesType:",CoordinatesType,"  ",CoordinatesTypeLit,"\n")
   
   callVarList$CoordinatesType <- CoordinatesType
   callVarList$Coordinates     <- CoordinatesTypeLit
   
   #
   # The coordinates values used by the caller can be anything.  Their
   # choice does effect the cluster calculation, but not the mapping.
   # All coordinates are replaced by the packages Albers Equal Area project
   # coordinates and centroids based on the LOC_ID values. 
   #
     
   #
   #  Coordinates ASSUMPTION:
   #
   #    All X,Y (Long,Lat) coordinates provided by SatScan are 
   #    based on the original *.geo input file.  In this situation
   #    the coordinates are generated by a package like ESRI
   #    and have been transformed into a "Equidistance Conic" projection 
   #    in meters, with lat1=33, lat2=45, lat0=39 and long0=96w 
   #    (lat0 and long0 are the center point of the projection.  
   #    The projection used in the package's centroids and 
   #    boundaries are essentially the same.  However,
   #    the caller may use other projections for the cartesian coordinates
   #    of the LOC_ID areas.  
   #
   #    To attempt to be able to map any and all results data from 
   #    SatScan (TM), the package keys on the LOC_ID values and not
   #    any coordinates values.  All of the centroid values and 
   #    distance values are converted to the Albers Equal Area projection 
   #    used by the boundary data information.
   #    The centroids are directly replaced.  The circle and ellipical 
   #    dimensions are a ratio of the centroids in the results
   #    to the replacement centroids in the boundary data information.
   #
   #    For example:  if a cluster consists of only one LOC_ID, the distances
   #    are Zero.  In the future this may be adjusted to the size of the LOC_ID area.
   #    for visualization.   If a cluster has three LOC_IDs, then the sum of 
   #    there X and Y coordinates will be divided by the sum of the new X and Y 
   #    coordinates and the ratio divided into the circle or elliptical metrics.
   #
   ###
   
   #SSShape <- 0    # 0 = circular
   #                # 1 = ellipical
   
   SSShape    <- NULL
   SType      <- "Invalid"
   
   if (Prm.ShapeType == 0) { 
      SSShape  <- 0
      SType    <- "Circular"
   }
   
   if (Prm.ShapeType == 1) {
      SSShape  <- 1
      SType    <- "Elliptical"
   }
   
   if (is.null(SSShape)) {
      # No ShapeType was found in the .prm file - invalid shape type
      xmsg      <- paste0("Invalid or no Shape Type found = ",Prm.ShapeType," in session saved file (.prm).\n")
      warning(xmsg,call.=FALSE)
      FNDError  <- TRUE
      NumErrors <- NumErrors + 1
   }
     
   #cat("SSShape:",SSShape,"  >>  ",SType,"\n")
   
   callVarList$SSShape <- SSShape
   callVarList$SType   <- SType
   
   #
   #  The package is hardcoded to run either Spatial or Spatial-Time output data with start and 
   #    stop dates.  SatScan generates several formats depending on the type of analysis.  The
   #    program supports the following three forms:
   #
   #    a) Retrospective Purely Spatial
   #    b) Retrospective Space-Time
   #    c) Prospective Space-Time
   #
   #   Set via function call arguments
   #
   
   #
   #if (SSForm < 0 && SSForm > 2)
   #  {
   #    if (SSForm != "S" && SSForm != "ST" && SSForm != "SD")
   #      {
   #        message("SSForm Argument out of range. Must be 0, 1, 2, S, ST, or SD")
   #        FNDError <- TRUE
   #      }
   #  }
   #
   # Old argument definition:
   #
   #  SSForm  <- 1    # 0 = Spatial Only with start and stop dates
   #                  # 1 = Spatial-Time with start and stop dates
   #                  # 2 = Spatial Only without start and stop dates (retired)
   #
   # Parameter file definition:
   #
   #  AnalysisType => 1=Purely Spatial (no dates or with dates?),
   #                  2=Purely Temporal (Not supported) 
   #                  3=Retrospective ST 
   #                  4=Prospective ST, 
   #                  5=Spatial Var. in Temporal Trends (Not supported), 
   #                  6=Prospective Purely Temporal (Not supported)
   #
   
   SForm         <- NULL
   AnalysisType  <- as.integer(AnalysisType)
   ModelType     <- as.integer(ModelType)
   
   if (AnalysisType == 1) {
       SSForm    <- 0
       SForm     <- "Spatial Only"
   }
   if (AnalysisType == 2) {
       xmsg      <- paste0("AnalysisType of Purely Temporal is not supported.  No spatial content.\n")
       warning(xmsg,call.=FALSE)
       FNDError  <- TRUE
       NumErrors <- NumErrors + 1
   }
   if (AnalysisType == 3) { 
       SSForm    <- 1
       SForm     <- "Retrospective Space-Time"
   }
   if (AnalysisType == 4) { 
       SSForm    <- 1
       SForm     <- "Prospective Space-Time"
   }
   if (AnalysisType == 5) {
       xmsg      <- paste0("AnalysisType of Spatial Variation in Temporal Trends is not supported.\n")
       warning(xmsg,call.=FALSE)
       FNDError  <- TRUE
       NumErrors <- NumErrors + 1
   }
   if (AnalysisType == 6) {
       xmsg      <- paste0("AnalysisType of Prospective Purely Temporal is not supported.\n")
       warning(xmsg,call.=FALSE)
       FNDError  <- TRUE
       NumErrors <- NumErrors + 1
   }
   if (AnalysisType > 6 || AnalysisType < 1) {
      # back value
       xmsg      <- paste0("AnalysisType found in the Parameter file is unknown and not supported.\n")
       warning(xmsg,call.=FALSE)
       FNDError  <- TRUE
       NumErrors <- NumErrors + 1
   }
   
   #cat("SSForm:",SSForm," >> ",SForm,"\n")
   
   callVarList$SSForm <- SSForm
   callVarList$AnalysisType <- AnalysisType
   
   #
   #   Model Type Check
   #
      
   #
   #  All critical arguments in .prm have been evaluted.
   #
   if (ModelType == 3) {
      xmsg   <- paste0("Model Type set to Ordinal.  Ordinal Models are not supported at this time.")
      warning(xmsg,call.=FALSE)
      FNDError <- TRUE
      NumErrors <- NumErrors + 1
   }
   if (ModelType == 5) {
      xmsg   <- paste0("Model Type set to Normal.  Normal Models are not supported at this time.")
      warning(xmsg,call.=FALSE)
      FNDError <- TRUE
      NumErrors <- NumErrors + 1
   }
   if (ModelType == 7) {
      xmsg   <- paste0("Model Type set to Multinormial.  Multinormial Models are not supported at this time.")
      warning(xmsg,call.=FALSE)
      FNDError <- TRUE
      NumErrors <- NumErrors + 1
   }
   
   callVarList$ModelType <- ModelType
   
   
   
   if (FNDError) {
      xmsg    <- paste0("Error(s) found in processing Saved Session (.prm) file.  Run terminated.")
      stop(xmsg,call.=FALSE)
   }
   
   #
   #  End of Parameter File Processing.
   #
   ###############################
   
   #cat("SaTScan session save file read and analyzed.\n")
   #printCallVarList(callVarList)
   
   ####
   ###
   ##
   #
   #  End of argument and parameter checks - time to get rolling.
   #
   ##
   ###
   ####
   
   FNDError <- FALSE

   #
   ####
   
   ######
   #
   #  Initialize Global Variables - 
   #
   Save_Width <- getOption("width")
      options(width=200)
   
   #
   #  FUTURE with palColors implemented - categ must be checked against the number
   #    of colors available.
   #
   #  FUTURE allow caller to specify the break point value list. 
   #
   
   #####
   #
   #  Break point lists for categorizations.
   #
   #BreakListOE <- c(-1, 0.33, 0.4, 0.5, 0.67, 0.9, 1.1, 1.5, 2.0, 2.5, 3.0, 999999)         #  US O/E breakpoint list for Obs_Exp (0 to > 2.0)
   #
   #  New O/E break list to handle 11 levels. (12 variables) 

   
   BreakListClus <- list(c(0),c(0))
   ColorListClus <- list(c(F),c(F))
   
   BreakListClus[[3]]  <- c(-Inf,                  0.9,1.1,                    Inf)
   ColorListClus[[3]]  <- c(F,F,      T,F,F,        T,     F,F,T,F,F)
   BreakListClus[[4]]  <- c(-Inf,                  0.9,1.1,               3,0, Inf)  
   ColorListClus[[4]]  <- c(F,F,      T,F,F,        T,       F,T,F,        T,F)
   BreakListClus[[5]]  <- c(-Inf,       0.5,       0.9,1.1,               3.0, Inf)
   ColorListClus[[5]]  <- c(F,   T,F,       T,F,    T,       F,T,F,        T,F)
   BreakListClus[[6]]  <- c(-Inf,       0.5,       0.9,1.1,     2.0,      3.0, Inf)
   ColorListClus[[6]]  <- c(F,   T,F,       T,F,    T,   F,T,       T,     T,F)
   BreakListClus[[7]]  <- c(-Inf, 0.33, 0.5,       0.9,1.1,     2.0,      3.0, Inf)
   ColorListClus[[7]]  <- c(F, T,   T,      T,F,    T,   F,T,       T,     T,F)
   BreakListClus[[8]]  <- c(-Inf, 0.33, 0.5,       0.9,1.1,1.5, 2.0,      3.0, Inf)
   ColorListClus[[8]]  <- c(F, T,   T,      T,F,    T,  T,   T,     T,     T,F)
   BreakListClus[[9]]  <- c(-Inf, 0.33, 0.5, 0.67, 0.9,1.1,1.5, 2.0,      3.0, Inf)
   ColorListClus[[9]]  <- c(F, T,   T,    T,    T,  T,  T,   T,     T,     T,F)
   BreakListClus[[10]] <- c(-Inf, 0.33, 0.5, 0.67, 0.9,1.1,1.5, 2.0, 2.5, 3.0, Inf)
   ColorListClus[[10]] <- c(F, T,   T,    T,    T,  T,  T,   T,   T,   T,   T)
   BreakListClus[[11]] <- c(-Inf,0.33,0.4,0.5,0.67,0.9,1.1,1.5, 2.0, 2.5, 3.0, Inf)
   ColorListClus[[11]] <- c(  T,   T,  T,  T,   T,  T,  T,   T,   T,   T,   T)
   #  Always get 11 colors..
   
   BreakListOE <- BreakListClus[[categ]]
   #
   #   Improvement - allow user to provide the Cluster and rate break point lists 
   #
   ######

   ####
   #
   #  Color schemes for Cluster and Rate maps.
   #
   #  Future - palColors parameter.
   
   #cat("categ:",categ,"\n")
   
   #  Used with O/E Ratio and Cluster Maps  ("n" levels plus outlines)
   
   ColorsB_Clus_Full <- rev(brewer.pal(11,"RdYlBu"))
   #cat("ColorsB_Clus_Full:",ColorsB_Clus_Full,"\n")
   
   # original way of getting the colors for the legend and areas.
   #ColorsB_Clus_Mid  <- ColorsB_Clus_Full[ColorListClus[[categ]]]
   
   ColorsB_Clus_Mid   <- rev(brewer.pal(categ,"RdYlBu")) 
   #cat("ColorsB_Clus_Mid :",ColorsB_Clus_Mid,"\n")
   
   ColorsB_Clus_Low  <- ColorsB_Clus_Full[1]
   ColorsB_Clus_High <- ColorsB_Clus_Full[11]
   
   #
   #  Posible change to use high and low values in maps.
   #
   
   #
   ####
   
   ####
   #
   #  Load information structures for State, State/HSA, State/County, State/County/Census Tract
   #   and complete build of any of the index tables.
   #
   #   Get state, HSA, county and census tract location ID information 
   #    for census years 2000 and 2010.  As of v1.0.0 we don't deal
   #    boundary files any longer.
   #
   st99_data <- NULL
   hs99_data <- NULL
   co99_data <- NULL
   tr99_data <- NULL
   pl99_data <- NULL
   
   #  Step 1 - load index structures.
   
   #  Initialize - make sure old structures are GONE.
   #
   if (debugFlag) {
      ldir <- "c:/projects/statnet/r code/satscanMapper-1.0.0/data/"
      st99FN <- paste0(ldir,"st99_data.rda")
      #hs99FN <- paste0(ldir,"hs99_data.rda")
      co99FN <- paste0(ldir,"co99_data.rda")
      pl99FN <- paste0(ldir,"pl99_data.rda")   # load only if tr level data 
      tr99FN <- paste0(ldir,"tr99_data.rda")   # load only if tr level data  
      #
      load(file=st99FN)
      #load(file=hs99FN)
      load(file=co99FN)
      load(file=pl99FN)
      load(file=tr99FN)
   } else {
      data(st99_data,envir = environment(),package="satscanMapper")
      #data(hs99_data,envir - environment(),package="satscanMapper")
      data(co99_data,envir = environment(),package="satscanMapper")
      data(pl99_data,envir = environment(),package="satscanMapper")
      data(tr99_data,envir = environment(),package="satscanMapper")
   }
   
   #cat("data sets loaded for satscanMapper Z-2117 .\n")
  
   #
   #  The above loads the three working table indexes:
   #        st99_data   (State Fips to literals)
   #        hs99_data   (HSA numbers to literals)
   #        co99_data   (State/County Fips to literals)
   #        pl99_data   (State/County/Place to literals, used with tract data)
   #        tr99_data   (State/County/Census Tract Fips to literals ???)
   #
   #    Both the county and census tract index will be edited to only the required states.
   #

   ########
   #
   #  The cartesian coordinates or Lat/Long in the results files cannot be 
   #    used. We have no way of knowing the units, projection, or origin
   #    of the coordinates.  Instead, the package will key on the 
   #    Location IDs in the results and match them to our own
   #    boundaries and centroid coordinates.  It also allows us to 
   #    change the projection on request.  Our default is Albers Equal Area
   #    centered on the continental US.
   #
   ########
   
   # 
   #    census year filter = 2000 -> 1    2010 -> 2
   #
   yF   <- 1
   
   if (censusYear == "2010")  yF <- 2
   
   #cat("Set up st99.index structure for run Z- .\n")
   
   st99.index       <- st99_data
   st99.index$ID    <- row.names(st99.index)
   st99.index$stID  <- as.character(st99.index$ID)
   #stAllList       <- st99.index$stID
      
   #cat("st99.index updated - Z-2167 \n")
   #print(st99.index)
   
   
   #
   # st99.index fields
   #    row.names -> 2 digit state FIPS
   #    $abbr
   #    $stName
   #    $rgID
   #    $rgName
   #    $dvID
   #    $dvName
   #    $loc
   #    $hsa_00
   #    $hsa_10
   #    $county_00
   #    $county_10
   #    $tract_00
   #    $tract_10
   #    $change10
   #    $c_X
   #    $c_Y
   #  added:    
   #    $ID
   #    $stID
   #    $scale
   #
   #    $use
   #
   
   
   #cat("Set up hs99.index structure for run.\n")
    
   #hs99.index         <- hs99_data
   #hs99.index$ID      <- row.names(hs99.index)
   #hs99.index$HSAID   <- hs99.index$ID
   #hs99.index$stID    <- as.character(hs99.index$stID)
   #hs99.index$stName  <- as.character(hs99.index$stName)
   #hs99.index$saID    <- as.character(hs99.index$saID)
   #hs99.index$hsaName <- as.character(hs99.index$hsaName)
   
   #hsAllList        <- hs99.index$ID

   #
   # hs99.index fields
   #    row.names  -> 3 digit HSA number
   #    $hsaName
   #    $stID
   #    $stName
   #    $county_00
   #    $county_10
   #    $tract_00
   #    $tract_10
   #    $change10
   #    $c_X_00
   #    $c_Y_00
   #    $c_X_10
   #    $c_Y_10
   #  added:  
   #    $ID
   #
   #

   #cat("Set up co99.index structure for run.\n")
    
   co99.index       <- co99_data
   co99.index$ID    <- row.names(co99.index)
   co99.index$stcoID<- co99.index$ID
   co99.index$stID  <- as.character(co99.index$stID)
   co99.index$stName<- as.character(co99.index$stName)
   co99.index$saID  <- as.character(co99.index$saID)
   
   co99No           <- str_sub(co99.index$ID,3,5) == "000"
   co99.index       <- co99.index[!co99No,]
   
   #coAllList        <- co99.index$ID
   
   #
   # co99.index fields
   #    $stID
   #    $stName
   #    $coName
   #    $saID
   #    $c_X_00
   #    $c_Y_00
   #    $c_X_10
   #    $c_Y_10
   #    $tracts_00
   #    $tracts_10
   #    $y
   #
   #    $ID
   #
   #    $ID
   #    $stID
   #    $stName
   #    $saID
   #    $hsaID
   #    $hsaName
   #    $stcoID
   #
   
   #
   #  edit tables to only have data for the selected census year.
   #   Setup co99.index based on StateListDAll from Clus, Gis, RR tables
   #
   
   #cat("Set up pl99.index structure for run.\n")
   
   pl99.index       <- pl99_data
   pl99.index$ID    <- row.names(pl99.index)
   pl99.index$stID  <- str_sub(pl99.index$stcoID,1,2)
   #plAllList        <- pl99.index$ID
   
   #
   #  pl99.index fields
   #     $stcoID
   #     $plName
   #     $tracts_00
   #     $tracts_10
   #
   #     $ID
   #     $stID 
   #
   
   #cat("Set up tr99.index structure for run.\n")
   
   tr99.index       <- tr99_data
   tr99.index$ID    <- row.names(tr99.index)
   tr99.index$stID  <- str_sub(tr99.index$ID,1,2)
   tr99.index$stcoID<- str_sub(tr99.index$ID,1,5)
   #trAllList        <- tr99.index$ID
   
   tr99.index$stAbbr<- st99.index[tr99.index$stID,"abbr"]
   tr99.index$stName<- st99.index[tr99.index$stID,"stName"]
   tr99.index$coName<- co99.index[tr99.index$stcoID,"coName"]
   
   tr99.index$plName <- str_sub(tr99.index$plKey,6)
   
   # 
   #  tr99.index fields
   #    $stID
   #    $plKey
   #    $c_X_00
   #    $c_Y_00
   #    $c_X_10
   #    $c_Y_10
   #    $y
   #    
   #    $ID
   #    $hsaID
   #    $stcoID
   #    $stAbbr
   #    $stName
   #    $hsaName
   #    $coName
   #    $plName
   #
   
   #str(st99.index)
   #str(co99.index)
   #str(pl99.index)
   #str(tr99.index)
   
   suppressWarnings(rm(st99_data))
   suppressWarnings(rm(co99_data))
   suppressWarnings(rm(pl99_data))
   suppressWarnings(rm(tr99_data))
   
   #
   #  In all tables (co99 and tr99) only the states with data are included and only
   #  the entries related to the specified census year.   xxxAllList variables 
   #  are taken before the census year and state edit. (all IDs in data tables.)
   #  The DI lists after the state and year edit are in xxxSelList.
   #
   
   #
   #  2000 and 2010 Census area information is loaded.
   #
   ####

   ########
   #
   #   FILES  Locations:
   #
   #   Set up location (directory) of SatScan ResultsFiles (output).
   #
   #   The output base filename for all output SatScan files should be the same.
   #
   #  Base Directory for location of all ResultsFile output files:  
   #           (CHANGE to your environment)
   #
   
   #print("Call parameters:\n")
   #print(paste("resDir           :",resDir,sep=""))
   #print(paste0("outDir           :",outDir))
   #
   #  Set SatScan's output base filename to be used by the program - where:
   #
   #print(paste("prmFile          :",prmFile,sep=""))
   #
   #print(paste("Prm (full)       :",wPrmFileName,sep=""))
   #print(paste("Cluster Info     :",wColFileName,sep=""))
   #print(paste("Location Info    :",wGisFileName,sep=""))
   #print(paste("RR Estimated Info:",wRRFileName,sep=""))
   #
   callVarList$ColFileName <- wColFileName              # duplicate effort.
   callVarList$GisFileName <- wGisFileName
   callVarList$RRFileName  <- wRRFileName
   #
   #
   ######
   

   ######
   #
   #  The output SatScan files are assumed to be located in the same directory
   #  as the SatScan input files and have the same base filename.
   #
   #    Base directory and filename for any output reports or maps.
   #
   OutputFNBase <- ""
   OutputFNpdf  <- ""
   OutputFNtxt  <- ""
   
   if (DoOutput) {
      OutputBase <- paste0(outDir,baseDataFN)
   
      if (runId != "") {
         # use runId
         OutputFNBase <- paste0(OutputBase,"-",runId)
         OutputFNpdf  <- paste0(OutputFNBase,".pdf")
         OutputFNtxt  <- paste0(OutputFNBase,".txt")
         if (file.exists(OutputFNpdf) || file.exists(OutputFNtxt)) {
            # one of the other exists
            xmsg      <- paste0("The pdf or txt output file using the runId of ",runId," already exists.  Change runId and retry.")
            stop(xmsg, call.=FALSE)
         }
      } else {
   
         OutNum  <- 1
   
         repeat {
            OutNumLit    <- formatC(OutNum,format="f",width=2,digits=0,flag="0")
            OutputFNBase <- paste0(OutputBase,"-R",OutNumLit)
            OutputFNpdf  <- paste0(OutputFNBase,".pdf")
            OutputFNtxt  <- paste0(OutputFNBase,".txt")
            if (!file.exists(OutputFNpdf) && !file.exists(OutputFNtxt)) {
               # either file exists - good to go.
               break
            }
            OutNum <- OutNum +  1   # look for another filename
            if (OutNum > 99) {
               xmsg     <- paste0("There already exists output files with -Rxx expanded names. Cannot generate graphs or reports.")
               stop(xmsg, call.=FALSE)
            }
         }
      }
      #cat("OutputFNpdf  :",OutputFNpdf,"\n")
      #cat("OutputFNtxt  :",OutputFNtxt,"\n")
   }
   
   callVarList$OutputFNpdf <- OutputFNpdf
   callVarList$OutputFNtxt <- OutputFNtxt
     
   #
   
   #####
   #
   #     vTitle contains the main title for the page.
   #
   #     xxxxTitle2 contains the type of map and any information.
   #      this title is the second title on the page and is setup by the program based on the type of map.
   #
   #####
   
   ################## Start running new data HERE to  END
   
   ########
   #
   # Next step is to read in the Cluster Info, Location Info, RR estimates and optional
   #  rate data and validate it.
   #
   ########
   
   ######
   ##
   #
   #   START to read the data and verify fields and values.
   #
   #  Read Table outputted by SatScan - "*.col.dbf", "*.gis.dbf", and "*.rr.dbf"
   #
   #  Read Rate Table if present
   #
   ##
   ######
   
   #####
   #
   #  Cluster Information File   *.col.dbf
   #
   #  The possible combinations for the *.col.dbf file fields are:
   #
   #  Fields:         Of interest
   #   CLUSTER           *-all    
   #   LOC_ID            *-all
   #                               Lat/Long or X,Y,Zn
   #   LATITUDE        Y
   #   LONGITUDE       Y
   #   X                    Y
   #   Y                    Y
   #   Z1                   n
   #   Z2                   n
   #   Z3                   n
   #                              Radius or Ellipical
   #   RADIUS          Y-->
   #   E_MINOR              Y
   #   E_MAJOR              Y
   #   E_ANGLE              Y
   #   E_SHAPE              Y
   #
   #   START_DATE        *-all?
   #   END_DATE          *-all?
   #   NUMBER_LOC        *-all
   #   LLR               *-all?
   #   TEST_STAT         n
   #   P_VALUE           *-all
   #   RECURR_INT        n           if temporal.
   #   OBSERVED              *-most
   #   EXPECTED              *-most
   #   ODE                   *-most
   #   REL_RISK              *-most
   #
   #   WEIGHT_IN         n
   #   MEAN_IN           n
   #   MEAN_OUT          n
   #   VARIANCE          n
   #   STD               n
   #   W_MEAN_IN         n
   #   W_MEAN_OUT        n
   #   W_VARIANCE        n
   #   W_STD             n
   #
   #  If Observed, Expected and Obs/Exp not present, analysis not supported.
   #
   #  Not all fields are present in every run.   The fields that must be present for this
   #  package are:
   #      CLUSTER, LOCATION_ID, 
   #      LATITUDE/LONGITUDE or X/Y (???), 
   #      RADIUS or E_MINOR, E_MAJOR, E_ANGLE, E_SHAPE,  (could outline using boundaries of areas.)
   #      P_VALUE, OBSERVED, EXPECTED, ODE.  (must for our version!)
   #
   #  Optional fields:
   #      STARTDATE, and ENDDATE  (required for Space-Time, but not for Spatial Only),
   #      NUMBER_LOC
   #
   ######
   
   ###
   #
   #  ### Need to determine if SatScan run was Spatial, Temporal, or Spatial/Temporal.
   #
   #         Plotting Spatial is all or nothing.
   #         Plotting Temporal - no plot - it's all time - not handled by this program.
   #         Plotting Spatial/Temp - plots are over period of time and space.
   #            Program setup to handle period in years.
   #
   #   Have this information from the PRM file.
   #
   #  Th Cluster Information fields and formats don't match the SatScan documentation - 
   #  an extra column appears to be added after "P_Value" -> RECURR_INT.  Best way 
   #  to handle these variations in the Cluster Information file is to use the DBF format
   #  that carries the field names. 
   #
   #  Load Cluster (Col) Table
   #
   ColTable    <- NULL     # clear the table 
   
   #cat("SaTScan Cluster Information file",wColFileName,"\n")
   
   ColTable      <- foreign::read.dbf(wColFileName, as.is=TRUE)
   
   #cat("Cluster Info - field names:",names(ColTable),"\n")
   
   ColTable      <- as.data.frame(ColTable)
   
   ColTableNames <- colnames(ColTable)
   ColTableLen   <- dim(ColTable)[1]
  
   #
   #    Cluster Locations
   #
   res              <- CleanLocID(ColTable)
   ColTable         <- res$Table
   ColLocList       <- res$LocList
   ColStateList     <- res$StateList
   ColIDType        <- res$IDType
  
   ColTableNames    <- colnames(ColTable)   # get column names - UPDATE
   
   #
   #  Cluster Info. - OBSERVED, EXPECTED and Obs/Exp fields.
   #
   ColHasODE   <- TRUE
   
   if (is.null(ColTable$OBSERVED) 
         || is.null(ColTable$EXPECTED) 
         || is.null(ColTable$ODE) ) {
     
       # The observed, expected, and O/E fields are missing from the Cluster Information Results File.
       
       xmsg = paste0("Cluster Information file does not contain 'Observed', 'Expected', or 'Observed/Expected' data.\n",
               "The SaTScan results cannot be mapped.\n",
               "Change SaTScan run parameters to an analysis that will generate these variables and re-run package.\n")
       warning(xmsg,call.=FALSE)
       FNDError  <- TRUE
       NumErrors <- NumErrors + 1
       ColHasODE <- FALSE
   }
   
   #
   # verify ColTable required fields are present.
   #
   
   #
   #  Cluster Info - RADIUS and Ellipical fields
   #
   # Check ColTable for Circle and Ellipical fields and convert to Ellipical.
   
   if (is.null(ColTable$RADIUS) && is.null(ColTable$E_MINOR)) {  # if Radius and Ellipical are missing
      # if Radius and Ellipical are missing...
      
      xmsg      <- paste0("Cluster Information file does not contain either cluster outline radius or ellipical data.\n",
                    "Outlines of clusters cannot be drawn.  Verify the correct files are being referenced.\n")
      warning(xmsg,call.=FALSE)
      FNDWarn   <- TRUE 
      NumWarns  <- NumWarns + 1 
      outline   <- FALSE     # no data for outlines - turn it off
   }
   
   #
   #  Cluster Info - Shape & process Radius and Ellipical paramaters.
   #
   
   ClusArea     <- "Elliptical"
   ClusSSShape  <- 1
   ClusAreaType <- TRUE

   if (!is.null(ColTable$RADIUS)) {
      # RADIUS field is in the data.  Convert to Ellipical information and added to data.
      # Use Radius for Ellipse parameters (major and minor)
          
      ColTable$E_MAJOR <- ColTable$RADIUS    # in KM
      ColTable$E_MINOR <- ColTable$RADIUS 
      ColTable$E_ANGLE <- 0
      ColTable$E_SHAPE <- 0
      ClusArea         <- "Circular"
      ClusAreaType     <- FALSE
      ClusSSShape      <- 0

   }
   
   #
   #  if Ellipical - maj and min are in "units" of projection.
   #
   
   ColTable$RADIUS <- 0   # reset Radius field -> it's now an Ellipse
               # or add RADIUS for consistancy if elliptical data.
   
   #
   #  Cluster SHAPE  vs.  information available in file.
   #
  
   if (ClusSSShape != SSShape) {
      # The run parameter for shape should match the data we see in the results cluster file.
      # It doesn't.
      xmsg      <- paste0("The saved session (.prm) file specifies ",SType,
                     " type clusters.  The Cluster Information file defines ",ClusArea,
                     " type clusters.\n",
                     "They do not match - check to ensure the correct files are being referenced.\n")
      warning(xmsg, call.=FALSE)
      FNDError  <- TRUE
      NumErrors <- NumErrors + 1
   }
  
   #
   #  Cluster Info - START_DATE
   #
   ColHasDates   <- TRUE

   if (is.null(ColTable$START_DATE)) {
      # dates don't exist - we have Spatial ONLY data add fake dates.
      
      ColHasDates          <- FALSE
      # set up fake.
      ColTable$START_DATE  <- "2000/01/01"
      ColTable$END_DATE    <- "2000/12/31"
      ColTable$StrDate     <- as.Date(ColTable$START_DATE,"%Y/%m/%d")
      ColTable$EndDate     <- as.Date(ColTable$END_DATE,"%Y/%m/%d")
      
   } else {
      # dates do exist - fix them up   - convert from character to "Date" class.
      ColHasDates          <- TRUE    
      ColTable$StrDate     <- as.Date(ColTable$START_DATE,"%Y/%m/%d")
      ColTable$EndDate     <- as.Date(ColTable$END_DATE,"%Y/%m/%d")
   }
   minDate <- min(ColTable$StrDate)
   maxDate <- max(ColTable$EndDate)
   
   #   ColHasDates indicates real or fake.   From now on, Cluster always has a date - real or fake.
   
   if (!ColHasDates && SSForm == 1) {
      # no dates in cluster records, but S-T type analysis.
      
      # This may not be a problem - have to handle as Spatial Only - so dates would only be for S-T.

      # attempting to do Space-Time, but cluster information does not have any dates (Spatial Only).
      xmsg      <- paste0("The parameter file specifies ",SForm,
                      " type analysis, but the Cluster Information file does not have any date information.\n",
                      "Check to ensure the correct files are being referenced.\n")
      warning(xmsg,call.=FALSE)
      FNDError  <- TRUE
      NumErrors <- NumErrors + 1
   }
    
   #
   #  Cluster Info. - pValue
   #
   ColHasPValue   <- TRUE
   
   if (is.null(ColTable$P_VALUE)) {
      # pValue field is missing.
      ColHasPValue  <- FALSE
       
   } else {    
      # have pValue column in cluster info.
      
      wPV_NA <- all(is.na(ColTable$P_VALUE))     # TRUE if all NAs.
      
      if (wPV_NA) {
         # all of the pValues in the cluster table are NA.
         xmsg     <- paste0("All P_VALUES in the cluster results file are set to NA.\n", 
                     "Check SaTScan run parameters.\n")
         warning(xmsg,call.=FALSE)
         FNDError     <- TRUE
         NumErrors    <- NumErrors + 1
         ColHasPValue <- FALSE

      } else {  
         # Have some pValues in the data.
         
         wPV <- any( ColTable$P_VALUE < pValue)   # TRUE - Cluster to map.
   
         if (!wPV)  {
            # None of the clusters have pValues < 0.05 (or the defined limit), so no 
            # clsuters will be mapped.
            xmsg      <- paste0("There are no clusters in the data that have a P_VALUE less than ",pValue,"\n",
                           " There is nothing to map.\n")
            warning(xmsg,call.=FALSE)
            FNDError  <- TRUE
            NumErrors <- NumErrors + 1
         
         }               
      }               
   }               
                  
   #
   #  Check if error flagged
   #
   if (FNDError) {
      xmsg   <- paste0("Error(s) found in reading and validating the Cluster information results file.\n",
                    "Run Terminated.\n")
      stop(xmsg,call.=FALSE)
   }   
   
   FNDError <- FALSE
   
   # we will not get this far if there are no dates (or fakes), Obs/Exp, etc.

   ####   
   #
   #  Cluster Info - Process Date and get the year for each date.
   #     real or fake.
   
   ColYrList     <- NULL
   ColTable$NYr  <- NULL
   
   # Cluster info has datas, clean them up and convert to YEARS.
   ColTable$SYr  <- getYear(ColTable$StrDate)
   ColTable$EYr  <- getYear(ColTable$EndDate)
  
   ColTable$NYr  <- ColTable$EYr - ColTable$SYr + 1
  
   MinYr         <- min(ColTable$SYr,ColTable$EYr)
   MaxYr         <- max(ColTable$SYr,ColTable$EYr)
   ColYrList     <- c(MinYr:MaxYr)
   
   #cat("Years referenced in the Cluster Info :")
   #cat(ColYrList,"\n",sep="  ")
   #cat("Years min:", MinYr, "  max:",MaxYr,"\n")
   
   #  
   #  Col Table data has been read.
   #
   ###
   
   LClusCol         <- sort(unique(ColTable$CLUSTER))
   NumClusCol       <- length(LClusCol)
   
   #cat("   ",NumClusCol," clusters found in the Col Info. results file.\n")
   #cat("    List of cluster numbers:\n")
   #cat("       ",LClusCol,"\n")
   
   ###
   #
   #  Make adjustments to ColTable and expand 
   #
  
   ColTable$HL      <- 0
   ColTable$Cat     <- 1
   ColTable$Col     <- 0
   ColTable$Lab     <- ""
   
   ColTableSize     <- length(ColTable$CLUSTER)
   
   if (ColHasODE) {
      #cat("Calculated ClusCategories.\n")
   
      #    High, Middle, Low cluster (1, 0, -1) 
      ColTable$HL   <-  sign(ColTable$ODE - 1.0)   # 1 > 1.0 (High), -1 < 1.0 (Low),  0 = 1.0
   
      #   Categorize Cluster and Assign Color.
   
      ColTable$Cat  <- cut(ColTable$ODE,breaks=BreakListOE)
  
      ColTableCat   <- levels(ColTable$Cat)
      ColTable$Col  <- ColorsB_Clus_Mid[as.integer(ColTable$Cat)]
  
      #  Expanded cluster label and number of years cluster has existed.
      
      #cat("label:",label,"  SSForm:",SSForm,"  MaxTemporalSize:",MaxTemporalSize,"\n")
     
      if (label) {
      # if (label && SSForm == 1 ) {
         # label requested and Space-Time analysis.
         if (MaxTemporalSize > 1 && SSForm == 1) {    # why?
            # 
            ColTable$Lab <- paste0(as.character(ColTable$CLUSTER),"-",as.character(ColTable$NYr))
         } else {
            ColTable$Lab <- as.character(ColTable$CLUSTER)
         }
      }
   }
   #
   ####
   
   #cat("ColTable Z-2312 :\n")
   #print(head(ColTable))
   #print(table(ColTable$Cat))
   #print(table(ColTable$Col))
   #cat("ColTable-ODE range:",range(ColTable$ODE),"\n")
   
   ####
   #
   #  Move check for valid P_VALUE and cluster numbers > 0 to here - centralize Cluster Information checks.
   #
   ####
      
   ####
   #
   # Future Code (NO GO on this)
   #
   # 1) Provide alternative information when no Radius or Ellipical parameters are provided by SatScan.
   # 2) Ensure the cluster outline is just outside of the cluster and allows the tracts to be seen.
   #
   ####
    
   ####
   #
   #  Read the GIS File from SatScan
   #
   #  *.gis.dbf   GIS file format (columns)  
   #
   #
   #  1 Location_ID
   #  2 Cluster_number  ***
   #  3 P-Value of Cluster ***
   #  4 Observed Cases in cluster
   #  5 Expected Cases in cluster
   #  6 Obs/Expect in cluster (OBS_EXP)
   #  7 Obs Cases at location
   #  8 Exp Cases at location
   #  9 Obs/Exp at location (OBS_EXP)
   # 10 Relative Risk (Optional)
   #
   #  Same format in all cases.
   #
   #  GIS - one record per location/cluster pair - no more - no time relationship 
   #     only key is through cluster's time period.
   #     Does not include all locations.  Only those in clusters.
   #
   
   FNDError  <- FALSE
   
   #cat("SaTScan Location Info file:",wGisFileName,"\n")
   
   #
   #  Later add code to look at ext of file and use TXT or DBF as needed to read file.
   #
   
   GisTable <- NULL
   
   # dbf format    - *.dbf
   GisTable <- foreign::read.dbf(wGisFileName, as.is=TRUE)
   
   #cat("Gis/Loc Info - orig field names:",names(GisTable),"\n")
   
   GisTable      <- as.data.frame(GisTable)
   
   GisTableNames <- colnames(GisTable)
  
   oldCNames     <- c("CLU_RISK","LOC_RISK","OBS/EXP")
   newCNames     <- c("CLU_RR","LOC_RR","OBS_EXP")
   xM     <- match(GisTableNames,oldCNames)
   if (any(!is.na(xM))) {
      # found at least one match
      GisTableNames[!is.na(xM)] <- newCNames[xM][!is.na(xM)]
      colnames(GisTable) <- GisTableNames
   }
      
   #cat("Gis/Loc Info - new field names:",names(GisTable),"\n")
   
   #
   # Fields possible:
   #    LOC_ID
   #    CLUSTER
   #    P_VALUE
   #
   #    RECURR_INT
   #    CLU_OBS
   #    CLU_EXP
   #    CLU_ODE
   #    CLU_RISK  -> CLU_RR (v9.4.x)
   #    LOC_OBS
   #    LOC_EXP
   #    LOC_ODE
   #    LOC_RISK  -> LOC_RR (v9.4.x)  
   #
   #    LOC_X   or LOC_LATITUDE
   #    LOC_Y   or LOC_LONGITUDE
   #
   #    GINI_CLUST
   #
   
   # Old --- GisTableNCol <- dim(GisTable)[2]  # get number of columns (variables)
   
   res              <- CleanLocID(GisTable)
   GisTable         <- res$Table
   GisLocList       <- res$LocList
   GisStateList     <- res$StateList
   GisIDType        <- res$IDType
   
   GisHasODE  <- TRUE
  
   if (is.null(GisTable$CLU_OBS) || is.null(GisTable$CLU_EXP) 
          || is.null(GisTable$CLU_ODE) || is.null(GisTable$LOC_OBS) 
          || is.null(GisTable$LOC_EXP) || is.null(GisTable$LOC_ODE))  {

      # If any of the following are missing:  Cluster observation, Cluster expected,
      #     Cluster ODE, location observation, location expected, location ODE
      #     are missing - report error to user and signal found error.
      #
   
      xmsg     <- paste0("The Location Information file does not contain one or more of the following:\n",
                         "  Observed, Expected and Observed/Expected data.\n",
                         "  Mapping requires these variables.\n")
      warning(xmsg,call.=FALSE)
      FNDError  <- TRUE
      NumErrors <- NumErrors + 1
      GisHasODE <- FALSE
   
   }
     
   #  location information generally has no dates.
 
   GisHasDates   <- TRUE
   
   if (is.null(GisTable$START_DATE)) {
      # start date field is not present in location results file.
      # cluster info has dates
         
      GisHasDates <- FALSE 
      
      GisTable$StrDate <- ColTable[GisTable$CLUSTER,"StrDate"]
      GisTable$EndDate <- ColTable[GisTable$CLUSTER,"EndDate"]
      GisTable$SYr     <- ColTable[GisTable$CLUSTER,"SYr"]
      GisTable$EYr     <- ColTable[GisTable$CLUSTER,"EYr"]
  
   } else {
      # we have dates.
      
      GisHasDates          <- TRUE    
      GisTable$StrDate     <- as.Date(GisTable$START_DATE,"%Y/%m/%d")
      GisTable$EndDate     <- as.Date(GisTable$END_DATE,"%Y/%m/%d")
      GisTable$SYr         <- getYear(GisTable$StrDate)
      GisTable$EYr         <- getYear(GisTable$EndDate)
   }   
   GisTable$NYr  <- GisTable$EYr - GisTable$SYr + 1
     
   GMinYr        <- min(GisTable$SYr,GisTable$EYr)
   GMaxYr        <- max(GisTable$SYr,GisTable$EYr)
   GisYrList     <- c(GMinYr:GMaxYr)

   #cat("Years referenced in the Location Info :")
   #cat("  ",ColYrList,"\n",sep="  ")

   #####
   #
   #   An assumption is made that is GIS has dates (not supposed to), then they will
   #   be within the dates of the Cluster Info.  Dates and Years are more important at 
   #   the cluster level.
   #
   #####

   if (FNDError) {
      xmsg   <- paste0("Error(s) found in reading and validating the Location information results file.\n",
                    "Run Terminated.\n")
      stop(xmsg,call.=FALSE)
   }   

   FNDError <- FALSE
  
   #
   #  GIS file has been read 
   #
   ###

   LClusGis       <- sort(unique(GisTable$CLUSTER))
   NumClusGis     <- length(LClusGis)
   
   #cat("   ",NumClusGis," clusters found in the GIS input file.\n")
   #cat("     List of cluster number in Location Info:\n")
   #cat("    ",LClusGis,"\n")
   
   ###  
   #
   #  Col and Gis Consistancy Test
   #
   #  None at the current time.
   #
   ###
     
   ###
   #  
   #  Extra values added to data.frame
   #
   #   Location and Cluster O/E  - Category/Colors
   #
   
   GisTable$LocCol   <- 0
   GisTable$LocCat   <- 1
   GisTable$ClusCol  <- 0
   GisTable$ClusCat  <- 1
   GisTableCat       <- ""
   
   if (GisHasODE) {
      #cat("Calculated Local Categories.\n")
      
      # We have ODE info, so can do the categorizing for location information
      
      #  Set category based on loc ODE
      GisTable$LocCat    <- cut(GisTable$LOC_ODE,breaks=BreakListOE)
   
      # build category list
      GisTableCat        <- levels(GisTable$LocCat)
 
      # convert table into characters for legend
      GisTableCat[1]     <-  str_replace(GisTableCat[1],"\\(-1,","\\[0,")
      GisTableCat[categ] <-  str_replace(GisTableCat[categ],"1e\\+07\\]","Inf\\]")
 
      # pick up appropriate color
      GisTable$LocCol    <- ColorsB_Clus_Mid[as.integer(GisTable$LocCat)]   # color based on location O/E Cat
   
      #   Cluster O/E - Category   - copy clusters color      
      GisTable$ClusCat   <- ColTable[GisTable$CLUSTER,"Cat"]
      GisTable$ClusCol   <- ColTable[GisTable$CLUSTER,"Col"]
   }    
   
   #
   ####
   
   #cat("Gis/Loc Table Z-2524 :\n")
   #print(head(GisTable))
   #print(table(GisTable$LocCat))
   #print(table(GisTable$LocCol))
   #print(table(GisTable$ClusCat))
   #print(table(GisTable$ClusCol))
   #cat("GisTable-ODE range:",range(GisTable$LOC_ODE),"\n")
   
   
 ##  ####
   #
   #  Check the number of clusters in the Cluster and Location data.
   #
   mapClusList  <- unique(LClusCol,LClusGis)
   mapClusNum   <- length(mapClusList)
   
   if (NumClusCol != NumClusGis) {
      # the number of clusters in the Cluster and Location results differ.
      # Situations:
      #    1 - More clusters in Cluster than Location.
      #        Don't have detailed locations to map cluster, only outline and centroid.
      #        Can't do locO_EMap.   Need RR is present to provide list of areas.
      #        
      #        Use RR to map all areas.  
      #        Create polygon of circle or ellipse.
      #        Find areas in each cluster space for clusters without location records.
      #        Build fake GisTable entries with no locO_E info only cluster info.
      #
      #        or
      #
      #        Warning to user and map only clusters with cluster and location records.
      #
      #    2 - More clusters in Location than Cluster.
      #        Have locations to map, but no cluster ODE data for Cluster Map.
      #        Don't have circle/ellipse data on cluster.
      #
      #        Create Cluster records for missing clusters using Location Info.
      #        For Cluster ODE, sum location Obs and Exp and calculate ODE.
      #        Disable outline for cluster.
      #        
      #        or 
      #
      #        Warning to user and map only clusters cluster and location records. 
      #
      #    3 - No location (gis) records for an cluster with P_VALUE = 1.
      #        Since we don't know how true and consistant this is.  Let the location records
      #        lead the way.  The cluster records are in P_VALUE decending order.
      #        So, if the location records match the "first n" records, it's OK.
      #
      #
      
      GisCList <- sort(unique(GisTable$CLUSTER))
      GisCLen  <- length(GisCList)
      if (GisCLen != GisCList[GisCLen]) {
         # The location records do not match up with the first "n" cluster records.
         # If the number of clusters in each are not the same, the must be at least the first "n"
         
         xmsg     <- paste0("The number of clusters in the Cluster and Location results files are not the same.\n",
                               "Only the clusters with both Cluster and Location records will be mapped.\n")
         warning(xmsg, call.=FALSE)
      }
      mapClusList <- intersect(LClusCol,LClusGis)
      mapClusNum  <- length(mapClusList)
      
   }
   #
   #
   ####
   
   ##
   #
   #  RR dbf File
   #
   #  Fields:
   #     Location ID         LOC_ID
   #     Observed Cases      OBSERVED
   #     Expected Cases      EXPECTED
   #     Observed/Expected   ODE
   #     Relative Risk       REL_RISK
   #
   #  The only reason for reading this file is to get a complete list of the areas covered by 
   #  the SaTScan (TM) analysis.
   #
   #
    
   RRLocTable   <- NULL
   RRIDType     <- GisIDType
   RRStateTable <- NULL
   RRTable      <- NULL
   
   if (NoRRFile) {
      xmsg     <- paste0("No RR File (dbf format) was found. Will attempt to run without it.\n")
      warning(xmsg,call.=FALSE)
      FNDWarn  <- TRUE
      NumWarns <- NumWarns + 1
    
   } else {
   
      #cat("SaTScan Relative Risk Estimates file:",wRRFileName,"\n")
  
      #
      #  Later add code to look at ext of file and use TXT or DBF as needed to read file.
      #
      # dbf format    - *.dbf
        
      RRTable <- foreign::read.dbf(wRRFileName, as.is=TRUE)
      
      #cat("RRTable read - field names:",names(RRTable),"\n")

      if (is.null(RRTable$LOC_ID))  {
         xmsg        <- paste0("RR Estimated Information does not contain LOCATION IDs.\n",
                        "File cannot be used to determine area covered by analysis.\n")
         warning(xmsg,call.=FALSE)
         FNDWarn     <- TRUE
         NumWarns    <- NumWarns + 1
      } else {
         res         <- CleanLocID(RRTable)
         RRTable     <- res$Table
         RRLocList   <- res$LocList
         RRStateList <- res$StateList
         RRIDType    <- res$IDType    
      }
   }    
   RRTable <- NULL    # release unneeded space
   rm(RRTable)
   
   #
   # based on RR file, get list of ALL locations  - That is the only reason
   #  to read the RR file - it has an entry for each location ID.
   #
   ####
   
   #####
   ####
   ###
   ##
   #
   #  At this point, all data has been read and processing into data.frames.
   #
   ##
   ###
   ####
   #####
   
   #####
   #
   #  Check to make sure all result files are using the same type of US.Fips codes.
   
   if ( ColIDType != GisIDType || ColIDType != RRIDType) {
      # The ID Types are not the same.
      xmsg      <- paste0("The LOC_IDs in the Cluster, Location and RR results files are not the same type or length.  Please correct and rerun,\n")
      stop (xmsg,call.=FALSE)
   }
   LocIDType <- GisIDType     # 2, 3, 5, or 11
   
   #
   #####
   
   ######
   #
   #  Initialize SeerMapper variable lists 
   #
   #require(SeerMapper)  # should not be used, the SeerMapper::XXX should load the package.
   
   rPM       <- NULL
   MV        <- NULL
   
   rPM       <- SeerMapper::SM_GlobInit()  # initialize Seer Mapper variable named lists.
   
   #cat("completed SM_GlobInit Z-3237 .\n")
   
   #
   #  Now we have the rPM named list to start filling in.
   #  Most of it must be done before calling SM_Build
   #
    
   rPM$proj4          <- proj4       # save projection proj4string character string for SeerMapper
   rPM$CRSproj4       <- CRSproj4
   
   rPM$palColors      <- palColors    # save palColors scheme, however the colors are passed along as the values.
   
   #   Set up the boundary options the way satscanMapper wants them.
   rPM$regionB        <- "NONE"
   regionB            <- "NONE"
   rPM$regionB_caller <- TRUE
   
   rPM$stateB         <- "NONE"
   stateB             <- "NONE"
   rPM$stateB_caller  <- TRUE
   
   rPM$seerB          <- "NONE"
   seerB              <- "NONE"
   rPM$seerB_caller   <- TRUE
   
   rPM$hsaB           <- "NONE"
   hsaB               <- "NONE"
   rPM$hsaB_caller    <- TRUE
   
   rPM$countyB        <- "NONE"
   countyB            <- "NONE"
   rPM$countyB_caller <- TRUE
   
   rPM$tractB         <- "NONE"
   tractB             <- "NONE"
   rPM$tractB_caller  <- TRUE
   
   # No hatching
   
   rPM$hatch          <- FALSE
   rPM$HatchFlag      <- FALSE   # dont do hatching
   rPM$hatch2         <- FALSE
   rPM$Hatch2Flag     <- FALSE
   
   # Legends
   
   rPM$mLegendFlag    <- FALSE   # dont do legends
   
   #cat("Basic SeerMapper callparms set - Z-3285 - LocIDType:",LocIDType," \n")
   
   mapB               <- "STATE"  #### temp until option added.
   
   #  parameter:   mapB <- c("DATA","COUNTY", "HSA", "STATE","REGION")   # default = "DATA"
   
   if (LocIDType == 2) {
      # state location IDs
      stateB  <-"DATA"       #   option for "REGION" or "ALL"
      if (mapB=="REGION") { 
          regionB <- "DATA"
          stateB  <- "REGION"
      }
      #  mapB = "STATE", "HSA", or "COUNTY" is invalid
   }
   if (LocIDType == 3) {
      # HSA location IDs (data)
      hsaB <- "DATA"
      if (mapB=="REGION") {
         regionB <- "DATA"
         stateB  <- "REGION"
         hsaB    <- "STATE"
      }
      if (mapB=="STATE") {
         regionB <- "NONE"
         stateB  <- "DATA"
         hsaB    <- "STATE"
      }
      # mapB = "HSA" or "COUNTY" is not valid.
   }
   if (LocIDType == 5) {
      # County location IDs
      countyB <-"DATA"
      if (mapB=="REGION") {
         regionB  <- "DATA"
         stateB   <- "REGION"
         hsaB     <- "STATE"
         countyB  <- "STATE"
      }
      if (mapB == "STATE") {
         regionB  <- "NONE"
         stateB   <- "DATA"
         hsaB     <- "STATE"
         countyB  <- "STATE"
      }
      if (mapB == "HSA") {
         regionB  <- "NONE"
         stateB   <- "NONE"
         hsaB     <- "DATA"
         countyB  <- "HSA"
      }
      # mapB = "COUNTY" is not valid
   }
   if (LocIDType == 11) {
      # Census Tract location IDs.
      tractB  <- "DATA"
      if (mapB=="REGION") {
         regionB  <- "DATA"
         stateB   <- "REGION"
         hsaB     <- "STATE"
         countyB  <- "STATE"
         tractB   <- "STATE"
      }
      if (mapB == "STATE") {
         regionB  <- "NONE"
         stateB   <- "DATA"
         hsaB     <- "STATE"
         countyB  <- "STATE"
         tractB   <- "STATE"
      }
      if (mapB == "HSA") {
         regionB  <- "NONE"
         stateB   <- "NONE"
         hsaB     <- "DATA"
         countyB  <- "HSA"
         tractB   <- "COUNTY"
      }
      if (mapB == "COUNTY") {
         regionB  <- "NONE"
         stateB   <- "NONE"
         hsaB     <- "NONE"
         countyB  <- "DATA"
         tractB   <- "COUNTY"
      }
   }
   
   # put the parameters back into the rPM list for SeerMapper
   
   rPM$regionB    <- regionB
   rPM$stateB     <- stateB
   rPM$seerB      <- seerB
   rPM$hsaB       <- hsaB
   rPM$countyB    <- countyB
   rPM$tractB     <- tractB
   
   #
   #####
       
   #####
   #
   #  Loc IDs must be validated against the boundary database.  The placename table for 
   #  census tracts should provide enough to do the job.
   #
   #  Also classify the mapping as - State, County (in state), or CensusTract (in state).
   #  Use *.geo to findout the area that was analyzed by SaTScan (TM).
   #
   #  since all results files are generated by SaTScan (TM), it should be impossible 
   #  for the LocIDs to have different lengths.
   #
   #  We have the ClusLocList, GisLocList, and RRLocList to compare.
   #
   #cat("Setup state and location lists - Z-3400 \n")
      
   StateListDAll  <- sort(unique(c(ColStateList,GisStateList,RRStateList)))
   StateListData  <- sort(unique(c(ColStateList,GisStateList)))
   
   LocListDAll    <- sort(unique(c(ColLocList,GisLocList,RRLocList)))
   LocListData    <- sort(unique(c(ColLocList,GisLocList)))
   
   LocNumNAAll    <- is.na(LocListDAll)                              # find NA values
   
   LocNumTestAll  <- (unlist(gregexpr("^[0-9]*$",LocListDAll)) < 0)  # set true if not a number
   LocNumTestAll[LocNumNAAll] <- TRUE     # make NA values TRUE
   
   #  test 1 = are the LOC_ID values numeric?
   
   if (any(LocNumNAAll)) {
      # some of the location IDs are NA
      xmsg   <-  paste0("Some of the LOC_ID values in the results are set to NA.  Correct and rerun.\n")
      warning(xmsg,call.=FALSE)
      FNDError <- TRUE
      NumErrors <- NumErrors + 1
   }
   if (any(LocNumTestAll)) {
      # some of the location IDs are not numbers ???
      xmsg   <-  paste0("Some of the LOC_ID values in the results are not numerical values.  Correct and rerun.\n")
      warning(xmsg,call.=FALSE)
      FNDError <- TRUE
      NumErrors <- NumErrors + 1
   }

   #
   # Build out dataMapDF  before call to SM_BUILD
   #
   
   idList          <- LocListDAll                   # all Loc_IDs from Cluster, Loc and RR files.
   lenIdList       <- length(idList)
   #cat("Number of Locations IDs found Z-3432 :",lenIdList,"\n")
   
   dataList        <- rep("white",length(idList))
   hDataList       <- rep(NA,length(idList))
    
   dataMapDF       <- data.frame(ID=idList, data=dataList, hData=hDataList, 
                         stringsAsFactors=FALSE)
   
   #  fill out dataMapDF
   cNA                <- as.character(NA)
   dataMapDF$good     <- TRUE      # all rows are valid at this time.
   # following variables get filled in based on level of Location ID.
   dataMapDF$rgID     <- cNA       # region ID
   dataMapDF$stID     <- cNA       # state  ID
   dataMapDF$saID     <- cNA       # Seer Registry ID
   dataMapDF$hsaID    <- cNA       # health district ID
   dataMapDF$stcoID   <- cNA       # state/county ID
   dataMapDF$stcotrID <- cNA       # state/county/tract ID
   dataMapDF$plName   <- cNA       # place Name
   dataMapDF$cat      <- 0         # data category #
   dataMapDF$col      <- "white"   # color
   dataMapDF$hDen     <- as.integer(NA)  # hatching density
   dataMapDF$hCol     <- "grey50"  # hatching color
   dataMapDF$hRes     <- FALSE     # hatching T/F  (test results) - does sub-area get hatched?
   dataMapDF$hRes2    <- FALSE     # hatching T/F  (test results) - does sub-area get hatched?
   
   row.names(dataMapDF) <- dataMapDF$ID

   #cat("dim(dataMapDF) Z-3460 :",dim(dataMapDF),"\n")
   
   #cat("dataMapDF - init-empty:\n")
   #print(str(dataMapDF))
   #print(head(dataMapDF))
   
   rPM$dataMapDF  <- dataMapDF
   
   #
   #  Must set xxxxB controls.
   #
   # have the rPM and dataMapDF ready for SM_Build.
   #
   ####
   
   ####
   #
   #  Named List from SM_GlobInit()
   #
   #    rPM <- SM_GlobInit()
   #
   #
   #    Setup before call:
   #          rPM$censusYear, cYear, cY = "2000", "00"
   #          rPM$stateSelDel = ""
   #          rPM$ndfName = <name of .prm file>
   #          rPM$categMode = 4
   #          rPM$dataMapDF with $ID.
   #          rPM$xxxxxB set 
   #
   #          rPM$palColors
   #          rPM$proj4
   #
   #    xRes <- SM_Build(rPM)
   #          rPM$censusYear, rPM$cYear, rPM$cY
   #          rPM$stateSelDel (?)
   #          rPM$OrigCRS, ProjCRS
   #          rPM$ndfName (>)
   #          rPM$idMode (returned)
   #          rPM$categMode (???)
   #          rPM$dataMspDF  (must have build before this routine.
   #
   #          rPM$proj4, rPM$CRSproj4    -- transformation of returned polygons.
   #
   #        Returns:
   #          rPM$idMode                 -->   type of map.
   #          MV$data_proj, data_data, dataListData
   #          MV$states_proj, states_data,StateListDAll, StateListData
   #          MV$SeerRegs_proj, SeerRegs_data, SeerRegListAll, SeerRegListData
   #          MV$regions_proj, regions_data, RegionsListAll, RegionsListData
   #          MV$hsa_proj, hsa_data, HsaListAll, HsaListData  # ******
   #          MV$co_proj,co_data,CountyListAll, CountyListData
   #          MV$tr_proj,tr_data,TractListAll, TractListData
   #          MV$co99_mapr
   #          Internal Calls:
   #             SM_SetDef(rPM) -> rPM
   #             SM_Impl_B(rPM,MV) -> MV
   #             SM_box_sel(rPM,MV) ->xRes$rPM and xRes$MV
   #          xRes$rPM and xRes$MV
   #          
   #
   #    xRes <- SM_Mapper(rPM,MV)
   #          variables used:
   #          rPM$ColorB_Data, ColorB_O_Tract... Region
   #          rPM$dataBCol,
   #          rPM$dataMapDF
   #
   #          rPM$proj4String,
   #          rPM$palColors,
   #   
   #          MV$data_proj_sel
   #          MV$data_data_sel    ($col holds the sub-area color.)
   #          MV$rg_proj_sel ... tr_proj_sel
   #          MV$xlPlot,ylPlot
   #          MV$regionPList...tractBList
   #          MV$rgGO ... trGO
   #          data_data$col (color of subarea)
   #
   #          rPM$dataBCol  (color of border)
   #          rPM$HatchFlag
   #          rPM$Hatch2Flag 
   #
   #          MV$trSP, ... MV$rgSP
   #          return(xRes$xyBox)
   #

   #####
   #
   #  Setup for Call SM_Build
   #
   rPM$categMode      <- 4    # "COLORS"
   rPM$categ          <- "COLORS"
   
   rPM$censusYear     <- censusYear
   rPM$cYear          <- str_sub(censusYear,-2,-1)
   rPM$cY	      <- ""
   if (rPM$cYear == "10")  rPM$cY <- "10"

   rPM$stateSelDel    <- ""
   rPM$ndfName        <- prmFile
   #cat("length of dataMapDF before build:",dim(rPM$dataMapDF),"\n")
   
   #cat("Calling SM_Build\n")
   xRes  <- SeerMapper::SM_Build(rPM)    # problem SeerMapper can't load data... when loaded in this manner.
   #str(xRes)
   #cat("Return from SM_Build\n")
   
   #cat("length of dataMapDF after build:",dim(rPM$dataMapDF),"\n")
   
   rPM   <- xRes$rPM
   MV    <- xRes$MV
   
   dataMapDF <- rPM$dataMapDF   # restore dataMapDF is changed in SM_Build
   
   #cat("completed SM_Build Z-3574 .\n")
   
   idList <- dataMapDF$ID   # get updated ID Lists after SM_Build edit.
   lenIdList <- length(idList)
   
   #cat("idMode:",rPM$idMode,"\n")
   #cat("len(st99):",length(MV$st99_d00),"\n")
   #cat("names(MV):",names(MV),"\n")
   
   #
   #  Have all boundaries loaded and selected to meet the map needs.
   #  This is the same for all maps to be drawn, only the 
   #  dataMapDF$data changes (the color)
   #  dataMapDF is indexed by location ID
   #
   #  If Location IDs are not found, their row in dataMapDF was deleted.
   #  We plot from dataMapDF pulling data from the other satscan tables.
   
   states_data    <- MV$states_data
   states_proj    <- MV$states_proj
   
   # pull over "scale" variable - needed to scale outline for states like Alaska.
   xM               <- match(st99.index$stID,states_data$stID)
   st99.index$scale[!is.na(xM)] <- states_data$scale[xM][!is.na(xM)]   # save scale value
   st99.index$moveX[!is.na(xM)] <- states_data$moveX[xM][!is.na(xM)]   # save X move
   st99.index$moveY[!is.na(xM)] <- states_data$moveY[xM][!is.na(xM)]   # save Y move
   
   xM               <- match(co99.index$stID,StateListDAll)
   #  must be stID match and in the right census year.
   xMGood           <- !is.na(xM) & ( bitwAnd(co99.index$y,yF)==yF )  # exist in state list and cur census year.
      
   co99.index       <- co99.index[xMGood,]  # keep only counties for census year and states with data.
   #coSelList        <- co99.index$ID
      
   xM               <- match(pl99.index$stID,StateListDAll)
   xMGood           <- !is.na(xM) 
   pl99.index       <- pl99.index[xMGood,]   # keep places for states with data.
   #plSelList        <- pl99.index$ID
   
   xM               <- match(tr99.index$stID,StateListDAll)
   xMGood           <- !is.na(xM) & ( bitwAnd(tr99.index$y,yF)==yF )
   tr99.index       <- tr99.index[xMGood,]   # only keep rows with states matching the data.
   #trSelList        <- tr99.index$ID
   
   #
   # differences between states_data and st99_data 
   #                    states_data  		st99_data -> st99.index
   #        row.names   -   Y        		  Y
   #	    ID          -   Y			  Added
   #        stID        -   Y   		  Added
   #        abbr        -   Y        		  Y
   #        stName      -   Y          		  Y
   #        rgID        -   Y          		  Y
   #        rgName      -   Y          		  Y
   #        dvID        -   Y          		  Y
   #        dvName      -   Y          		  Y
   #        loc         -   Y			  Y
   #        DoAdj       -   Y          		  Yes
   #        moveX       -   Y                     Added
   #        moveY       -   Y                     Added
   #        scale       -   Y   		  Added      
   #        proj        -   Y         		  Yes
   #        hsa         -   hsa00, hsa10          hsa_00, hsa_10   # ******
   #        county      -   county00, county10	  county_00, county_10
   #        tracts      -   tracts00, tracts10	  tracts_00, tracts_10
   #        change10    -   Y			  Y
   #        c_X         -   Y         		  Y 
   #        c_Y         -   Y         		  Y 
   #                               c_X and c_Y can be computed from states_proj
   #                                cXY <- coordinates(states_proj)
   #                                colnames(cXY) <- c("c_X", "c_Y")
   #                                states_data <- cbind(states_data,cXY)
   #
   #        use         -              ?    internal
   #        USAPlot     -              us48Only similar
   #
   #     
   #   comparison of        co99_data  and  co99.ind00
   #            coName          Y             Y 
   #            tracts          Y            empty
   #            ID              Y            row.names
   #            stID            Y             Y
   #            stcoID          Y           = ID
   #            stName                        Y     (can get from states_data)
   #            stcotrID        Y*NA         
   #            saID            Y*???
   #            c_X and c_Y                   Y     (can build from co99_proj)
   #
   #
   #    comparison of       tr99_data and tr99.ind00
   #      ID/row.names         ID/rn         rn
   #        stID                 Y           Y
   #        stcoID               Y          stcoID 
   #        stcotrID             Y           rn    (same as ID and row.name)
   #        saID                 Y           ?     (can get from county_data)
   #        stName                           Y     (can get from states_data)
   #        coName                           Y     (can get from county_data)
   #        c_X & c_Y            *           Y     (can get from tr99_proj)
   #        Place                            Y     (can't get)
   #
   #   comparison of pl99.ind00  ???
   #     problem - places assigned to tracts at time of build.
   #     now have a co99 for counties.   But no tr99...
   #     looks like minimal tract information needed that can't be
   #     rebuilt is "place".
   #
   #   Except for the centroid values (X, Y), coName, Place Name, and 
   #   Seer Registry index, it looks like the NCI_Ind00 set
   #   of indexes (st99.ind00, co99.ind00, tr99.ind00 could 
   #   be replaced by the states_data, county_data, tract_data
   #   structures from SeerMapper.   satscanMapper could carry
   #   parallel items to these, integrate the info into SeerMapper,
   #   or keep all index inform in satscanMapper and not need 
   #   SeerMapper.   Issue is keeping them synchronized and 
   #   the ability to do 2000 and 2010 census.
   #
   
   #
        
   ####
   #
   #   Collect the low level Location ID List.  (What we see in the data.)
   #

   #   Based on the Cluster Information, get the list of Years to map.
   #

   YearList    <- sort(unique(ColYrList))
   YrList      <- YearList
   #cat("Years to be plotted:",YearList,"\n")

   #
   ####

   #####
   ##   Preformed by SM_Build in SeerMapper
   ##
   ##   Validate the LOC_IDs against the boundary information for identificed states
   ##
   ##   This verification is done in SeerMapper = SM_Build - should it be done here??
   ##
   ##
   #LocFipsAll   <- str_sub(tr99.index$ID,1,LocIDType)  # get list of all valid fips at the correct level
   #LocFipsAll   <- unique(LocFipsAll)
   #
   #fipsMatch    <- match(LocListDAll,LocFipsAll)
   #fipsMatchNA  <- is.na(fipsMatch)    # true = no match
   #
   #if (any(fipsMatchNA)) {
   #   BadList    <- LocListDAll[fipsMatchNA]
   #   xmsg       <- paste0("The LOC_IDs contain FIPS codes that do not exist in the boundary info.\n",
   #                          "They are (",length(BadList),"):\n")
   #   warning(xmsg,call.=FALSE)
   #   xmsg       <- paste0(BadList,collapse=", ")
   #   warning(xmsg,call.=FALSE)                       
   #   xmsg       <- paste0("Please review the location IDs used and correct. Verify the right census year was specified.\n")
   #   warning(xmsg,call.=FALSE)
   #   FNDError   <- TRUE
   #   NumErrors  <- NumErrors + 1
   #   # don't shut down - ignore area.
   #}
   #
   #####
   #
   #if (FNDError) {
   #   xmsg     <- paste0("Error(s) found in the LOC_ID values in the results.  Run Terminated.\n")
   #   stop(xmsg, call.=FALSE)
   #}
   #
   ##
   ######
   
   #####
   #
   #   Identify level of data being used for the location ID
   #
   #####
   
 
   #####
   #
   #  Identify te states covered by the data.
   #
   #   From the LocListDAll, get the list of States involved.
   #

   #cat("List of states involved in the analysis Z-3282 :",StateListDAll,"\n")
   
   st99.index$use <- FALSE
   st99.index[StateListDAll,]$use <- TRUE     # mark used states in st99.index table
   
   #####
   #
   #  get centroids for the right level and census year.
   #        st99_data, co99_data, tr99_data contain mapping centroids.
   #        st99_M_data, co99_M_data, tr99_M_data contain the adjusted calculation centroids for 
   #        the states that were moved on the visual map - Alaska, Hawaii, PR  
   #        these files are only used by CreateGeo4SS function.
   #
   # county and/or tract data
   cyXY  <- c("c_X_00","c_Y_00")  # centroids
   cyHs  <- c("hsa_00")           # counties
   cyCo  <- c("county_00")        # counties
   cyTr  <- c("tracts_00")        # tracts

   if (censusYear == "2010")  {
      cyXY  <- c("c_X_10","c_Y_10")
      cyHs  <- c("hsa_10")
      cyCo  <- c("county_10")
      cyTr  <- c("tracts_10")
   }

   if (LocIDType == 2) {  # if state data..
      # state level
      cyXY  <- c("c_X","c_Y")    # same for both census years
   }

   #cat("cyXY:",cyXY,"\n")  # list the labels
   
   #####
   #
   #  The index are loaded and provide enough information to process the data and verify it.
   #  Only the centroids are needed for the circle or ellipse overlays with labels.
   #
   #####
   
   #st99.index$hsa     <- st99.index[,cyHs]     # pick up state hsa count for year
   st99.index$county   <- st99.index[,cyCo]     # pick up state county count for year
   st99.index$tracts   <- st99.index[,cyTr]     # pick up state tract  count for year
   
   #hs99.index$county  <- hs99.index[,cyCo]     # pick up hsa county count for year
   #hs99.index$tracts  <- hs99.index[,cyTr]     # pick up hsa tract count for year
   
   co99.index$tracts   <- co99.index[,cyTr]     # pick up county tract count for year

   pl99.index$tracts   <- pl99.index[,cyTr]     # pick up place tract count for year
   
   #hs99.index         <- hs99.index[hs99.index$county != 0,]   # remove entries with no counties
   #hs99.index         <- hs99.index[hs99.index$tracts != 0,]   # remove entries with no tracts

   co99.index          <- co99.index[co99.index$tracts != 0,]   # remove entires with no tracts
   pl99.index          <- pl99.index[pl99.index$tracts != 0,]
   
   us99.index          <- NULL
   us99.index$states   <- dim(st99.index)[1]
   #us99.index$hsa      <- sum(st99.index$hsa)
   us99.index$county   <- sum(st99.index$county)
   us99.index$tracts   <- sum(st99.index$tracts)
   
   
   #cat("us99.index:\n")
   #print(us99.index)
   
   #####
   #
   #  The package can now handle State, County and Census Tract location IDs.
   #  The original package handled tracts, only 
   #  
   
   #
   #  Build working index table for the approviate lavel.  
   #    Key items are: names, centroids, and link to SeerMapper mapping.
   #    In our case, the IDs will always be the Fips codes (2, 5, and 11 digits)
   #

   #
   #  At the state level, we have state names, abbreviations and centroids
   #  At the county level, we have state names, abbr. and county names, abbrs.
   #  At the census tract level, we have state names/abbr., county names/abbrs,
   #      placenames, and tract IDs.
   #
   #  resulting data.frame has:
   #      Name of area
   #      ID of area
   #      stID associated with area
   #      stcoID associated with area
   #      plID associated with area
   #      centroid X, Y coordinates
   #
   #  SeerMapper Call Parameter
   #
   vRegionB <- "NONE"
   vStateB  <- "NONE"
   vSeerB   <- "NONE"
   vHsaB    <- "NONE"
   vCountyB <- "NONE"
   vTractB  <- "NONE"
   
   if (outline) {
      vDataB  <- "DATA"
   } else {
      vDataB  <- "NONE"
   }
   
   vDataBCol <- bndyCol 
   
   # Init and Collect working data.frame (wk99)
   wk99 <- NULL
   
   #cat("LocIDType Z-3841 :",LocIDType,"  idList:\n")
   #print(head(idList))
   
   if (LocIDType == 2) {
      # state level data 
      #cat("getting state level wk99 table with centroids.\n")
      fstr           <- c("ID","stID","stName","abbr","county","tracts",cyXY)
      fstr2          <- c("ID","stID","stName","stAbbr","county","tracts","t_X","t_Y")
      wk99           <- st99.index[,fstr]    # all states.  NO difference from 2000 and 2010
      colnames(wk99) <- fstr2
      
      #  Only keep the states with data.
      xM             <- match(wk99$stID,StateListDAll)
      xKeep          <- !is.na(xM)
      wk99           <- wk99[xKeep,]
      
      wk99$Name      <- wk99$stName
  
      wk99$stcoID    <- NA
      wk99$coName    <- NA
      wk99$stcotrID  <- NA
      wk99$plName    <- NA
  
      IDName         <- "States"
      rPM$idMode     <- 1
      
      vRegionB       <- "DATA"
      vStateB        <- vDataB
      
   }
   if (LocIDType == 3) {
      cat("getting HSA level wk99 table with centroids.\n")
      # HSA level data
      #cat("colnames(hs99.index):",colnames(hs99.index),"\n")
      #
      #  Complete Coding 
      #
   
   }
   if (LocIDType == 5) {
      #cat("getting county level wk99 table with centroids.\n")
      # county level data
      #cat("colnames(co99.index):",colnames(co99.index),"\n")
      
      fstr           <- c("ID","stID","stName","stcoID","coName","tracts",cyXY)
      fstr2          <- c("ID","stID","stName","stcoID","coName","tracts","t_X","t_Y")
      wk99           <- co99.index[,fstr]   # all counties
      colnames(wk99) <- fstr2
   
      wk99$stAbbr    <- st99.index[wk99$stID,"abbr"]
      wk99$Name      <- wk99$coName
      wk99$stcoID    <- wk99$ID
      wk99$county    <- 1
      
      wk99$stcotrID  <- NA
      wk99$plName    <- NA
   
      IDName         <- "Counties"
      rPM$idMode     <- 2
      
      vStateB        <- "ALL"
      vCountyB       <- vDataB
      if (vCountyB != "NONE")  vCountyB <- "STATE"
      
   }
   if (LocIDType == 11) {
      #cat("getting tract level wk99 table with centroids.\n")
      # census tract level data
      fstr           <- c("ID","stID","stcoID","stName","coName","plName",cyXY)
      fstr2          <- c("ID","stID","stcoID","stName","coName","plName","t_X","t_Y")
      wk99           <- tr99.index[,fstr]
      colnames(wk99) <- fstr2
   
      wk99$Name      <- wk99$ID
      wk99$stcotrID  <- wk99$ID
      wk99$tracts    <- 1
      wk99$county    <- NA
      
      xM             <- match(wk99$stID,st99.index$stID)
      wk99$stAbbr    <- st99.index[xM,"abbr"]    # back fill state abbrevation.
   
      IDName         <- "Census Tracts"
      rPM$idMode     <- 3
   
      vCountyB       <- "DATA"
      vTractB        <- vDataB
      if (vDataB != "NONE")  vTractB <- "COUNTY"
  
   }
   
   #cat("initial wk99 table Z-3525 build from boundary tables not data.\n")
   #print(str(wk99))
   #print(head(wk99,15))
   #cat("typeof:",typeof(wk99),"  class:",class(wk99),"\n")
   
   
   idMode           <- rPM$idMode
   
   ##
   # save values in rPM list
   # 
   rPM$regionB      <- vRegionB
   rPM$stateB       <- vStateB
   rPM$seerB        <- vSeerB
   rPM$hsaB         <- vHsaB
   rPM$countyB      <- vCountyB
   rPM$tractB       <- vTractB
   
   dataMapDF$stAbbr <- wk99[dataMapDF$ID,"stAbbr"]
   dataMapDF$col    <- wk99[dataMapDF$ID,"LocCol"]
   dataMapDF$cat    <- wk99[dataMapDF$ID,"LocCat"]
   
   if (idMode == 3) {
      # census tract  - extra fields
      dataMapDF$plKey  <- tr99.index[dataMapDF$ID,"plKey"]
      dataMapDF$plName <- str_sub(dataMapDF$plKey,6,-1)
      
   }
   
   rPM$dataMapDF    <- dataMapDF
   
   #cat("dataMapDF after wk99 initialization.\n")
   #print(str(dataMapDF))
   #print(head(dataMapDF))
 
   ####
   #
   #  Last adjustments of X,Y coordinates - Get new center coordinates 
   #
   #    Build label points (centroids) of the area (state, county, or census tract) 
   #    based on the centroid information the package has that matches the 
   #    associate boundry data.  
   #
   #  The wk99 data.frame is setup based on the IDType and LOC_IDs.  This is done in case
   #  the centroids in the SaTScan (TM) results are not based on the SeerMapper boundary data.
   #  This links the location IDs to the correct centroids for the mapping.
   #
   #  If cartesian, calculate ratio to adjust circle and ellipses
   #  
   #  Used later with ColTable and GisTable  
   #
   # scaling range
   # 
   #  The user data is either X, Y (cartesian of some metric) or LATITUDE and LONGITUDE
   #  if Lat/Long.  
   #  Circles have RADIUS;  Ellipse have Min/Maj Axis, angle, shape.
   #    units => lat/long - circle - radius is "km".  (we use "m")
   #             castesian - circle/elliptical - units of projection
   #
   #  shape = circular = lat/long used - radius is KM and is divided by 1000 and used.
   #  shape = circular = castesian   - range of x,y coordinates before and after
   #          are used to get ratio
   #  shape = ellipical = cartesian - range of x,y coordinates before and after
   #          are used to get ratio
   #
   #
   
   # New X/Y centroid   - RISK  ->   LOC_ID does not exist.
   ColTable$t_X <- wk99[ColTable$LOC_ID,"t_X"]
   ColTable$t_Y <- wk99[ColTable$LOC_ID,"t_Y"]
   
   GisTable$t_X <- wk99[GisTable$LOC_ID,"t_X"]
   GisTable$t_Y <- wk99[GisTable$LOC_ID,"t_Y"]

   XRangeT   <- range(ColTable$t_X,GisTable$t_X)
   YRangeT   <- range(ColTable$t_Y,GisTable$t_Y)
   XDiffT    <- diff(XRangeT)
   YDiffT    <- diff(YRangeT)
   
   # see if we need to calculate a ratio..
   
   if (CoordinatesType == 0) {   # cartesian (circle or ellipical)

      xM <- match ("X", ColTableNames)  # is "X" present
      if (is.na(xM)) { 
         # no "X"
         stop("Cartesian Coordinates specified, but X/Y values missing from Cluster table.")
      } else {
         # Original X/Y centroid
         OrigXRange    <- range(ColTable$X,GisTable$LOC_X)
         OrigXDiff     <- diff(OrigXRange)
         OrigYRange    <- range(ColTable$Y,GisTable$LOC_Y)
         OrigYDiff     <- diff(OrigYRange)
   
         ColXRatio    <- XDiffT/OrigXDiff
         ColYRatio    <- YDiffT/OrigYDiff
         ColRatio     <- mean(ColXRatio,ColYRatio)
 
         #cat("Col - xDiff:",OrigXDiff,"  xDiff2:",XDiffT,"  XRatio:",ColXRatio,"\n")
         #cat("Col - yDiff:",OrigYDiff,"  yDiff2:",YDiffT,"  YRatio:",ColYRatio,"\n")

         # adjust ellipse or circle metrics      
         ColTable$E_MAJOR <- ColTable$E_MAJOR * ColRatio
         ColTable$E_MINOR <- ColTable$E_MINOR * ColRatio
         
      }
   
   } else {
      # Lat/Long    - circle only.
         ColRatio         <- 1000    #   km to m    
         # adjust ellipse or circle metrics      
         ColTable$E_MAJOR <- ColTable$E_MAJOR * ColRatio
         ColTable$E_MINOR <- ColTable$E_MINOR * ColRatio
   }
 
   #
   ####
   #cat("Mapping Level Z-4065 - IDName:",IDName,"\n")   
   #cat("SeerMapper Parameters - regionB:",vRegionB,"  stateB:",vStateB,"  seerB:",vSeerB,"  countyB:",vCountyB,"  tractB:",vTractB,"\n")
   
   ####
   #
   # Check the rough x,y box needed by the plot to see what paper to request 
   #   on the PDF.
   #
   #  REDO - must be a better way to set page, width and height - mis-using most 
   #    of the space on the page.
   #
   #wWide   <- diff(range(wk99$t_X))
   #wHeight <- diff(range(wk99$t_Y))
   #aspect  <- wHeight/wWide
   aspect   <- XDiffT/YDiffT   #  width over height
   
   # default
   vW       <- 7.5
   vH       <- 10
   if (aspect > 0.1 && aspect < 0.577) {
      vW <- 7.5
      vH <- 13
   }
   if (aspect >= 0.577 && aspect <= 1) {
      vW <- 7.5
      vH <- 10
   }
   if (aspect > 1 && aspect <= 1.733) {
      vW <- 10
      vH <- 7.5
   }
   if (aspect > 1.733 ) {
      vW <- 13
      vH <- 7.5
   }
   #cat("aspect:",aspect,"  vW:",vW,"  vH:",vH,"\n")
   
   
   
   #wPaper  <- "letter"
   #if (aspect>1.8) {
   #   wPaper <- "legal"   # vertical long page
   #}
   #if (aspect>1.4) {
   #   wPaper <- "letter"  # vertical page
   #}
   #if (aspect<0.72) {
   #   wPaper <- "USr"     # letter rotated landscape
   #}
   #if (aspect<0.55) {     # legal rotated.
   #  wPaper <- "USr"
   #}
   #cat("wPaper Z-3632 size:",wPaper,"  aspect:",aspect,"\n")
   
   #
   #
   ####

   #########################Start Here #####
   
   ####
   #
   #  Main loop to plot the cluster maps.   Single pass if Spatial - Multiple pass by year if Space/Time
   #
   ####
   #cat("Main loop to plot clusters per year.\n")
   
   ####
   ##
   ##  Next step is to plot the patterns and clusters
   ##  (Spatial Only and one per month for the Space/Time results)
   ##
   ####
   #
   
   #  Create output file names for graphics and report
   
   YrList <- sort(unique(ColYrList))              # how to do loop - NULL year = "0000" 
   YearList <- ""

   if (SSForm == 0) {  
      # for Spatial Only = year is zero 
      TypeForm <- "Spatial Only Report"
      #
   }
   
   if (SSForm == 1) {  # Space/Time Report  
      #
      #  Total by year -- then plot maps.   (Pop and Cas data is tagged with year of collection.
      #
      
      #  Setup report line of list of years.
  
      YearText <- paste0(YrList,collapse=" ")
  
      YearList <- paste0("For the following Years: ",YearText)
      #
      TypeForm <- "Space / Time Report"
      #
   }
   #cat("Type of Report:",TypeForm,"\n")
   
   #  Report the output filenames - Graphics and Report
   
   #cat("Output PDF and Report Filenames:\n")
   #cat("   OutputFNpdf (graphics):",OutputFNpdf,"\n",sep="")
   #cat("   OutputFNtxt (report)  :",OutputFNtxt,"\n",sep="")
   
   ####
   #
   #  Verify that any census tract referenced in the *.col.dbf, and *.gis.dbf 
   #  tables exists in the boundary data census tract tables.
   #  If one is missing, generate warning to user.
   #
   # With the new logic - we have collected the LOC_ID from the files to determine 
   # what area we are mapping.  Then pull the state, county and tract as required.
   #
   # We can use this list to check if the boundaries exists in the files and report
   # any issues.
   #
   # wk99 index DF is only for the states in the StateList which was derived from 
   #  the IDs in the cluster, Location and RR results files.  The other county or tract
   #  information is not needed and ignored.
   # wk99 based on boundary data.
   # LocListDAll based on data in results.
   # dataMapDF based on data in results reduced is no match LOC_ID.
   #
   
   # Bad LOC_ID in gis - ignore
   # Bad LOC_ID in col - (centroid area)
   
   xx <- setdiff(ColTable$LOC_ID,wk99$ID)  # find items in cluster loc_id list that are not in xx99.index
   #cat("compare cluster ID lists Z-3760 - ",length(ColTable$LOC_ID), " vs. ", length(wk99$ID),"\n")
   
   # should be identical - since wk99$ID was derived from the data and LocListDAll ... 
   
   if (length(xx) > 0) {  # will not be NULL - length = 0 if nothing left out.
   
      xmsg    <- paste0("SaTScan cluster results data contains one or more ",IDName," LOC_IDs that are not in the boundary database. Unable to map results.")
      warning(xmsg, call.=FALSE)
      cat("SaTScan and rate data LOCATION IDs not found:\n")
      print(xx)
      stop("Fix and rerun - Execution stopped.")
   }
   rm(xx)
   
   xM             <- match(GisTable$LOC_ID,wk99$ID)   # find items in location (gis) loc_id list that are ot in xx99.index
   #cat("compare location ID lists Z-3727 - ",length(GisTable$LOC_ID)," vs. ",length(wk99$ID),"\n")
   
   # should be identical - since wk99$ID was derived from the data and LocListDAll ... 
   
   if (any(is.na(xM))) {  # NA means not found it boundary list.
   
      xmsg    <- paste0("SaTScan location(GIS) results data contains one or more ",IDName," LOC_IDs that are not in the boundary database. Locations will not be mapped.")
      warning(xmsg, call.=FALSE)
      cat("SaTScan and rate data LOCATION IDs not found:\n")
      xx       <- GisTable$LOC_ID[is.na(xM)]
      print(xx)
      GisTable <- GisTable[!is.na(xM)]   # keep good entries.   Drop entries without boundaries.
      rm(xx)
   }
  
   
   #
   #  Reverse GisTable so lowest cluster seen first for mapping.
   #
   xGord        <- order(GisTable$CLUSTER,decreasing=TRUE)
   GisTable$Seq <- as.integer(seq(1:length(GisTable$CLUSTER)))       # save original sequence number
   GisTable     <- GisTable[xGord,]   # reorder GisTable by cluster.
   
   #
   #
   ####
   
   ####
   #
   #  Everything is basically setup - now.
   #
   #  If testing and re-running from here - otherwise change filenames and run from 
   #  the beginning.
   #
   
   #cat("Start of mapping Z-3908 \n")

   ########
   ###
   ###   Have not had time to repeat this process for each year of the
   ###   Space/Time reports.
   ###
  
   ####
   ###
   ##    These plots are not printed at the moment, so I put a sleep statement
   ##    between each one to see it on the screen for a while.
   ##
   ##    If Spatial Only, only one year "0000", so loop executed once.
   ###
   ####
  
   ####
   #
   #   wk99 setup to build the data.frame for SeerMapper call.
   #
   #   idCol    = "ID"
   #   dataCol  = either "LocCol" or "ClusCol"
   #   categ    = "COLORS"
   #   dataBCol = bndyCol
   #   outline  -> impacts the data Level setting to "???" or "NONE"
   #
   #   if state:
   #       region =TRUE
   #       stateB =DATA    or "NONE"
   #       seerB  =NONE
   #
   #   if hsa:
   #       region =FALSE
   #       stateB =NONE
   #       seerB  =NONE
   #       hsaB   =DATA
   #
   #   if county:
   #       region =FALSE
   #       stateB =ALL
   #       seerB  =NONE
   #       hsaB   =NONE
   #       countyB=STATE   or "NONE"
   #
   #   if tract:
   #       region =FALSE
   #       stateB =NONE
   #       seerB  =NONE
   #       hsaB   =NONE
   #       countyB=DATA
   #       tractB =COUNTY  or "NONE"
   #
   
   ####
   ###
   ##
   #
   #     Setup the Ellipse/Circle information for later overlaying the generated maps.
   #
   # Find Ellipses and plot over US Map and location/rate images.
   # Get list of Clusters with P_Value <= PValue (0.05)
   # Build form Cluster Table entries where:  
   #       a) Only use entries with valid P_Values,  b) P_Value < PValue (def = 0.05)
   
   LabelBox  <- matrix(c(range(wk99$t_X),range(wk99$t_Y)),ncol=2,byrow=TRUE)
   colnames(LabelBox) <- c("min","max")
   rownames(LabelBox) <- c("x","y")
   
   xDist    <- diff(LabelBox["x",])
   yDist    <- diff(LabelBox["y",])
   
   xLabOff  <- xDist * 0.005  # Label offset from centroid = 0.5% of width of map centroids.
   
   #cat("xDist:",xDist,"  yDist:",yDist,"  Labelbox-x:",LabelBox["x",],"  y:",LabelBox["y",],"  xLabOff:",xLabOff,"\n")
   
   EUpValue <- FALSE   # indicator if any Ellipses/circle have pValues < 0.05 to map.
  
   if (all(is.na(ColTable$P_VALUE))) {
      # all P_VALUEs in the cluster table are NA.  No comparison - print all.
      EllipseSetBase <- ColTable     
      print("No pValues found in the Cluster record, all are NA.")
   } else {
      EllipseSetBase <- ColTable[!is.na(ColTable$P_VALUE) & ColTable$P_VALUE <= pValue,]   # copy all entries with values and that are < 0.05
      #cat("Found ",dim(EllipseSetBase)[1]," of ",dim(ColTable)[1], " clusters have a pValue <= 0.05 and will be outlined on the map.\n")
      EUpValue <- TRUE  # yes there are values.
   }
   
   ElState                <- str_sub(EllipseSetBase$LOC_ID,1,2)
   EllipseSetBase$scale   <- st99.index[ElState,"scale"]

   EllipseSetBaseLen      <- dim(EllipseSetBase)[1]  # Get number of clusters with P < pValue (0.05)
   
   #cat("Number of ellipse/circles:",EllipseSetBaseLen,"\n")
   #cat("EllipseSetBase:\n")
   #print(EllipseSetBase)
   #cat("EllipseSetBase - typeof:",typeof(EllipseSetBase),"  class:",class(EllipseSetBase)," Z-3987 \n")
   
   #
   # Calculate the X, Y coordinates for the label, and back fill for no circle dimensions
   # for single sub-area cluster.
   #
   
   if (EllipseSetBaseLen > 0) {
   
      for (ine in c(1:EllipseSetBaseLen)) {
         if (EllipseSetBase[ine,"E_MAJOR"] == 0) {
            if (EllipseSetBase[ine,"NUMBER_LOC"] == 1) {
              # no outline parameters
               x <- MV$data_proj[EllipseSetBase[ine,"LOC_ID"],]
               xa <- sum(sapply(slot(x,"polygons"), function(x) x@area))
               xRadius <- sqrt(xa / pi)
               EllipseSetBase[ine,"E_MAJOR"] <- xRadius
               EllipseSetBase[ine,"E_MINOR"] <- xRadius
            } else {
               cat("Found cluster with more than one location with no circle/ellipse values: ",ine,"\n")
            }
         }
      }
            
      EllipseSetBase$E_MAJOR <- EllipseSetBase$E_MAJOR * EllipseSetBase$scale
      EllipseSetBase$E_MINOR <- EllipseSetBase$E_MINOR * EllipseSetBase$scale
     
      #  angle of ellipse - radians ---  circle should be 0 degrees
      EllipseSetBase$RAng    <- (EllipseSetBase$E_ANGLE + 90) * pi / 180    # convert degrees to radians
      
      #  get radius distance 
      RDis <- sqrt ( (EllipseSetBase$E_MAJOR * sin(EllipseSetBase$RAng) ) ^ 2  
                               + (EllipseSetBase$E_MINOR * cos(EllipseSetBase$RAng) ) ^ 2)  # distance to left or right edge of circle/ellispe
   
      #  calculate the label off on x axis (off the centroid.
      EllipseSetBase$Lab_X   <- EllipseSetBase$t_X +  RDis + xLabOff    # center (x) + distance + 0.5% of width = label offset x
       
      # EllipseSetBase$Lab_X <- EllipseSetBase$t_X + (EllipseSetBase$E_MAJOR + EllipseSetBase$E_MINOR)/2 + 12500
      EllipseSetBase$Lab     <- as.character(EllipseSetBase$CLUSTER)
      
      #  Y lab offset is the same as the cY value.
      
      #  Assign color to the outline - default is black, if cluster is high or low, change to highest or lowest color reserved.
      EllipseSetBase$OutCol   <- "black"                                   # default color
      EllipseSetBase$OutCol[EllipseSetBase$HL == 1]  <- ColorsB_Clus_High  # if high, new color
      EllipseSetBase$OutCol[EllipseSetBase$HL == -1] <- ColorsB_Clus_Low   # if low, different color       
      EllipseSetBase$OutCol   <- as.character(EllipseSetBase$OutCol)       # make sure it's character
      
      # $Col already assigne previously  (color of sub-area)
   
      #  Get list of cluster numbers from each table
      #cat("sorting ellipse and giscluster list\n")
      
      EllipseClusList <- sort(unique(EllipseSetBase$CLUSTER))   # actually the Cluster table with pValue < 0.05
      GisClusList     <- sort(GisTable$CLUSTER)                 # location table
     
      #  Ellipse Cluster must have mate in GisClusList
      EllipseCnt <- sapply(EllipseClusList, function(x)  length(GisClusList[GisClusList==x]))  # number of location records per cluster.
 
      if (!all(EllipseSetBase$NUMBER_LOC == EllipseCnt))  {  # number of locations.
         # the number of locations per cluster don't match
         xmsg     <- paste0("The number of GIS entries for each interesting Cluster do not agree. Check SatScan files.")
         warning(xmsg, call.=FALSE)
      }
   } else {
      # no clusters in the results files.
      xmsg     <- paste0("The *.col.dbf file does not contain any clusters with P_VALUE < 0.05.")
      warning(xmsg,immediate=TRUE, call.=FALSE)
      stop(call.=FALSE)
   }

   ###
   
   #cat("EllipseSetBase Z-4400 \n")
   #print(str(EllipseSetBase))
   #print(head(EllipseSetBase))
   
   TL2Cex <- 1

   # Open output file. *** for testing comment out.

   #cat("opening output files...\n")
   
   #
   #pdf(file=OutputFNpdf,paper=wPaper)
   #   set page margins (in lines.)
   
   pdf(file=OutputFNpdf,width=vW,height=vH)
   
   par( omi=c(0.5,0.5,0.5,0.5) )    # outer margin of 0.5 border
   #           B   L   T   R
   rPM$omi <- c(0.5,0.5,0.5,0.5) 
   
   par( mar=c(4.1,2.1,5.1,2.1) )  # default inner margins
   #           B   L   T   R
   rPM$mar <- c(4.1,2.1,4.1,2.1)  # top = 5 lines, bottom = 4 lines , left & right = 2 (not used)

   
   #   Titles
   #     vTitle - user provided MAIN title  (vTitle)
   #     SMTitle - Page title               (vGLine)
   #     R1Title or R2Title - type of data on map Loc or Clu ODE (vPLine) or legends
   #
   vTLine <- 2            # (Line 3)
   vGLine <- 1            # (Line 2)  
   vPLine <- 0.01         # (Line 1)
   if (length(vTitle) > 1) {   # number of user title lines.
      vTLine <- vTLine + 0.5
      vGLine <- vGLine - 0.2
   } 
   #
   #  Legends Page:   vTLine & space (.5 line) & vPLine 
   #
   #  Map pages:      vTLine & vGLine & vPline 
   #     
   #
   ###
   #
   # print Cluster legends on Page # 1
   #
   #cat("Drawing the Legends.\n")
   
   plot.new()

   #   run title at line 1.5; page title at 0.01 line 0
   #
   title(main=vTitle,            line=vTLine-0.5,  cex.main=TL2Cex)     # cluster title (Top Line)
   
   #
   #
   title(main="Mapping Legends", line=vPLine,      cex.main=TL2Cex*0.93)
   
   #
   #  Adjust label off right
   #
   #EllipseSetBase$OffRight <- !( EllipseSetBase$Lab_X < LabelBox[1,2] )
   #
   #RightEdge <- LabelBox[1,2] - strwidth("00")
   #EllipseSetBase[OffRight,"Lab_X"] <- RightEdge
  
  
   #
   #  x,y coordinates are in user units - x range 0-1, y range 1-0
   # 
  
   #    left side
   
   #
   #   Observed/Expected Legends
   #
   #  Cluster Info Table Observed/Expected ratio legend
   xx <- legend(0.1,0.95,ColTableCat,
              title = c("Cluster Obs/Exp Ratio"),
              cex   = 0.75,
              pch   = 22, pt.cex=1.5, pt.bg=ColorsB_Clus_Mid,
              xjust = 0,
              text.width=0.3
         )
   
   topY  <- xx$rect$top - xx$rect$h - 0.05
   
   #  Location Observed/Expected ratio legend
   xx <- legend(0.1,topY,GisTableCat,
              title = c("Location Obs/Exp Ratio"),
              cex   = 0.75,
              pch   = 22, pt.cex=1.5, pt.bg=ColorsB_Clus_Mid,
              xjust = 0,
              text.width=0.3
         )
   
   #print(head(xx,10))
   topY   <- xx$rect$top - xx$rect$h - 0.05  
   #
   #  Cluster Low and High cluster 
   #
   xx <- legend(0.1,topY,c("Low Cluster","High Cluster"),
              title = "Cluster Outline Colors",
              cex   = 0.75, 
              pch   = 22, pt.cex = 1.5, pt.bg =c(ColorsB_Clus_Low,ColorsB_Clus_High),
              text.width = 0.3
         )
   
   #
   #    right side 
   #
   #    upper left
   
   #  By Cluster Classification Key   
   #
   #  How to handle more clusters then Space ?????
   #
   xL <- paste0("Cluster # ",as.character(EllipseSetBase$CLUSTER))
   lenLeg <- length(xL)
   if (lenLeg > 80) {
      xmsg=paste0("Legend page will only list the first 80 clusters.")
      warning(xmsg,call.=FALSE)
      xL <- xL[1:80]
      lenLeg <- 80
   }
   #maxL <- max(strwidth(xL)) 
   nLCol <- 1
   nLCol <- as.integer((lenLeg-1)/40) + 1  #  1-40 > 1 col  41-80 > 2 col...
   
   xx <- legend(0.55, 0.95, paste0("Cluster # ",as.character(EllipseSetBase$CLUSTER)),
              title = "By Cluster Obs/Exp Category",
              cex   = 0.75,
              pch   = 22,   pt.cex=1.5, pt.bg=EllipseSetBase$Col,
              # text.width = maxL,
              ncol  = nLCol
         )
   
   #cat(" res legend:\n")
   #print(xx)
   
   #
   ###
   
   ###
   #
   #  page 2 and forward.
   #
   #cat("page 2 - maps Z-4565 \n")

   data_data_sel <- MV$data_data_sel
   wk99IDs       <- wk99$ID
   #cat("YrList:",YrList,"\n")
   
   for (WorkYr in YrList) {
   
      #cat("Processing Year:",WorkYr,"\n")
      
      #
      # Set up the common variables for this YEAR for Cluster 
      #
      #  Set up for the type of mapping - Spatial or Space/Time
      #  Place in work variable ALL or only data for specific YEAR
      #
      #  For spatial only YrList = "0000", only one logical year.
      #
      #  Get ellipical data for the year.
      #  Set up titles and pull the ellipical/circle data 
      #
      if (SSForm == 0) {  # Spatial Only 
         #cat("Spatial Only.\n")
      
         SMTitle    <- paste0( "Cluster(s) and ",IDName," - Spatial Only")
             
         # one pass - if no years or one, then all of the rate table is used.
         #            if multiple years, only the last year is used. 
           
         wEllipseSet   <- EllipseSetBase       # spatial only - all of the cluster 
         CurYear       <- " "                  # for title.
         # Ignore year stuff.         
      }
      if (SSForm == 1) { # Space/Time  - by year 
         #cat("Space-Time by year.\n")
         #       
         CurYear       <- as.character(WorkYr)
         SMTitle       <- paste0("Cluster(s) and ",IDName," during ",CurYear)
         
         wEllipseSet   <- EllipseSetBase[(WorkYr >= EllipseSetBase$SYr & WorkYr <= EllipseSetBase$EYr) ,]
      
      }
       
      wEllipseSetLen   <- dim(wEllipseSet)[1]
      
      # End of Setup 
      #cat(SMTitle,"\n")
      rPM$mTitle <- SMTitle
         
      #cat("Number of clusters found in this group/year:",wEllipseSetLen,"\n")
      # cat("wEllipseSet:\n")
      #print(head(wEllipseSet))
      
      ####
      #
      #  Now map just the cluster data - outline of cluster and location points within 
      #   the cluster - color coded ("n" levels) from low to high or high to low based
      #   on whether the cluster is ranked a high or low cluster (O/E ratio > or < 1).
      #
      ####
     
      ####
      #
      #  Two forms of cluster maps - census tracts colored by local O/E  and
      #  census tracts colored by the Cluster's O/E.
      #
      #  Allow both forms as options.  leave code for original form.
      #
      
      rPM$categMode <- 4
           
      #
      # Start of Cluster Maps - for current year
      #
      #  Have to change the logic.  Before we mapped by cluster - boundaries and overlay.
      #  then on to the next cluster.
      #
      #  Now we have to do all the area mapping and coloring before drawing the ellipses.
      #
     
      #  Do this for each year.
      
      #  clear mapping colors to white
      data_data_sel$col     <- "white"
      
      #
      # Reset the LocCat, LocCol, CluCat, CluCol and CluNum for each area.
      #
      #cat("Setting Cat and Col into wk99 loc and clust.\n")
      
      wk99$LocCat  <- NA
      wk99$LocCol  <- "white"
      wk99$ClusCat <- NA
      wk99$ClusCol <- "white"
      wk99$ClusNum <- NA
      
      # Get working Cluster (Col) and Location (Gis) table for the active year - only.
      
      wGisTable  <- GisTable[(WorkYr >= GisTable$SYr & WorkYr <= GisTable$EYr),]
      wColTable  <- ColTable[(WorkYr >= ColTable$SYr & WorkYr <= ColTable$EYr),]
      
      # Fill in the wk99 working data.frame for mappinng.
      wk99[wGisTable$LOC_ID,"ClusNum"] <- wGisTable$CLUSTER
      wk99[wGisTable$LOC_ID,"ClusCat"] <- wGisTable$ClusCat
      wk99[wGisTable$LOC_ID,"ClusCol"] <- as.character(wGisTable$ClusCol)
      wk99[wGisTable$LOC_ID,"LocCat"]  <- wGisTable$LocCat
      wk99[wGisTable$LOC_ID,"LocCol"]  <- as.character(wGisTable$LocCol)
      wk99[wGisTable$LOC_ID,"SYr"]     <- wGisTable$SYr
      wk99[wGisTable$LOC_ID,"EYr"]     <- wGisTable$EYr
           
           
      xM     <- match(wk99$ClusNum,wColTable$CLUSTER)   # get cluster dates.
   
      wk99$ClusSYr  <- wColTable$SYr[xM]
      wk99$ClusEYr  <- wColTable$EYr[xM]
      
      #  final adjustments
      
      #cat("table cat and col Z-4113 - cluster  dim(wk99):",dim(wk99),"\n")
      #print(table(wk99$ClusCat))
      #print(table(wk99$ClusCol))
      #cat("table cat and col - local.\n")
      #print(table(wk99$LocCat))
      #print(table(wk99$LocCol))
  
      #cat("Col and Cat set.\n")
  
      #rPM$dataMapDF$ID  <- wk99$ID    # NO NO - wk99 is derived from the data and dataMapDF.
      
      #cat("wk99 before mapping.\n")
      #print(str(wk99))
      #print(head(wk99))
      
      #cat("dataMapDF before mapping.\n")
      #print(str(dataMapDF))
      #print(head(dataMapDF))
      
      #  dataMapDF represents the locations in the Cluster and Location tables
      #  wk99      represents the total area to be mapped (at same level be may have more areas
      #  If dataMapDF areas don't map to boundary - they are deleted from dataMapDF
      #  There may not be a match from wk99 to dataMapDF due to: bad loc ID, 
      #        or wk99 contain more areas then the data.
      #  We are mapping dataMapDF with boundary options to fill in rest of region, state, seer.
      #  There should always be an entry in wk99 for an entry in dataMapDF.
      #
      #  If at state level  - fill into regions (regionB=DATA, stateB=REGION)
      #  If at county level - fill into state (regionB=NONE, stateB=DATA, countyB=STATE
      #  If at tract level  - fill into state (regionB=NONE, stateB=DATA, countyB=STATE, tractB=STATE
      #

      idList2             <- dataMapDF$ID
      wk99IDs             <- wk99$ID

      
      if (locO_EMap) {
         #cat("locO_EMap generation.\n")
         #
         #  Now we do the Location ODE 
         #
         R1Title        <- c("Location Obs/Exp Ratio Map")
         
         #cat("R1Title:",R1Title,"\n")
         
         # Next page.  SeerMapper does plot.new() for new page
         
         if (debugFlag) {
            ListColors <-  sort(unique(wk99$LocCol))
            cat("ListColors:",ListColors,"\n")
            cat("vDataBCol :",vDataBCol,"\n")
            cat("R1Title   :",R1Title,  "\n")
            cat("regionB   :",regionB,  "\n")
            cat("stateB    :",vStateB,  "\n")
            cat("seerB     :",vSeerB,   "\n")
            cat("hsaB      :",vHsaB,    "\n")
            cat("countyB   :",vCountyB, "\n")
            cat("tractB    :",vTractB,  "\n")
            cat("LocIDType :",LocIDType,"\n")
            cat("vTitle    :",vTitle,   "\n")
         }
         #cat("colnames(wk99):",colnames(wk99),"\n")
         #print(str(wk99))
         #print(head(wk99,10))
         
         #cat("colors from wk99:\n")
         #print(table(wk99$LocCol))
         #print(table(wk99$ClusCol))
         #print(wk99[,c("ID","LocCol","ClusCol")])
         
         #cat("len dataMapDF:",dim(dataMapDF)[1],"   len wk99:",dim(wk99)[1],"\n")
   
         #cat("assign LocCol to dataMapDF.\n")
         dataMapDF[idList2,"data"]    <- wk99[idList2,"LocCol"]         # set data to wk99 local color.
         
         #print(str(dataMapDF))
	 #print(head(dataMapDF,10))
	 #print(table(dataMapDF$data))
         
         #cat("assign colors to data_data_sel.\n")
         data_data_sel[wk99IDs,"col"] <- wk99[wk99IDs,"LocCol"]       # full geo data list.
         
         #print(str(data_data_sel))
         #print(head(data_data_sel,10))
         
         #cat("saving data into rPM.\n")
         rPM$dataMapDF       <- dataMapDF
         MV$data_data_sel    <- data_data_sel
         
         #cat("Call SM_Mapper for locO_EMap.\n")
         
         resBBox <- SeerMapper::SM_Mapper(rPM,MV)
         
         #cat("Mapping X-Y box:\n")
         #print(resBBox)
         
         if (outline) {
             #cat("doing overlay shapes\n")
             PlotClusterOutline(wEllipseSetLen, wEllipseSet, vLabel, ClusLabCol)
         }
         
         title(main=vTitle,line=vTLine,cex.main=TL2Cex)     # cluster title (Top Line        
         title(main=SMTitle,line=vGLine,cex.main=TL2Cex*0.9)
         title(main=R1Title,line=vPLine,cex.main=TL2Cex*0.9)
      }
      
      #
      #  draw outlines.
      #
       
      if (clusO_EMap) {
         #cat("clusO_EMap drawing.\n")
         #
         #  Now we do the cluster ODE 
         #
         R2Title    <- c("Cluster Obs/Exp Ratio Map")
         rPM$mTitle <- R2Title
   
         data_data_sel$col  <- "white"
         data_data_sel[wk99IDs,"col"] <- wk99[wk99IDs,"ClusCol"]
         
         dataMapDF[idList2,"data"]    <- wk99[idList2,"ClusCol"]
         
         #cat("Set cluster colors into dataMapDF$data.\n")
         #cat("dim(dataMapDF):  ",dim(dataMapDF),      "  dim(wk99):  ",dim(wk99),"\n")
         #cat("names(dataMapDF):",names(dataMapDF),"\n")
         #cat("     names(wk99):",names(wk99),"\n")
         #cat("table(dataMapDF$data):\n")
         #print(table(dataMapDF$data))
        
         MV$data_data_sel    <- data_data_sel
         rPM$dataMapDF       <- dataMapDF
         
         # Next page.  SeerMapper does plot.new() for new page
   
         #cat("Call SM_Mapper for clusO_EMap.\n")
         
         resBBox <- SeerMapper::SM_Mapper(rPM,MV)
         
         #print(resBBox)
         
         if (outline) {
            #cat("drawing overlays on clusO_Emap.\n")
            PlotClusterOutline(wEllipseSetLen, wEllipseSet, vLabel, ClusLabCol)
         }
   
         
         title(main=vTitle,line=vTLine,cex.main=TL2Cex)     # cluster title (Top Line)
         title(main=SMTitle,line=vGLine,cex.main=TL2Cex*0.9)
         title(main=R2Title,line=vPLine,cex.main=TL2Cex*0.9)
  

      } # end of cluster maps for this year.
        
      #print("Completed Cluster map for year.")
 
   }   # end of YEAR Loop FOR

   #  End of main Loop
   
   dev.off()    # close PDF file..
   
   ############################
   ############################

   #print("Starting text report")

   #
   # If call TextRepHdr_Files(TxtCon,TypeForm,YearList,PopInput,CasInput,ColInput,GisInput, RateFile) 
   #
   
   TractN <- "tracts_00"
   CountyN<- "county_00"
   HsaN   <- "hsa_00"
   if (censusYear == "2010") {
      TractN <- "tracts_10"
      CountyN<- "county_10"
      HsaN   <- "hsa_10"
   }
   
   DataTypeLit <- switch(as.character(LocIDType),
                   "2" = "State Location IDs",
                   "3" = "State/HSA Location IDs",
                   "5" = "State/County Location IDs",
                   "11" = "State/County/Census Tract Location IDs",
                   "Unknown Location IDs"
               )
   
  
   
   ########
   ###
   ###   Repeat loops to generate a text report of the detailed data
   ###   used in the mapping above.
   ###

   # Open output text file.
   #cat("Text: Overview Parameters\n")

   TxtCon <- file(OutputFNtxt,"w")    # output to text file
   
   #TxtCon <- stdout()     # output to screen instead of file 

   writeLines("SatScan-R Mapping Program - Report",con=TxtCon)
   writeLines(" ",con=TxtCon)
   writeLines(paste("   software version -> ",SSMVersion,sep=""),con=TxtCon)
   writeLines(" ",con=TxtCon)
   writeLines(date(),con=TxtCon)
   writeLines(" ",con=TxtCon)
   writeLines(" ",con=TxtCon)
   writeLines(DataTypeLit,con=TxtCon)
   writeLines(" ",con=TxtCon)
   writeLines(paste0("Census Year:",censusYear),con=TxtCon)
   writeLines(" ",con=TxtCon)
   writeLines(TypeForm,con=TxtCon)
   writeLines(" ",con=TxtCon)
   writeLines(YearList,con=TxtCon)
   writeLines(" ",con=TxtCon)
   writeLines(paste0("Shape type:",SType),con=TxtCon)
   writeLines(" ",con=TxtCon)
   
    
   writeLines(" ",con=TxtCon)
   writeLines("Data Source Files:",con=TxtCon)
   writeLines(" ",con=TxtCon)
   writeLines("  SatScan Results files:",con=TxtCon)
   writeLines(paste("     ",wColFileName,sep=""),con=TxtCon)
   writeLines(paste("     ",wGisFileName,sep=""),con=TxtCon)
   writeLines(paste("     ",wRRFileName,sep=""),con=TxtCon)
   
   #writeLines(" ",con=TxtCon)
   writeLines(" ",con=TxtCon)

   # If call TextRepOptions(TxtCon,LocBoundary,StateList) 

   writeLines(" ",con=TxtCon)
   writeLines("General Options:",con=TxtCon)
   writeLines(" ",con=TxtCon)
   writeLines(paste("   Run Identifier for Output Files:",runId,sep=""),con=TxtCon)
   writeLines(" ",con=TxtCon)
   writeLines("   Data references the following state(s) (Fips & Name): ",con=TxtCon)
   for (ind in 1:length(StateListData)) {
      xstr <- as.character(st99.index[StateListData[ind],"stName"])   # print list of State Names.
      writeLines(paste("      ",StateListData[ind],"  ",xstr,sep=""),con=TxtCon)
   }
   #writeLines(" ",con=TxtCon)
   writeLines(" ",con=TxtCon)
   
   writeLines(paste("Number of cluster in *.col.dbf file:", NumClusCol,"."),con=TxtCon)
   writeLines(paste("Number of cluster in *.gis.dbf file:", NumClusGis,"."),con=TxtCon)
   writeLines(" ",con=TxtCon)
   writeLines(paste("Number of locations referenced in the data:", length(LocListData),"."),con=TxtCon)
   writeLines(" ",con=TxtCon)
   writeLines(" ",con=TxtCon)

   # Cluster Mapping Options:   
   writeLines("Cluster Mapping Options:",con=TxtCon)
   writeLines(" ",con=TxtCon)
   writeLines(paste("   Number of Categories     :",categ,sep=""),con=TxtCon)
   writeLines(paste("   Cluster Map Title        :",vTitle,sep=""),con=TxtCon)
   writeLines(paste("   Census Tract Border Color:",bndyCol,sep=""),con=TxtCon)
   writeLines(paste("   Draw Cluster Outlines    :",outline,sep=""),con=TxtCon)
   writeLines(paste("   Label Cluster Outlines   :",label,sep=""),con=TxtCon)
   writeLines(paste("   Location Obs/Exp Ratio Map   :",locO_EMap,sep=""),con=TxtCon)
   writeLines(paste("   Cluster Obs/Exp Ratio Map    :",clusO_EMap,sep=""),con=TxtCon)
   
   writeLines(" ",con=TxtCon)
   writeLines(" ",con=TxtCon)
   
   writeLines("Obs/Exp Ratio Categories:",con=TxtCon)
   xx <- capture.output(print(GisTableCat,row.names=FALSE,digits=8,right=FALSE))
   writeLines(xx ,con=TxtCon)
   #writeLines(" ",con=TxtCon)
   writeLines(" ",con=TxtCon)
   
   writeLines(" ",con=TxtCon)
   
   # If call  TextRepClus(TxtCon, categ, ClusLocMap, ClusClusMap, ClusClusLab, GisTableCat)  
   
   #  Build US header variables
   USRes              <- NULL
   USRes$TStates      <- us99.index$states
   #USRes$THsas        <- us99.index$hsas
   USRes$TCounty      <- us99.index$county
   USRes$TTracts      <- us99.index$tracts
   # format
   USRes$LTStates     <- sprintf("%8i",USRes$TStates)
   #USRes$LTHsas       <- sprintf("%8i",USRes$THsas)
   USRes$LTCounty     <- sprintf("%8i",USRes$TCounty)
   USRes$LTTracts     <- sprintf("%8i",USRes$TTracts)
   
   #cat("USRes - typeof:",typeof(USRes),"  class:",class(USRes)," Z-4947 \n")
 
   ####
   #
   #   Loop through clusters and report statistics on Cluster Info, Location Info
   #
   ####
   #cat("Loop through cluster and build report stats Z-4954 .\n")
   
   ESet       <- EllipseSetBase
   ESet_L     <- dim(ESet)[1]
   
   if (EUpValue) {
   
      writeLines(paste0("Number of clusters found with P_VALUE < 0.05 = ",ESet_L)) 
      writeLines("Only clusters are reported when location (gis) records are available for the cluster.")
      writeLines(" ") 
   } else {
      writeLines("The Cluster Results P_VALUE data is missing. All cluster mapped.")
      writeLines(paste0("Number of clusters mapped = ",ESet_L)) 
      writeLines(" ") 
   }
   ####
   #
   #  Changes in St/Co, St/Co/Place, and Location reports.
   #
   #  To support the tiered report, all St/Co, St/Co/Place and Location data frames for 
   #  each cluster are collected and saved for one cycle.
   #
   #  As each tier level is requested and generated, the matching rows in the data frames
   #  will be pulled out and printed.
   #
   #  If the older reports are also requested, then the 
   #  report will generated using the same data.frames for the cluster.
   
   ####
   
   ####
   #
   #  Loop through the cluster data, print it and pull the location information
   #  again and print it.
   #
   ####
   
   TotalHdrUS        <- "  US --     TotSts  TotCnty  TotTrts"
   writeLines(TotalHdrUS,con=TxtCon)
   USRes$Lead   <- "         "
   #str(USRes)
 
   cat(unlist(USRes[c("Lead","LTStates","LTCounty","LTTracts")]),"\n",sep=" ",file=TxtCon)
   writeLines(" ",con=TxtCon)
   
   #
   #  Essentually - EllipseSetBase is the main list of Clusters and Years.
   #
   
   if ( ESet_L > 0 ) {
      # have clusters..
      
      for (ind in c(1:ESet_L))  {
         
         # loop through cluster one by one.
         # pull together - Cluster Data
         #cat("Text : Cluster Number : ",ind,"\n")
         
         # Pull off one cluster record (into three sub-records for printing.)
         
         ESWork1          <- NULL
         ESWork1          <- ESet[ind,c("StrDate","EndDate","NUMBER_LOC","LLR","P_VALUE","NYr","HL")]
         
         ESWork2          <- NULL
         ESWork2          <- ESet[ind,c("OBSERVED","EXPECTED","ODE","REL_RISK","Cat")]
    
         ESWork3          <- NULL
         if (CoordinatesType == 0) {
            ESWork3          <- ESet[ind,c("LOC_ID","E_MAJOR","E_MINOR","E_ANGLE","X","Y","t_X","t_Y")]
         } else {
            ESWork3          <- ESet[ind,c("LOC_ID","E_MAJOR","E_MINOR","E_ANGLE","LATITUDE","LONGITUDE","t_X","t_Y")]
         }
         # pull together - Location Data for Cluster
         ClusNum          <- ColTable$CLUSTER[ind]
          
          # Format numerical values - numeric -> character
         ESWork1$LLR      <- sprintf("%.3f",ESWork1$LLR)
         ESWork1$P_VALUE  <- sprintf("%.6f",ESWork1$P_VALUE)
        
         ESWork2$EXPECTED <- sprintf("%.3f",ESWork2$EXPECTED)
         ESWork2$ODE      <- sprintf("%.3f",ESWork2$ODE)
         ESWork2$REL_RISK <- sprintf("%.3f",ESWork2$REL_RISK)
         

         # Get list of locations for this cluster area
         GisList          <- GisTable[GisTable$CLUSTER == ClusNum,]  
         
         #cat("GisList colnames:",colnames(GisList)," Z-4682 \n")

         GLGOrder         <- seq(1:length(GisList$Cluster))  # no sort order.
        
         # set HIGH or LOW indicator for Cluster in GisList (convert from numeric to characters)
         if (ESWork1$HL < 0) {
            # Indicate Cluster is a LOW cluster and order locations Low to High
            ESWork1$HL    <- as.character("LOW ")
            GLGOrder      <- order(GisList$LocCat,decreasing=T) 
         } else {
            if (ESWork1$HL > 0) {
               # Indicate Cluster is a HIGH cluster and order locations High to Low.
               ESWork1$HL <- as.character("HIGH")
               GLGOrder   <- order(GisList$LocCat,decreasing=F)
            }
         } 
         
         #cat("GLGOrder:",GLGOrder,"\n")
          
         ####
         #
         # Start text report for this cluster
         #
         ####
          
         #  Cluster Statistics.
         #cat("Text : Cluster Statistics.\n")
          
         writeLines(paste("CLUSTER # ",ClusNum," Details:",sep=""),con=TxtCon)
         writeLines("",con=TxtCon)
       
         names(ESWork1) <- c("Start Date","End Date","Number of Locs","LLR","P_Value","Num_Years","High/Low")         
         names(ESWork2) <- c("Observed","Expected","Obs/Exp","Relative Risk","Category")         
         names(ESWork3) <- c("Centroid_ID","E_Major","E_Minor","E_Angle","Org_X","Org_Y","Tran_X","Tran_Y")         
      
         #xx <- capture.output(ESWork1)
         xx <- paste("        ",capture.output(print(ESWork1,row.names=FALSE)),sep="")
         writeLines(xx,con=TxtCon)
        
         #xx <- paste("       ",capture.output(ESWork2),sep="")
         xx <- paste("            ",capture.output(print(ESWork2,row.names=FALSE)),sep="")
         writeLines(xx,con=TxtCon)
        
         #xx <- paste("       ",capture.output(ESWork3),sep="")
         xx <- paste("            ",capture.output(print(ESWork3,row.names=FALSE)),sep="")
         writeLines(xx,con=TxtCon)
        
        
         #  Catagorize locations based on Obs_Exp ratio value
         #
         #  End of cluster header and data
         #
         #####
   
         #####
         #
         #  Now gather the data.frames for the location, state/county, and state/county/place census tract
         #   reports.
         #
         #####
                   
         ####
         #
         #  Get Location rows related to cluster
         #
         #  Add state, county, place names information from tr99.index (wk99)
         #    possible layouts:  state only,  state & county only,  state county place tract.
         #    written for s c t only, add extra set ups.
         #
         ####
         
         ####
         #
         #  Get and set up detail record - based on level
         #
         ####
         
         #cat("Get GisLocList Z-5118 - colnames:",colnames(GisList),"\n")
         
         #
         #  Modify if NORMAL supported.
         #
         
         # Get the list of locations belonging to this cluster ONLY  (GisLocList)
         # Reorder detail list.
         
         GisLocList <- GisList[GLGOrder,c("LOC_ID","CLUSTER","LOC_OBS","LOC_EXP","LOC_ODE","LOC_RR","LocCat")]
         lenGisLoc  <- length(GLGOrder)
         
         #cat("GisLocList - len:",lenGisLoc,"  typeof:",typeof(GisLocList),"  class:",class(GisLocList),"\n")
                  
         if (lenGisLoc > 0) { # there are locations in the cluster.
            #cat("Text:Gathering location records:",GisLocList$LOC_ID,"\n")
            
            GisLocList$LOC_ID <- as.character(GisLocList$LOC_ID)
            #  Bug - ID was used as factor not character.
              
            names(GisLocList) <- c("LOC_ID","Cluster","Observed","Expected","Obs_Exp","RelRisk","Category")
            
            # based on type of data - LocIDType = 2, 5, 11 (State, St/Co, St/Co/Tr)
            
            if (LocIDType == 2) {
               leadBlk            <- "    "
               #  add additional fields for State Location for aggregation
               
               GisLocList$stID    <- GisLocList$LOC_ID
               GisLocList$State   <- st99.index[GisLocList$stID,"stName"]
               #GisLocList$HsaID   <- ""
               #  Key # 1 - stcoID   ->  State/County code
               GisLocList$stcoID  <- ""
               GisLocList$County  <- ""
               GisLocList$plKey   <- ""
               GisLocList$Place   <- ""
               
               #  Counts
               GisLocList$TState   <- 1
               GisLocList$StateIn  <- 1
               #GisLocList$THsas    <- st99.index(GisLocList$stID,"hsas"]
               #GisLocList$THsasIn  <- st99.index(GisLocList$stID,"hsas"]
               GisLocList$TCounty  <- st99.index[GisLocList$stID,"county"]
               GisLocList$CountyIn <- st99.index[GisLocList$stID,"county"]
               GisLocList$TTracts  <- st99.index[GisLocList$stID,"tracts"]
               GisLocList$TractsIn <- st99.index[GisLocList$stID,"tracts"]
            }
            
            #if (LocIDType == 3) {   # State/HSAs
            #   leadBlk            <- "        "
            #   #  add additional fields for State/Hsas Location
            #   GisLocList$HsaID   <- str_trim(GisLocList$LOC_ID)
            #   GisLocList$stID    <- hs99.index[GisLocList$LOC_ID,"stID"] # go to HSA table and look up state.
            #   GisLocList$State   <- st99.index[GisLocList$stID,"stName"]
            #   GisLocList$stcoID  <- ""
            #   GisLocList$County  <- ""
            #   GisLocList$plKey   <- ""
            #   GisLocList$Place   <- ""
            #   
            #   # Counts
            #   GisLocList$TState   <- 0
	    #   GisLocList$StateIn  <- 0
	    #   GisLocList$THsas    <- 1
	    #   GisLocList$HsasIn   <- 1
	    #   GisLocList$TCounty  <- hs99.index[GisLocList$HsaID,"county"]
	    #   GisLocList$CountyIn <- hs99.index[GisLocList$HsaID,"county"]
	    #   GisLocList$TTracts  <- hs99.index[GisLocList$HsaID,"tracts"]
	    #   GisLocList$TractsIn <- hs99.index[GisLocList$HsaID,"tracts"]
            #   
            #}
            
            if (LocIDType == 5) {
               leadBlk            <- "        "
               #  add additional fields for State/County Location
               GisLocList$stID    <- str_sub(GisLocList$LOC_ID,1,2)
               GisLocList$State   <- st99.index[GisLocList$stID,"stName"]
               #  Key # 1 - stcoID   ->  State/County code
               GisLocList$stcoID  <- GisLocList$LOC_ID  # 5 digit fip code for state county
               GisLocList$County  <- co99.index[GisLocList$stcoID,"coName"]
               #GisLocList$HsaID   <- co99.index(GisLocList$stcoID,"HsaID"]
               GisLocList$plKey   <- ""
               GisLocList$Place   <- ""
               
               #  Counts
               GisLocList$TState   <- 0
               GisLocList$StateIn  <- 0
               #GisLocList$THsas    <- 0
               #GisLocList$HsasIn   <- 0
               GisLocList$TCounty  <- 1
               GisLocList$CountyIn <- 1
               GisLocList$TTracts  <- co99.index[GisLocList$stcoID,"tracts"]
               GisLocList$TractsIn <- co99.index[GisLocList$stcoID,"tracts"]
            
            }
            
            if (LocIDType == 11) {
               leadBlk            <- "            "
               #  add additional fields for State/County/Tract location
               GisLocList$stID    <- str_sub(GisLocList$LOC_ID,1,2)
               GisLocList$State   <- st99.index[GisLocList$stID,"stName"]
               #  Key # 1 - stcoID   ->  State/County code
               GisLocList$stcoID  <- substr(GisLocList$LOC_ID,1,5)  # 5 digit fip code for state county
               GisLocList$County  <- co99.index[GisLocList$stcoID,"coName"]
               #GisLocList$HsaID   <- co99.index(GisLocList$stcoID,"HsaID"]
              
               GisLocList$plKey   <- tr99.index[GisLocList$LOC_ID,"plKey"]
               GisLocList$Place   <- str_sub(GisLocList$plKey,6)
  
               #  Counts
               GisLocList$TState   <- 0
               GisLocList$StateIn  <- 0
               #GisLocList$THsas    <- 0
               #GisLocList$HsasIn   <- 0
               GisLocList$TCounty  <- 0
               GisLocList$CountyIn <- 0
               GisLocList$TTracts  <- 1
               GisLocList$TractsIn <- 1
            
            }
            
            #
            #  Design:
            #     US summary record
            #       Cluster summary record (CLUS)
            #          State Detailed (LOC)
            #
            #     US summary record
            #       Cluster summary record (CLUS)
            #          State Accum 
            #             County Detailed (LOC)
            #
            #  Alternate # 1
            #
            #     US summary record
            #       Cluster summary record (CLUS)
            #          State Accum 
            #             County Detailed (LOC)
            #
            #     US summary record
            #       Cluster summary record (CLUS)
            #          State Accum
            #             County Accum
            #                Place Accum
            #                   Tract Detailed (LOC)
            #
            #  Alternate # 2
            #
            #     US summary record
            #       Cluster summary record (CLUS)
            #          State Accum 
            #             Hsa Accum 
            #                County Detailed (LOC)
            #
            #     US summary record
            #       Cluster summary record (CLUS)
            #          State Accum
            #             Hsa Accum 
            #                County Accum
            #                   Place Accum
            #                      Tract Detailed (LOC)
            #
            #
            #
            #  Now we have the literal names.  We can calculate the widths of each element.
            #
          
            #  Get maximum width of each field - State, County, Place names
            MaxState      <- max(sapply(GisLocList$State,    function(x) nchar(as.character(x))))
            #MaxHsa        <- max(sapply(GisLocList$Hsa,      function(x) nchar(as.character(x))))
            MaxCounty     <- max(sapply(GisLocList$County,   function(x) nchar(as.character(x))))
            MaxPlace      <- max(sapply(GisLocList$Place,    function(x) nchar(as.character(x))))
            MaxCategory   <- max(sapply(GisLocList$Category, function(x) nchar(as.character(x)))) + 2  
            if (MaxCategory <= 14)  MaxCategory = 14   # minimum category width
              
            MaxColumn     <- 11 + 12 + 2   # standard column max.
            if (MaxColumn < (MaxState + 2))  MaxColumn <- MaxState + 2
            #if (MaxColumn < (MaxHsa + 2))    MaxColumn <- MaxHsa + 2
            if (MaxColumn < (MaxCounty + 5)) MaxColumn <- MaxCounty + 5
            if (MaxColumn < (MaxPlace + 10)) MaxColumn <- MaxPlace + 10
              
            # ssscccppplll = 12 characters
            # Longest  State + 2, Hsa + 2, County + 3, Place + 6, Location + 12 (11+12)  - all plus 2
            # 
            
            #
            #  Now finish creating any fields that are based on column width.
            #
            #  Create format strings using max string widths.  (space from left padding sections)
            
            SCPFormat       <- paste0("%-",MaxState,"s %-",MaxCounty,"s %-",MaxPlace,"s")    # State/County/Place
            SCFormat        <- paste0("%-",MaxState,"s %-",MaxCounty,"s")                    # State/County
            SFormat         <- paste0("%-",MaxState,"s")                                     # State

            #  Additional Fields based on column sizes.
            
            GisLocList$RIName        <- str_pad(GisLocList$LOC_ID,MaxColumn,"right")
            # FIPS code.
            GisLocList$FICol         <- str_pad(paste0("         ",GisLocList$LOC_ID),MaxColumn,"right")
            
             
            if (idMode == 1) {
               # setup State labels
               GisLocList$Trailer <- sprintf(SFormat,GisLocList$State) 
               GisLocList$Name    <- str_pad(paste0("  ",GisLocList$State),MaxColumn,"right")
            }
            if (idMode == 2) {
               # setup State/County labels
	       GisLocList$Trailer <- sprintf(SCFormat,GisLocList$State,GisLocList$County) 
	       GisLocList$Name    <- str_pad(paste0("      ",GisLocList$County),MaxColumn,"right")
            }
            if (idMode == 3) {
               # setup State/County/Place/Tract labels
	       GisLocList$Trailer <- sprintf(SCPFormat,GisLocList$State,GisLocList$County,GisLocList$Place) 
	       GisLocList$PName   <- str_pad(paste0("          ",GisLocList$Place),MaxColumn,"right")
               GisLocList$Name    <- str_pad(paste0("                ",GisLocList$LOC_ID),MaxColumn,"right")
            }
    
            #print(str(GisLocList))
            
            #
            #  State Detailed
            #  Hsa Detailed
            #  County Detailed
            #  Tract Detailed
            #
            #  US Accum
            #  State Accum
            #  County Accum
            #
            #  US Name
            #  State Name
            #  State/County Name
            #  State/County/Place Name
            #
            
            # Numbers sections:
            #    State Accum
            #    Hsa Accum
            #    County Accum
            #    Place Accum
            
            #    Titles
            
            #
            #  Update to include HSA counts - may need to have as option.
            #
            
            #   Basic Pattern for ODE, States, Counties, Tracts
            TObsExpODE        <- "  Observed  Expected  Obs/Exp"
            #                     12345678901234567890123456789
            #                         10        10        9
            
            TotalHdrUS        <- "   TotStates  TotCounty  TotTracts"
            
            TAccSt            <- "    TotSts  InClsSts PerStInCls"
            TAccCo            <- "   TotCnty InClsCnty PerCoInCls"
            TAccTr            <- "   TotTrts InClsTrts  PerTInCls"
            #                     1234567890123456789012345678901
            #                         10        10       11
           
            TTotSt            <- "   TotCnty   TotTrts"   # for state detail record
            #                     12345678901234567890
            #                         10        10
            
            catStr            <- str_pad("Category",MaxCategory,"left")
           
            DetailExtHdr      <- paste0("   RelRisk ",catStr,"  St/Co/Place")
            DetailExtHdr2     <- paste0("   RelRisk ",catStr)
            #                            1234567890
            #                                10
            
            #
            # US Summary header for cluster
            #
            AccumHdrUSAll     <- paste0(TObsExpODE,TAccSt,TAccCo,TAccTr)

            #
            # for State Data
            #
            # TotalHdrUS
            AccumHdrUS        <- paste0(TObsExpODE,TAccSt,TAccCo,TAccTr)        # collection of states in cluster
            DetailHdrSt       <- paste0(TObsExpODE,TTotSt,DetailExtHdr2)  # represents one state
            
            # for County Data
            AccumHdrUS        <- paste0(TObsExpODE,TAccSt,TAccCo,TAccTr)        # collection of states in cluster
            AccumHdrSt        <- paste0(TObsExpODE,TAccCo,TAccTr) 
            DetailHdrCo       <- paste0(TObsExpODE,TAccTr,DetailExtHdr2)  # represents one county
            
            # for Tract Data
            AccumHdrUS        <- paste0(TObsExpODE,TAccSt,TAccCo,TAccTr)        # collection of states in cluster
            AccumHdrSt        <- paste0(TObsExpODE,TAccCo,TAccTr) 
            AccumHdrCo        <- paste0(TObsExpODE,TAccTr)            
            AccumHdrPl        <- paste0(TObsExpODE,TAccTr)                # represents one place 
            DetailHdrTr       <- paste0(TObsExpODE,DetailExtHdr2)
            
            xBasicHeaderL      <- "   Observed   Expected    Obs/Exp  TotTracts  InClsTrts  PerInClus"
            xBasicHeaderS      <- "   Observed   Expected    Obs/Exp"
                      
            ####
            # 
            #    DETAILED RECORD numbers (same for each level).
            #
            ####
              
            #cat(" Z-4998 - GisLocList:\n")
            #print(str(GisLocList))
            # format literals - detail records should not change.
            
            # Calculate needed values and format  (All levels)
            GisLocList$LObserved    <- sprintf("%10i",GisLocList$Observed) 
            GisLocList$LExpected    <- sprintf("%10.3f",GisLocList$Expected)
            GisLocList[,"LObs_Exp"] <- sprintf("%9.3f",GisLocList[,"Obs_Exp"])
            
            GisLocList$LRelRisk     <- sprintf("%10.2f",GisLocList$RelRisk)
            GisLocList$LCategory    <- str_pad(as.character(GisLocList$Category),MaxCategory,"left")
            
            # Used at State Level
            
            GisLocList$LTCounty     <- sprintf("%10i", GisLocList$TCounty)
            GisLocList$LCountyIn    <- sprintf("%10i", GisLocList$CountyIn)
            xTF                     <- GisLocList$TCounty == 0
            #cat("check for zero T County:",sum(xTF),"\n")
            if (any(xTF)) {
               GisLocList[xTF,"TCounty"]  <- NA
               GisLocList[xTF,"CountyIn"] <- NA
            }
            GisLocList$PCCoIn       <- GisLocList$CountyIn/GisLocList$TCounty * 100
            GisLocList$LPCCoIn      <- sprintf("%10.2f%%", GisLocList$PCCoIn)
            
            # Used at State, County and Place Levels
            GisLocList$LTTracts     <- sprintf("%10i", GisLocList$TTracts)
            GisLocList$LTractsIn    <- sprintf("%10i", GisLocList$TractsIn)
            xTF                     <- GisLocList$TTracts == 0
            #cat("check for zero T tracts:",sum(xTF),"\n")
            if (any(xTF)) {
               GisLocList[xTF,"TTracts"]  <- NA
               GisLocList[xTF,"TractsIn"] <- NA
            }
            GisLocList$PCTrIn       <- GisLocList$TractsIn/GisLocList$TTracts * 100
            GisLocList$LPCTrIn      <- sprintf("%10.2f%%", GisLocList$PCTrIn)
  	    
            #GisLocHeader            <- paste0("   Census Tract  ",DetailHdrTr,DetailExtHdr)
            #GisLocHeaderT           <- paste0(str_pad("        Census Tract",MaxColumn,"right"),DetailHdrTr,DetailExtHdr2)
        
            #cat("GisLocList - after Z-4906:\n")
            #print(str(GisLocList))
            ####
            #
            #  Build report based on level of data
            #
            ####
            #cat("idMode:",idMode,"\n")
            
            if (idMode == 1)  {
            
               #  Design:
               #      Cluster
               #         US Level - Cluster Summary
               #            State Detail
               #
               ##  US level - Cluster Summary
               #USRes$Lead       <- "                   "
               #
               #USRes$StatesIn   <- length(unique(GisLocList$stID))
               #USRes$LStatesIn  <- sprintf("%10i",USRes$StatesIn)
               #USRes$PCStIn     <- USRes$StatesIn/USRes$TStates * 100
               #USRes$LPCStIn    <- sprintf("%10.2f%%",USRes$PCStIn)
               #
               #USRes$CountyIn   <- sum(GisLocList$NCounty)
               #USRes$LCountyIn  <- sprintf("%10i",USRes$CountyIn)
               #USRes$PCCoIn     <- USRes$CountyIn/USRes$TCounty * 100
               #USRes$LPCCoIn    <- sprintf("%10.2f%%",USRes$PCCoIn)
               #
               #USRes$TractsIn   <- sum(GisLocList$NTracts)
               #USRes$LTractsIn  <- sprintf("%10i",USRes$TractsIn)
               #USRes$PCTrIn     <- USRes$TractsIn/USRes$TTracts * 100
               #USRes$LPCTrIn    <- sprintf("%10.2f%%",USRes$PCTrIn)
               #
               #USRes$Observed   <- sum(GisLocList$Observed)
               #USRes$Expected   <- sum(GisLocList$Expected)
               #USRes$ODE        <- USRes$Observed/USRes$Expected
               #USRes$LObserved  <- sprintf("%10i",USRes$Observed)
               #USRes$LExpected  <- sprintf("%10.3f",USRes$Expected)
               #USRes$LODE       <- sprintf("%9.3f",USRes$ODE)
               #
               #writeLines(paste0(" US Cluster Summary:",AccumHdrUS),con=TxtCon)   
               #cat(unlist(USRes[c("Lead","LObserved","LExpected","LODE",
               #                      "LTStates","LStatesIn","LPCStIn",
               #                      "LTCounty","LCountyIn","LPCCoIn",
               #                      "LTTracts","LTractsIn","LPCTrIn")]),"\n",sep=" ",file=TxtCon)
            
               #
               #  State Detail
               #
               # order list in decending ODE order.
               
               GisLOrd <- order(GisLocList[,"Obs_Exp"],decreasing=TRUE)
               
               #  State Header
               writeLines(" ",con=TxtCon)
               writeLines(paste0(str_pad("  State",MaxColumn,"right"),DetailHdrSt),con=TxtCon)
               
               lenStList  <- dim(GisLocList)[1]
               #cat("length of state list:",lenStList,"\n")
               
               #  State Detail Records
               for (inSt in c(1:lenStList)) {
                  
                  inSSt   <- GisLOrd[inSt]
                  wStr    <- unlist(GisLocList[inSSt,c("Name","LObserved","LExpected","LObs_Exp",
                                   "LCountyIn","LTractsIn","LRelRisk","LCategory")],
                                   use.names=FALSE)
                  cat(wStr,"\n",sep="",append=TRUE, file=TxtCon)
      
               }
               writeLines(" ",con=TxtCon)
            }  # end of idMode = 1
            
            if (idMode == 2 ) {
               #  State/County data
               #
               #  Design:
               #     Cluster 
               #        US cluster Summary (sum)
               #           State Summary (agg)
               #              County Detail...
               #           State Summary (agg)
               #              County Detail...
               #
               #  US level - Cluster Summary
               #USRes$Lead       <- "                   "
               #
               #USRes$StatesIn   <- length(unique(GisLocList$stID))
               #USRes$LStatesIn  <- sprintf("%10i",USRes$StatesIn)
               #USRes$PCStIn     <- USRes$StatesIn/USRes$TStates * 100
               #USRes$LPCStIn    <- sprintf("%10.2f%%",USRes$PCStIn)
               #
               #USRes$CountyIn   <- sum(GisLocList$CountyIn)
               #USRes$LCountyIn  <- sprintf("%10i",USRes$CountyIn)
               #USRes$PCCoIn     <- USRes$CountyIn/USRes$TCounty * 100
               #USRes$LPCCoIn    <- sprintf("%10.2f%%",USRes$PCCoIn)
               #
               #USRes$TractsIn   <- sum(GisLocList$TractsIn)
               #USRes$LTractsIn  <- sprintf("%10i",USRes$TractsIn)
               #USRes$PCTrIn     <- USRes$TractsIn/USRes$TTracts * 100
               #USRes$LPCTrIn    <- sprintf("%10.2f%%",USRes$PCTrIn)
               #
               #USRes$Observed   <- sum(GisLocList$Observed)
               #USRes$Expected   <- sum(GisLocList$Expected)
               #USRes$ODE        <- USRes$Observed/USRes$Expected
               #USRes$LObserved  <- sprintf("%11i",USRes$Observed)
               #USRes$LExpected  <- sprintf("%10.3f",USRes$Expected)
               #USRes$LODE       <- sprintf("%9.3f",USRes$ODE)
               #
               #writeLines(paste0(" US Cluster Summary:",AccumHdrUS),con=TxtCon)   
               #cat(unlist(USRes[c("Lead","LObserved","LExpected","LODE",
               #                      "LTStates","LStatesIn","LPCStIn",
               #                      "LTCounty","LCountyIn","LPCCoIn",
               #                      "LTTracts","LTractsIn","LPCTrIn")]),"\n",sep=" ",file=TxtCon)
            
               ####
               #
               #  States - step through list
               #
               ####
               #
               # State List - and aggregate needed fields
               #   State Accum
               #
               ####
           
               stxx1 <- aggregate(Observed ~ State + stID, data=GisLocList, FUN=length)   # get number of counties
               stxx2 <- aggregate(cbind(Observed, Expected, CountyIn, TractsIn) ~ State + stID , data=GisLocList, FUN=sum)  # get sum of Obs and Exp
                 
               # the assumption is the rows in xx1 and xx2 are in the same order and match - one to one after
               #  the two aggregations.
               
               #  State Level Summary 
               ResSt           <- NULL
               #cat("dim(stxx1):",dim(stxx1),"\n")
               ResSt           <- stxx1                  # by state - number of counties involved
               
               colnames(ResSt) <- c("State","stID","CountyIn")
               ResSt$stID      <- as.character(ResSt$stID)
               ResSt$State     <- as.character(st99.index[ResSt$stID,"stName"])
               
               ResSt$Name      <- str_pad(paste0("  ",ResSt$State),MaxColumn,"right")
               #ResSt$FCol      <- str_pad(ResSt$State,MaxColumn,"right")
          
               ResSt$Observed  <- stxx2$Observed
               ResSt$LObserved <- sprintf("%10.3f",ResSt$Observed)
               ResSt$Expected  <- stxx2$Expected
               ResSt$LExpected <- sprintf("%10.3f",ResSt$Expected)
               ResSt$Obs_Exp   <- ResSt$Observed / ResSt$Expected
               ResSt[is.nan(ResSt$Obs_Exp),"Obs_Exp"] = 0
               ResSt$LObs_Exp  <- sprintf("%9.3f",ResSt$Obs_Exp)
                
               ResSt$TCounty   <- st99.index[ResSt$stID,"county"]
               ResSt$LTCounty  <- sprintf("%10i",ResSt$TCounty)
               ResSt$CountyIn  <- stxx2$CountyIn
               ResSt$LCountyIn <- sprintf("%10i",ResSt$CountyIn)
               
               ResSt$PCCoIn    <- ResSt$CountyIn/ResSt$TCounty * 100
               ResSt$LPCCoIn   <- sprintf("%10.2f%%",ResSt$PCCoIn)
               
               ResSt$TTracts   <- st99.index[ResSt$stID,"tracts"]
               ResSt$LTTracts  <- sprintf("%10i",ResSt$TTracts)
               ResSt$TractsIn  <- stxx2$TractsIn
               ResSt$LTractsIn <- sprintf("%10i",ResSt$TractsIn)
               ResSt$PCTrIn    <- ResSt$TractsIn/ResSt$TTracts * 100
               ResSt$LPCTrIn   <- sprintf("%10.2f%%",ResSt$PCTrIn)
             
               #  State Header with County Detail
               
               ResSOrd <- order(ResSt$Obs_Exp,decreasing=TRUE)
               
               lenStList   <- length(ResSOrd)  # get number of states from aggregate
               #cat("Length of state Aggreg List:",lenStList,"\n")
               
               for (inSt in seq_len(lenStList)) {
              
                  #  state title header
                  writeLines(" ",con=TxtCon)
                  writeLines(paste0(str_pad("State",MaxColumn,"right"),AccumHdrSt),con=TxtCon)
              
                  inSSt  <- ResSOrd[inSt]      # get order of states within cluster
                  
                  # State Summary Aggregate 
                  wStr <- unlist(ResSt[inSSt,c("Name","LObserved","LExpected","LObs_Exp",
                                "LTCounty","LCountyIn","LPCCoIn",
                                "LTTracts","LTractsIn","LPCTrIn")],use.names=FALSE)
                               
                  cat(wStr,"\n",sep="",append=TRUE, file=TxtCon)
                  
                  #
                  #  Now do the detail county records for this state entry.
                  #
                  
                  # Written State Summary - now for the county detail
	          
                  CurStID          <- ResSt[inSSt,"stID"]
                  
                  # get part of list for one states counties.
                  CoLocList        <- GisLocList[GisLocList$stID == CurStID,]
                  
                  # get sorted order (high ODE to low)
                  CoLocOrd         <- order(CoLocList[,"Obs_Exp"],decreasing=TRUE)
                  
                  # nubmer of counties for this state within cluster 
                  lenCoList        <- dim(CoLocList)[1]
                  #cat("Length of County details list:",lenCoList,"\n")
                  
                  # County Details Header
                  writeLines(" ",con=TxtCon)
                  writeLines(paste0(str_pad("    County",MaxColumn,"right"),DetailHdrCo),con=TxtCon)
   
                  # County detail records
                  for (inCo in c(1:lenCoList)) {

                     inSCo    <- CoLocOrd[inCo]
                     cat(unlist(CoLocList[inSCo,c("Name","LObserved","LExpected","LObs_Exp",
                                           "LTTracts","LTractsIn","LPCTrIn",
                                           "LRelRisk","LCategory")]),"\n",
                                           sep="",file=TxtCon,append=TRUE)
                   
                  }  # end of county loop within state within cluster
                  
              } # end of state loop within cluster
           } # end of county cluster
           
           if (idMode == 3) {
              # Census tract.
              
              #  Design:
              #    US Statistics
              #       Cluster
              #          US Summary
              #             State Summary
              #                County Summary
              #                   Place Summary
              #                      Trace Detail
              #     ...
              #
              
              #  US level - Cluster Summary
              #USRes$Lead       <- "                   "
              #
              #USRes$StatesIn   <- length(unique(GisLocList$stID))
              #USRes$LStatesIn  <- sprintf("%10i",USRes$StatesIn)
              #USRes$PCStIn     <- USRes$StatesIn/USRes$TStates * 100
              #USRes$LPCStIn    <- sprintf("%10.2f%%",USRes$PCStIn)
              #
              #USRes$CountyIn   <- sum(GisLocList$NCounty)
              #USRes$LCountyIn  <- sprintf("%10i",USRes$CountyIn)
              #USRes$PCCoIn     <- USRes$CountyIn/USRes$TCounty * 100
              #USRes$LPCCoIn    <- sprintf("%11.2f%%",USRes$PCCoIn)
              #
              #USRes$TractsIn   <- sum(GisLocList$NTracts)
              #USRes$LTractsIn  <- sprintf("%10i",USRes$TractsIn)
              #USRes$PCTrIn     <- USRes$TractsIn/USRes$TTracts * 100
              #USRes$LPCTrIn    <- sprintf("%10.2f%%",USRes$PCTrIn)
              #
              #USRes$Observed   <- sum(GisLocList$Observed)
              #USRes$Expected   <- sum(GisLocList$Expected)
              #USRes$ODE        <- USRes$Observed/USRes$Expected
              #USRes$LObserved  <- sprintf("%10i",USRes$Observed)
              #USRes$LExpected  <- sprintf("%10.3f",USRes$Expected)
              #USRes$LODE       <- sprintf("%9.3f",USRes$ODE)
              #
              ##writeLines(paste0(" US Cluster Summary:",AccumHdrUS),con=TxtCon)   
              ##cat(unlist(USRes[c("Lead","LObserved","LExpected","LODE",
              ##                      "LTStates","LStatesIn","LPCStIn",
              ##                      "LTCounty","LCountyIn","LPCCoIn",
              ##                      "LTTracts","LTractsIn","LPCTrIn")]),"\n",sep=" ",file=TxtCon)
            
              ####
              #
              #  States - step through list of tracts...
              #
              ####
              #cat("GisLocList for cluster:\n")
              #print(GisLocList)
              
              #
              #  Tract to Place.  Aggregate all tracts to places for state.
              #
              #cat("Tract GisLocList to Place records\n")
	      plxx2 <- aggregate(cbind(Observed, Expected,  TractsIn) ~ State + stID + County + plKey + Place, data=GisLocList, FUN=sum)  # get sum of Obs and Exp
	      
	      #cat("plxx2 - aggregate SUM.\n")
	      #print(plxx2)
	      #print(str(plxx2))
	      #
	      #cat("Build Place table for tracts\n")
	      
	      #  State/County/Place Level Summary - Aggregate place to county.
	      
	      ResSCP           <- plxx2                  # number of tracts aggre up to Place.
	      colnames(ResSCP) <- c("State","stID","County","plKey","Place","Observed","Expected","TractsIn")
	      #ResSCP$stID     <- as.character(ResSCP$stID)
	      #ResSCP$State    <- as.character(st99.index[ResSCP$stID,"stName"])
	      ResSCP$stcoID    <- str_sub(ResSCP$plKey,1,5)
	      ResSCP$Name      <- str_pad(paste0("            ",ResSCP$Place),MaxColumn,"right")
	      #ResSCP$FCol     <- str_pad(ResSCP$Place,MaxColumn,"right")
	      
	      # sum tract Obs, Exp and ODE values within place.
	      #ResSCP$Observed  <- plxx2$Observed                      # sum tract values

	      ResSCP$LObserved <- sprintf("%10.3f",ResSCP$Observed)
	      #ResSCP$Expected  <- plxx2$Expected
	      ResSCP$LExpected <- sprintf("%10.3f",ResSCP$Expected)
	      ResSCP$Obs_Exp   <- ResSCP$Observed / ResSCP$Expected
	      ResSCP[is.nan(ResSCP$Obs_Exp),"Obs_Exp"] = 0
	      ResSCP$LObs_Exp  <- sprintf("%9.3f",ResSCP$Obs_Exp)
	      
	      # sum of county values within place - NONE.
	      
	      # sum tract values
	      ResSCP$TTracts   <- pl99.index[ResSCP$plKey,"tracts"]   # total tracts for state
	      ResSCP$LTTracts  <- sprintf("%10i",ResSCP$TTracts)
	      #ResSCP$TractsIn  <- plxx2$TractsIn                      # in cluster tracts for state
	      ResSCP$LTractsIn <- sprintf("%10i",ResSCP$TractsIn)
	      ResSCP$PCTrIn    <- ResSCP$TractsIn/ResSCP$TTracts * 100
	      ResSCP$LPCTrIn   <- sprintf("%10.2f%%",ResSCP$PCTrIn)
	      
	      #cat("ResSCP - Z-5211 :\n")
	      #print(str(ResSCP))
	      #print(head(ResSCP,30))
	      
	      #
	      #   Place to County
	      #
	      
              #cat("Place aggregate to County \n")
              
              #  County Summary  (Aggregates by County of Place from Tract)
              coxx2 <- aggregate(cbind(Observed, Expected, TractsIn) ~ State + stID + County + stcoID, data=ResSCP, FUN=sum)  # get sum of Obs and Exp
              
              #cat("coxx2 - Z-5364 :\n")
              #print(coxx2)
              
              #cat("Build County Aggre Table for tract\n")
              
              #  State/County Level Summary 
              ResSC           <- coxx2                # number of counties
              colnames(ResSC) <- c("State","stID","County","stcoID","Observed","Expected", "TractsIn")
              #ResSC$stID      <- as.character(ResSC$stID)                       #>>>
              #ResSC$State     <- as.character(st99.index[ResSC$stID,"stName"])  #???
              ResSC$Name      <- str_pad(paste0("        ",ResSC$County),MaxColumn,"right")
              #ResSC$FCol      <- str_pad(ResSC$County,MaxColumn,"right")
              
              #cat("ResSC:\n")
              #print(ResSC)
              #print(str(ResSC))
          
              # sum tract values
              #ResSC$Observed  <- coxx2$Observed                      # ??? sum place-tract values
              ResSC$LObserved <- sprintf("%10.3f",ResSC$Observed)
              #ResSC$Expected  <- coxx2$Expected                      # ???
              ResSC$LExpected <- sprintf("%10.3f",ResSC$Expected)
              ResSC$Obs_Exp   <- ResSC$Observed / ResSC$Expected
              ResSC[is.nan(ResSC$Obs_Exp),"Obs_Exp"] = 0
              ResSC$LObs_Exp  <- sprintf("%9.3f",ResSC$Obs_Exp)
              
              # Now at county level - set TCounty and CountyIn = 1
              ResSC$TCounty   <- 1                                   # total counties for state
              ResSC$LTCounty  <- sprintf("%10i",ResSC$TCounty)
              ResSC$CountyIn  <- 1                                   # in cluster counties for state
              ResSC$LCountyIn <- sprintf("%10i",ResSC$CountyIn)
              
              ResSC$PCCoIn    <- ResSC$CountyIn/ResSC$TCounty * 100
              ResSC$LPCCoIn   <- sprintf("%10.2f%%",ResSC$PCCoIn)
              
              # sum tract values
              ResSC$TTracts   <- co99.index[ResSC$stcoID,"tracts"]   # total tracts for state  (should sum, but this feels better)
              ResSC$LTTracts  <- sprintf("%10i",ResSC$TTracts)
              
              #ResSC$TractsIn  <- coxx2$TractsIn                      # ??? in cluster tracts for state
 
              ResSC$LTractsIn <- sprintf("%10i",ResSC$TractsIn)
              ResSC$PCTrIn    <- ResSC$TractsIn/ResSC$TTracts * 100
              ResSC$LPCTrIn   <- sprintf("%10.2f%%",ResSC$PCTrIn)
              
              #cat("ResSC aggr from ResSCP\n")
              #print(str(ResSC))
              #print(head(ResSC,20))
              
              #cat("ResSC$Obs_Exp:",ResSC$Obs_Exp,"\n")
              #cat(typeof(ResSC$Obs_Exp),"\n")
              
              #  
              #  County aggregates to State level
              #
	      
              #cat(" County aggregate To State (from tract->Place-> County) \n")
           
              stxx2 <- aggregate(cbind(Observed, Expected, CountyIn, TractsIn) ~ State + stID , data=ResSC, FUN=sum)  # get sum of Obs and Exp
                
              # the assumption is the rows in xx1 and xx2 are in the same order and match - one to one after
              #  the two aggregations.
              #cat("stxx2:\n")
              #print(stxx2)
              
              #cat("st99.index\n")
              #print(st99.index[,c("ID", "county_00", "county_10", "tracts_00", "tracts_10")])
              
              #cat("build state aggre table for tract.","\n")
              
              #  State Level Summary (Aggregates of County to State)
              ResSt           <- NULL
              #cat("dim(stxx2):",dim(stxx2),"\n")
              
              ResSt           <- stxx2
              colnames(ResSt) <- c("State", "stID", "Observed", "Expected", "CountyIn", "TractsIn")
              #ResSt$stID      <- as.character(ResSt$stID)    # ???
              #ResSt$State     <- as.character(st99.index[ResS$stID,"stName"])  # >>>
              ResSt$Name      <- str_pad(paste0("    ",ResSt$State),MaxColumn,"right")
              #ResSt$FCol      <- str_pad(ResSt$State,MaxColumn,"right")
          
              #cat("ResSt:\n")
              #print(str(ResSt))
              #print(ResSt)
              
              # sum tract values
              #ResSt$Observed  <- stxx2$Observed                      # sum  tract values
              ResSt$LObserved <- sprintf("%10.3f",ResSt$Observed)
              #ResSt$Expected  <- stxx2$Expected
              ResSt$LExpected <- sprintf("%10.3f",ResSt$Expected)
              ResSt$Obs_Exp   <- ResSt$Observed / ResSt$Expected
              ResSt[is.nan(ResSt$Obs_Exp),"Obs_Exp"] = 0
              ResSt$LObs_Exp  <- sprintf("%9.3f",ResSt$Obs_Exp)
              
              # New state values
              ResSt$TState    <- 1
              ResSt$StatesIn  <- 1
              
              # sum county values
              ResSt$TCounty   <- st99.index[ResSt$stID,"county"]   # total counties for state
              ResSt$LTCounty  <- sprintf("%10i",ResSt$TCounty)
              #ResSt$CountyIn  <- stxx2$CountyIn                      # in cluster counties for state
              ResSt$LCountyIn <- sprintf("%10i",ResSt$CountyIn)
              
              ResSt$PCCoIn    <- ResSt$CountyIn/ResSt$TCounty * 100
              ResSt$LPCCoIn   <- sprintf("%10.2f%%",ResSt$PCCoIn)
              
              # sum tract values
              ResSt$TTracts   <- st99.index[ResSt$stID,"tracts"]   # total tracts for state
              ResSt$LTTracts  <- sprintf("%10i",ResSt$TTracts)
              #ResSt$TractsIn  <- stxx2$TractsIn                      # in cluster tracts for state
              ResSt$LTractsIn <- sprintf("%10i",ResSt$TractsIn)
              ResSt$PCTrIn    <- ResSt$TractsIn/ResSt$TTracts * 100
              ResSt$LPCTrIn   <- sprintf("%10.2f%%",ResSt$PCTrIn)
              
              #  State Header with County
              
              #cat("ResS Z-5913 :\n")
              #print(str(ResS))
              #print(head(ResS,10))
              
              ResSOrd <- order(ResSt$Obs_Exp,decreasing=TRUE)
              #cat("ResSOrd:",ResSOrd,"\n")
              
              #cat("starting Tract Report\n")
              # report 
              
              lenStList   <- length(ResSOrd)
              
              # state Loop
              for (inSt in seq_len(lenStList)) {
                 writeLines(" ",con=TxtCon)
                 # state header
                 writeLines(paste0(str_pad("  State",MaxColumn,"right"),AccumHdrSt),con=TxtCon)
                 
                 inSSt        <- ResSOrd[inSt]
                 CurStID      <- ResSt[inSSt,"stID"] # Get state ID for this loop.
                 
                 # State Summary Aggregate
                 wStr <- unlist(ResSt[inSSt,c("Name","LObserved","LExpected","LObs_Exp",
                               "LTCounty","LCountyIn","LPCCoIn",
                               "LTTracts","LTractsIn","LPCTrIn")],use.names=FALSE)
                              
                 cat(wStr,"\n",sep="",append=TRUE, file=TxtCon)
                   
	         #  Aggregate to County
                 
                 #  Select counties within current state
                 
                 CoWrkRes     <- ResSC[ResSC$stID == CurStID,]    # get aggre for all counties in state.
                 CoWrkOrd     <- order(CoWrkRes$Obs_Exp,decreasing=TRUE)
                 lenCoWrkRes  <- length(CoWrkOrd)      # number of counties for state in cluster.

                 #cat("County agg within State:",CurStID,"\n")
                 #cat("Length of County Work:",lenCoWrkRes,"\n")
                 #cat("CoWrkRes:\n")
                 #print(str(CoWrkRes))
                 #print(head(CoWrkRes,20))
                 
                 #cat("CoWrkOrd:",CoWrkOrd,"\n")
                  
                 # county loop
                 for (inCo in seq_len(lenCoWrkRes)) {
                    # county header 
                    writeLines(" ",con=TxtCon)
                    writeLines(paste0(str_pad("      County",MaxColumn,"right"),AccumHdrCo),con=TxtCon)
                 
                    inSCo       <- CoWrkOrd[inCo]
                    CurStcoID   <- ResSC[inSCo,"stcoID"]  # get current county id
                    
                    # County Aggre Record
                    wStr      <- unlist(CoWrkRes[inSCo,c("Name","LObserved","LExpected","LObs_Exp",
		                            "LTTracts","LTractsIn","LPCTrIn")],use.names=FALSE)
		                            
		    cat(wStr,"\n",sep="",file=TxtCon,append=TRUE)
		    		    
                    # Aggregate to Place
                    #  get list of place summary records for state and county
                    plWrkRes     <- ResSCP[ResSCP$stID == CurStID & ResSCP$stcoID == CurStcoID,]
                    
                    plWrkOrd     <- order(plWrkRes$Obs_Exp,decreasing=TRUE)
                    
                    lenPlWrkRes  <- length(plWrkOrd)       # number of entries

                    #cat("Place agg within State/County:",CurStcoID,"\n")
                    #cat("Length of Place Work:",lenPlWrkRes,"\n")
                    #cat("plWrkRes:\n")
                    #print(str(plWrkRes))
                    #print(head(plWrkRes,20))
                    
                    #cat("plWrkOrd:",plWrkOrd,"\n")
                    
                    # Place Loop
                    for (inPl in seq_len(lenPlWrkRes)) {
                       writeLines(" ",con=TxtCon)
                       writeLines(paste0(str_pad("          Placename",MaxColumn,"right"),AccumHdrPl),con=TxtCon)
                 
                       inSPl     <- plWrkOrd[inPl]    # ordered index to placename data line
                       
                       CurPlKey  <- plWrkRes[inSPl,"plKey"]   # current entry
                       #cat("CurPlKey:",CurPlKey,"\n")
                       
                       # Place Aggre Record - header.
                       wStr     <- unlist(plWrkRes[inSPl,c("Name","LObserved","LExpected","LObs_Exp",
                                          "LTTracts","LTractsIn","LPCTrIn")],use.names=FALSE)
                       
                       cat(wStr,"\n",sep="",file=TxtCon,append=TRUE)
                 
                       # Tract detail within place
                       
                       #cat("Tract Detail for place :",CurPlKey,"\n")
                       
                       trLocList  <- GisLocList[GisLocList$plKey == CurPlKey,]
                       #cat("trLocList DF:\n")
                       #print(str(trLocList))
                       #print(trLocList)
                       
                       trLocOrd   <- order(trLocList$Obs_Exp,decreasing=TRUE)
                       lenTrLoc   <- length(trLocOrd)
                     
                       #cat("Tract Detail within State/County/Place:",CurPlKey,"\n")
                       #cat("Length of Tract Work:",lenTrLoc,"\n")
                       #cat("trLocList:\n")
                       #print(str(trLocList))
                       #print(head(trLocList,20))
      
                       #cat("trLocOrd:",trLocOrd,"\n")
                       
                       #writeLines(" ",con=TxtCon)
                       writeLines(paste0(str_pad("              Census Tract",MaxColumn,"right"),DetailHdrTr),con=TxtCon)
                       
                       # Tract Loop
                       for (inTr in seq_len(lenTrLoc)) {
                          inSTr  <- trLocOrd[inTr]
                          wStr   <- unlist(trLocList[inSTr,c("Name","LObserved","LExpected","LObs_Exp",
                                                 "LRelRisk","LCategory")],use.names=FALSE)
                          cat(wStr,"\n",sep="",file=TxtCon,append=TRUE)
                       }  # end of trace detail loop with in place 
                    
                    }  # end of place loop with in county
	         
	         }  # end of county loop within state within cluster
              } # end of state loop within cluster
           } # end of tract cluster  (idMode = 3)
  
     #  #  State Header
     #          writeLines(paste0(str_pad("State",MaxState+2,"right"),DetailHdrSt),con=TxtCon)
     #          
     #          lenStList  <- dim(GisLocList)[1]
     #          #  State Detail Records
     #          for (inSt in c(1:lenStList)) {
     #             
     #             inSSt <- GisLOrd[inSt]
     #             wStr <- unlist(GisLocList[inSSt,c("LName","LObserved","LExpected","LObs_Exp",
     #                                 "LNCounty","LNTracts","LRelRisk","LCategory")],
     #                                 use.names=FALSE)
     #             cat(" ",wStr,"\n",sep=" ",append=TRUE, file=TxtCon)
     # 
     #          }
     #          writeLines(" ",con=TxtCon)
     #          
     #         
     #         }
     #          
     #       }
     #         
     #         
     #          SHeader        <- paste0(str_pad("State",MaxState+2,"right"),BasicHeaderL)
     #          SHeaderT       <- paste0(str_pad("State",MaxColumn,"right"),BasicHeaderL)
     #            
     #          
     #            
     #         # Build US Accum record and output
      #         USResSt$NStates   <- sprintf("%3i",   xx1)
      #         USResSt$NCounties <- sprintf("%5i",   xx2$county) 
      #         USResSt$NTracts   <- sprintf("%7i",   xx2$tract)
      #         USResSt$Observed  <- sprintf("%10i",  xx2$Observed)
      #         USResSt$Expected  <- sprintf("%10.3f",xx2$Expected)
      #         USResSt$Obs_Exp   <- sprintf("%10.3f",xx2$Observed / xx2$Expected)
      #         USResSt[is.nan(USResSt$Obs_Exp),"Obs_Exp"] = 0
      #         
      #         # report of State Detail - LOC records
      #         
      #         GisLocHeader      <- paste0("        States      ",DetailHdrTr,DetailExtHdr)
      #         GisLocHeaderT     <- paste0(str_pad("        States      ",MaxColumn,"right"),DetailHdrTr,DetailExtHdr2)
#
      #         # Tier 1 - State Header 
      #          
      #         writeLines(SHeaderT,con=TxtCon)
      #        
      #         for (iS in 1:NumS) {   # Loop at Tier 1
      #            # Tier 1 - State Data
      #            #cat("Text:Rep-State.\n")
      #            
      #            cat(unlist(ResSt[OrdS[iS],c("FCol","Observed","Expected","Obs_Exp","TotTracts","Tracts_In","PrcTracts")]),"\n",sep=" ",file=TxtCon,append=TRUE)
      #      
      #     
      #         }
   #
   #
   #
      #      } 
      #         
      #      
      ##      #
       #     # State List - and aggregate needed fields
       #     #   State Accum
       #     #
       #     ####
       #    
       ##     xx2 <- aggregate(cbind(Observed, Expected) ~ State + stID , data=GisLocList, FUN=sum)  # get sum of Obs and Exp
        #         
        #    # the assumption is the rows in xx1 and xx2 are in the same order and match - one to one after
        #    #  the two aggregations.
        #    xLen           <- MaxState+2
        #      
        #    ResS           <- xx1   # number of tracts
        #    colnames(ResS) <- c("State","stID","Tracts_In")
        #    ResSt$stID      <- as.character(ResSt$stID)
        #    ResSt$State     <- as.character(ResSt$State)
        #    ResSt$RName     <- str_pad(ResSt$State,xLen,"right")
        #    ResSt$FCol      <- str_pad(ResSt$State,MaxColumn,"right")
        # 
        #    ResSt$Observed  <- sprintf("%10i",xx2$Observed)
        #    ResSt$Expected  <- sprintf("%10.3f",xx2$Expected)
        #    ResSt$Obs_Exp   <- sprintf("%10.3f",xx2$Observed / xx2$Expected)
        #    ResSt[is.nan(ResSt$Obs_Exp),"Obs_Exp"] = 0
        #    
        #    ResSt$TotTracts <- st99.index[ResSt$stID,TractN]
        #    
        #    ResSt$PrcTracts <- sprintf("%10.3f%%",ResSt$Tracts_In/ResSt$TotTracts * 100)
        #    ResSt$Tracts_In <- sprintf("%10i",ResSt$Tracts_In)
        #    
        #    ResSt$TotTracts <- sprintf("%10i",ResSt$TotTracts)
        #      
        #    SHeader        <- paste0(str_pad("State",MaxState+2,"right"),BasicHeaderL)
        #    SHeaderT       <- paste0(str_pad("State",MaxColumn,"right"),BasicHeaderL)
        #    
        #    
        #    if (LocIDType >= 5) {
        #       ####
        #       #
        #       # State/County List - and aggregate needed fields
        #       #    County Accum
        #       #
        #       ####
        #      
        #       xx1 <- aggregate(Observed ~  St_County + State + County + stcoID, data=GisLocList, FUN=length)   # get number of tracts
        #       xx2 <- aggregate(cbind(Observed, Expected) ~ St_County + State + County + stcoID , data=GisLocList, FUN=sum)  # get sum of Obs and Exp
        #         
        #       # the assumption is the rows in xx1 and xx2 are in the same order and match - one to one after
        #       #  the two aggregations.
        #       xLen            <- MaxState+MaxCounty+2
        #      
        #       ResSC           <- xx1   # number of tracts
        #       colnames(ResSC) <- c( "State_County","State","County","stcoID","Tracts_In")
        #       ResSC$State     <- as.character(ResSC$State)
        #       ResSC$State_County <- as.character(ResSC$State_County)
        #       ResSC$County    <- as.character(ResSC$County)
        #     
        #       ResSC$RName     <- str_pad(ResSC$State_County,xLen,"right")
        #       ResSC$FCol      <- str_pad(paste0("   ",ResSC$County,sep=""),MaxColumn,"right")
        #  
        #       ResSC$Observed  <- sprintf("%10i",xx2$Observed)
        #       ResSC$Expected  <- sprintf("%10.3f",xx2$Expected)
        #       ResSC$Obs_Exp   <- sprintf("%10.3f",xx2$Observed / xx2$Expected)
        #       ResSC[is.nan(ResSC$Obs_Exp),"Obs_Exp"] = 0
        #     
        #       ResSC$TotTracts <- co99.index[ResSC$stcoID,TractN]
        #       ResSC$PrcTracts <- sprintf("%10.3f%%",ResSC$Tracts_In/ResSC$TotTracts * 100)
        #       ResSC$Tracts_In <- sprintf("%10i",ResSC$Tracts_In)
        #       ResSC$TotTracts <- sprintf("%10i",ResSC$TotTracts)
        #       
        #       SCHeader        <- paste0(str_pad("State/County",xLen,"right"),BasicHeaderL)
        #       SCHeaderT       <- paste0(str_pad("   County",MaxColumn,"right"),BasicHeaderL)
        #    
        #       if (LocIDType >= 11) {
        #    
        #          ####
        #          #
        #          # State/County/Place List - and aggregate needed fields
        #          #   Place Accum
        #          #
        #          ####
        #           
        #          xx1 <- aggregate(Observed ~   St_Co_Place + St_County + Place + plKey, data=GisLocList, FUN=length)               # count number of locations
        #          xx2 <- aggregate(cbind(Observed, Expected) ~ St_Co_Place + St_County + Place + plKey, data=GisLocList, FUN=sum)  # sum Observed and Expected values
        #            
        ##          #  the two aggregations.
         #           
         #         # problem with below code (???).
         #         xLen               <- MaxState+MaxCounty+MaxPlace+2
         #         
         #         ResSCP             <- xx1                                              # number of locations
         #         colnames(ResSCP)   <- c("State_County_Place", "State_County", "Place", "plKey","Tracts_In")
         #         ResSCP$State_County_Place <- as.character(ResSCP$State_County_Place)
         #         ResSCP$State_County <- as.character(ResSCP$State_County)
         #         ResSCP$Place       <- as.character(ResSCP$Place)
         #         ResSCP$RName       <- str_pad(ResSCP$State_County_Place, xLen, "right")
         #         ResSCP$FCol        <- str_pad(paste0("      ",ResSCP$Place,sep=""),MaxColumn,"right")
         # 
         #         # Get tracts in cluster from xx2.
         #         ResSCP$Observed    <- sprintf("%10i",xx2$Observed)
         #         ResSCP$Expected    <- sprintf("%10.3f",xx2$Expected)
         #         ResSCP$Obs_Exp     <- sprintf("%10.3f",xx2$Observed / xx2$Expected)
         #         ResSCP[is.nan(ResSCP$Obs_Exp),"Obs_Exp"] = 0
         #   
         #         # Get total tracts in St/Co/Place area from PL99 table
         #    
         #         ResSCP$TotTracts   <- 0
         #         ResSCP$TotTracts   <- pl99.index[ResSCP$plKey,TractN] 
         #       
         #         ResSCP$PrcTracts   <- sprintf("%10.3f%%",ResSCP$Tracts_In/ResSCP$TotTracts*100)
         #         ResSCP$Tracts_In   <- sprintf("%10i",ResSCP$Tracts_In)
         #       
         #         ResSCP$TotTracts   <- sprintf("%10i",ResSCP$TotTracts)
         #                 
         #         SCPHeader          <- paste0(str_pad("State/County/Place",xLen,"right"),BasicHeaderL)
         #         SCPHeaderT         <- paste0(str_pad("      Place",MaxColumn,"right"),BasicHeaderL)
         #
         #      }
         #   } 
         #   #
         #   #  data.frames are all setup for the report.
         #   #
         #   #  All list completed - now select and print                          
         #                
         #   # Tier report - State -> State/County -> State/County/Place -> Census Tract
         #     
         #   #  Tiered Summary 
         #   #
         #   #  Updated report to loop through entries instead of printing structure.
         #   #  After each State entry (Level=1), if the user wants more levels (tiers)
         #   #    of information (state/county, state/county/place and census tract), 
         #   #    they can specify the report level:  
         #   #        1=State,  
         #   #        2=State, State/County 
         #   #        3=State, State/County, State/County/Place,
         #   #        4=State, State/County, State/County/Place, Census Tract.
         #   #    "RepSt" and "RepTierLevel".
         #   #
         #   RepTier  <-  4   # all.
         ##   if (LocIDType == 2) {
          #      # minimum tier = 1
          #      RepTier = 1
          #W  }
          #  if (LocIDType == 5) {
          #      if (RepTier > 2) {
          #         # reduce report tier appropriatelys
          #         RepTier = 2
          #      }
          #  }
          #  if (LocIDType == 11) {
          #      if (RepTier > 4) {
          #         # reduce report tier appropriately
          #         RepTier = 4
          #      }
          #  }
          #  
          #  writeLines("   Tiered Summary of locations in Cluster list alphabetically by Name:",con=TxtCon)        
          #  writeLines(" ",con=TxtCon)
          #    
          #3  NumS       <- dim(ResS)[1]          # number of rows - state
          #  OrdS       <- order(ResSt$State)
          #  ResSt$Key   <- as.character(ResSt$State)
          #    
          #  # Tier 1 - State Header 
          #    
          #  writeLines(SHeaderT,con=TxtCon)
          #    
          #  for (iS in 1:NumS) {   # Loop at Tier 1
          #     # Tier 1 - State Data
          #     #cat("Text:Rep-State.\n")
          #        
          #     cat(unlist(ResSt[OrdS[iS],c("FCol","Observed","Expected","Obs_Exp","TotTracts","Tracts_In","PrcTracts")]),"\n",sep=" ",file=TxtCon,append=TRUE)
          #  
          #     if (RepTier >=2) {
          #   
          #        #cat("Text:Rep-State/County.\n")
          #     
          #        Key_S     <- ResSt[iS,"Key"]
          #     
          #        ###
          #        #  Find State/County rows that match this state
          #        ###
          #       
          #        wResSC    <- ResSC[ResSC$State == Key_S,]
          #               
          #        NumSC     <- dim(wResSC)[1]         # number of rows
          #        OrdSC     <- order(wResSC$State_County)
          #       
          #        # Tier 2 - State/County Header
          #       
          #        writeLines(SCHeaderT,con=TxtCon)
          #  
          #        for (iSC in 1:NumSC)  {   # loop at Tier 2
          #           # Tier 2 - State/County DATA
          #           #cat("Text:Rep-State/County Records:",iSC,"\n")
          #           
          #           cat(unlist(wResSC[OrdSC[iSC],c("FCol","Observed","Expected","Obs_Exp","TotTracts","Tracts_In","PrcTracts")]),"\n",sep=" ",file=TxtCon,append=TRUE)
          #        
          #           # have State/County - now find all State/County/Places within containing census tracts.
            
           #          if (RepTier >= 3) {
           #             # Tier 3 - State/County/Tract
           #          
           #             Key_S_C    <- wResSC[OrdSC[iSC],"State_County"]
           #     
           #             ###
           #             #  Find State/County/Place rows that match this state/county
           #             ###
           #      
           #             wResSCP    <- ResSCP[ResSCP$State_County == Key_S_C,]   # have rows - should be > 1
           #      
           #             # Tier 3 - State/County/Place Header
           #          
           #             writeLines(SCPHeaderT,con=TxtCon)
           #     
           #             NumSCP     <- dim(wResSCP)[1]
           #             OrdSCP     <- order(wResSCP$State_County_Place)
           #          
           #             for (iSCP in  1:NumSCP) { # Loop at Tier 3
           #                #  Tier 3 - State/County/Place DATA
           #                #cat("Text:S/C/P:",iSCP,"\n")
           #          
           #                writeLines(" ",con=TxtCon)  # blank at start of Places
           #          
           #                cat(unlist(wResSCP[iSCP,c("FCol","Observed","Expected","Obs_Exp","TotTracts","Tracts_In","PrcTracts")]),"\n",sep=" ",file=TxtCon,append=TRUE)   # get a SCP entry 
           #         
           #                ###
           #                #  Check to see if location record is wanted.
           #                ###
           #             
           #                if (RepTier >= 4) {
           #                   # Tier 4 - State/County/Place/Tract
           #               
           #                   Key_S_C_P <- wResSCP[iSCP,]$State_County_Place
           #           
           #                   wLoc      <- GisLocList1[GisLocList1$St_Co_Place == Key_S_C_P,]   # get associate rows to state/county/place
           #               
           #                   # Tier 4 - State/County/Place/Tract Header
           #               
           #                   writeLines(GisLocHeaderT,con=TxtCon)
           #               
           #                   NumLoc    <- dim(wLoc)[1]
           #               
           #                   for (iLoc in 1:NumLoc) {  # Loop at Tier 4
           #                      # Tier 4 - State/County/Place/Trace DATA
           #                      cat(unlist(wLoc[iLoc,c("FCol","Observed","Expected","Obs_Exp","RelRisk","Category")]),"\n",sep=" ",file=TxtCon,append=TRUE)                             
           #                   }
           #                }      
           #             }  # End of loop for level 3 - St/Co/Place
           #           
           #             writeLines(" ",con=TxtCon)  # blank at start of Places
           #          }       
           #       } # end of loop for Level 2 - St/Co
           #    }
           # } # end of level 1 loop    
            writeLines(" ",con=TxtCon)
            writeLines(" ",con=TxtCon)
            
            #  Wrap up reports.            
            writeLines(" ",con=TxtCon)
            writeLines(" ",con=TxtCon)
   
         } else {  
            # no locations in cluster
            writeLines("*** There are no locations found for this cluster ***",con=TxtCon)
            writeLines(" ",con=TxtCon)
            writeLines(" ",con=TxtCon)
         }             
           
      } # end of For Cluster loop  
   
   } else {
      xmsg <- paste0("*** No Cluster Data with P_VALUE < ",pValue," available in data.  ***")
      writeLines(xmsg,con=TxtCon)
      writeLines(" ",con=TxtCon)
      writeLines(" ",con=TxtCon)
   }
   
   ##########################
   close(TxtCon)
 
   print("End of text report")

   ############################ end of report ##########################

   options(width=Save_Width)


    
}  # end of function call.

#######  Function has been created and loaded.


