######
#
#  CreateGeo4SS  - Create Geo file for SatScan (TM) based on locations
#     in the .pop and .cas files and the boundary information 
#     used by satscanMapper and SeerMapper packages from NCI.
#
#  Call is geoData <- CreateGeo4SS(cas=<filename>)
#
#    The call has four call parameters:
#        path =  the directory path to locate the cas, pop and place the geo file.
#           (optional)
#        cas  =  the path/filename or filename (used with path= value) to locate the case file
#        pop  =  the path/filename or filename (used with path= value) to locate the population file
#        geo  =  the path/filename or filename (used with path= value) to specify the 
#            output location of the .geo file.
#
#        Either the cas= or pop= must be specified.  Both can be specified to provide
#        the most complete set of location IDs to build the .geo coordinate file.
#        
#        The geo file structure is returned as the value of the call.  If the geo=
#        parameter is specified, the .geo coordinate information is written to the
#        specified file.
#
#        The path= parameter is optional.  If present, it will be prepended to the 
#        filenames provided in the cas=, pop= and geo= parameters.  It the path=
#        parameter is omitted, the cas=, pop= and geo= pareameters must contains
#        any relevent path information needed.
#
#   The file extensions used on the case and population files are not check.
#   The file extension used on the geo file is specified in the geo= parameter.
#
#   The cas= and pop= filenames are check to see if they exists.  One must be valid
#   for the function to create a geo file.  A warning is generated and the bad file 
#   is ignored.   If the geo= file already exists, it will be overwritten, if overwrite=TRUE
#   in the function call.
#
#   Programming note:  All of the standard coordinate and centroid files
#   are based on the mapping coordinates for the graphic page.  This includes
#   the move of Alaska, Hawaii and PR to new geographic position on the map.
#   Since the coordinates used by SaTScan (TM) must be the real coordinates,
#   this function is provided with a set of files at the state, county and 
#   census tract level with the original coordinates (not moved) for the 
#   state, county and census tract centroids for ALASKA, HAWAII, and PR.
#   These files are st99_M_data, co99_M_data, and tr99_M_data.
#   When this function is called, the standard and non-moved files 
#   are loaded and the standard mapping files are modified to 
#   reflect the TRUE equal area projection coordinates.
#

CreateGeo4SS <- function(path  = NULL,
                    pop        = NULL,
                    cas        = NULL,
                    geo        = NULL,
                    overwrite  = FALSE,
                    header     = TRUE,
                    censusYear = NULL
                  )  {
   
   #####
   #
   #  Local Functions
   #
   #####
   #
   #
   strReverse <- function(x)
        sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
   #
   #####

   #####
   #
   #
   ParsePFE <- function(x) {
      results <- c("","","","")        # set results to "", Drive, Path, Filename-Base, extention
      if (!is.character(x)) {
         # if not character vector - error - stop
         stop("The path/file/ext string must be a character vector.")
      } else {
         # check first character of string is "." or "~" 
         fchr      <- str_sub(x,1,1)
         if (fchr == "." | fchr == "~") {
            x      <- normalizePath(x,mustWork=FALSE)
         }
         # reverse strings to allow search from right to left.
         xr        <- strReverse(x)
         
         # peel off extention  (terminators of ".", slashes, and ":"
         xes       <- regexpr("[\\.]|[//]|[\\]|[:]",xr)[1]  # get first answer (position)
         if (xes >= 0) {
           # found extention "."
           mchr       <- str_sub(xr,xes,xes)   # get matched character
           if (mchr==".") {
              #  Yes, we have an extension.
              results[4] <- str_sub(xr,1,xes)
              xr  <- str_sub(xr,xes+1)
           }
         }
         # peel off filename
         xes      <- regexpr("[\\]|[//]|[:]",xr)[1]
         if (xes >= 0) {
            # found filename
            results[3] <- str_sub(xr,1,xes-1)
            xr         <- str_sub(xr,xes)
         } else {
            # no slash or : found.  No path.
            results[3] <- xr
            xr         <- ""
         }
         # peel off path from drive
         xes           <- regexpr(":",xr)[1]
         if (xes >= 0 ) {
            # found separation of path and drive
            results[2] <- str_sub(xr,1,xes-1)
            xr         <- str_sub(xr,xes)
            results[1] <- xr  # drive
            xr         <- ""
         } else {
            # no separator of drive
            results[2] <- xr
            xr         <- ""
            results[1] <- ""
         }
      }
      results    <- strReverse(results)             # flip the results around again.
      #cat("results:",results,"\n")
      
      xPath   <- paste0(results[1],results[2],results[3],results[4])
      #cat("xPath:",xPath,"\n")
      
      if (dir.exists(xPath)) {
         #cat("The string is a path not file.  Put filename and extention back in path.\n")
         # not a file, put filename and extention back under path.
         results[2] <- paste0(results[2],results[3],results[4])
         results[3] <- ""
         results[4] <- ""
      }
      if (nchar(results[2])>0) {
         lchr       <- str_sub(results[2],-1,-1)
         #cat("results[2]:",results[2],"  lchr:",lchr,"\n")
         
         if (lchr != "/" && lchr != "\\") {
            #cat("add a slash..\n")
            results[2] <- paste0(results[2],"/")
         }
      }
      return(results)
   }
   # 
   #  end of the  ParsePFE function
   #
   #####
   
   #####
   #
   #  Initialize function variables
   #
   
   ErrFnd      <-  FALSE    # if TRUE - can't run.
   
   RV          <- NULL
   
   RV$basePath <- ""
   RV$v_Path   <- FALSE
   
   RV$popFN    <- ""
   RV$v_Pop    <- FALSE
   RV$casFN    <- ""
   RV$v_Cas    <- FALSE
   RV$geoFN    <- ""
   RV$v_Geo    <- FALSE
   
   RV$overwrite <- FALSE
   #
   #  Check call parameters
   #
   cat("Call Parameter - RAW:\n")
   cat("   path=",path,"\n")
   cat("   pop =",pop, "\n")
   cat("   cas =",cas, "\n")
   cat("   geo =",geo ,"\n")
   cat("   hdr =",header,"\n")
   cat("   censusYear =",censusYear,"\n")
   
   
   #
   #  # 1 - geo output - overwrite permission parameter
   #
   overwrite_def  <- FALSE
   RV$overwrite   <- overwrite_def
      
   if (is.null(overwrite) || is.na(overwrite) ) {
      RV$overwrite   <- overwrite_def
   
   } else {
      overwrite <- overwrite[[1]][1]
      
      if (!is.logical(overwrite)) {
         xmsg <- paste0("The overwrite call parameter is not a logical variable. The default value of FALSE will be used.")
         warning(xmsg, call.=FALSE)
         RV$overwrite   <- overwrite_def
      } else {
         # it's a logical variable
         RV$overwrite   <- overwrite
      }
   }
   
   cat("overwrite set to ",RV$overwrite,"\n")
   #
   #  # 2 - path parameter (optional)
   #        If provided used with the pop or cas filenames if the do not 
   #          stand alone as filename.
   #
   path_def <- ""
   
   if (is.null(path) || is.na(path)) {
      # no path provides
      RV$basePath <- path_def
      RV$v_Path   <- FALSE    # no path specified.
      #cat("no path seen\n")
   
   } else {
      # have a path value assigned by caller
      path <- path[[1]][1]
      
      if (is.character(path)) {
        # have a valid string, but is is valid path?
        # must be directory and not filename.
        
        if (dir.exists(path)) {
           # have good directory-path
           RV$basePath <- path
           RV$v_Path   <- TRUE
           
           lastChr <- str_sub(RV$basePath,-1,-1)    # get last character of path
           
           # see if "/" should be added
           if (!(lastChr == "\\" || lastChr == "/")) {
              RV$basePath <- paste0(RV$basePath,"/")   # no slash on end of path - add one
           }
           #cat("Good path=",RV$basePath,"\n")
        } else {
           # not a directory/path, did caller specify a file?
           if (file.exists(path)) {
              
              # user specified a file not a directory/path
              # throw error and ignore.
              xmsg <- paste0("The path parameter value provided is a file and not a directory/path. It cannot be used for a path to cas and pop files.")
              save_FN  <- path                  # save it for later?
              warning(xmsg,call.=FALSE)
              
           } else {
              # not a path or file - throw error and ignore.
              xmsg <- paste0("The path parameter value is not a valid directory/path.")
              warning(xmsg,call.=FALSE)
              #cat("path not found:",path,"\n")
              #cat("dir.exists(path):",dir.exists(path),"   file.exists(path):",file.exists(path),"\n")
              
           }
           #cat("path not found:",path,"\n")
           #cat("dir.exists(path):",dir.exists(path),"   file.exists(path):",file.exists(path),"\n")
           #  ignore parameter
           RV$basePath <- path_def
           RV$v_Path   <- FALSE
        }
      } else {      
        # not a character vector - throw error and ignore.
        xmsg <- paste0("The path parameter provided is not a character vector.  The path parameter will be ignored.")
        warning(xmsg, call.=FALSE)
        RV$basePath  <- path_def
        RV$v_Path    <- FALSE
      }
   }
   
   #
   #  # 3 - cas - case file parameter (can also be used for pop file location)
   #
   casFN_def     <- ""
  
   if (is.null(cas) || is.na(cas)) {          # Was parameter set by caller?
      # parameter not provided, set to "" (no name)
      RV$casFN   <- casFN_def                 # no - set defaults
      RV$v_Cas   <- FALSE
      cat("no cas parameter seen.\n")
   } else {
      cas        <- cas[[1]][1]
   
      # yes, caller provide cas=; is it a character vector?
      if (is.character(cas)) {
         #
         # test to see if cas points to an existing file. 
         #    That means the filename provide works without the use of the "path"
         #
         if (file.exists(cas)) {
            # cas contains a valid filename - it can stand on it's own.
            RV$casFN  <- cas
            RV$v_Cas  <- TRUE
            #cat("cas good on it own:",RV$casFN,"\n")
         
         } else {
            # not valid existing filen - check to see if it works with path= parameter added?
            # if path= exist?
            if (!RV$v_Path) {
               # no path provided - throw error - and ignore.
               #cat("The cas= filename does not exist. The cas= parameter is ignored.\n")
               RV$casFN  <- casFN_def
               RV$v_Cas  <- FALSE
            } else {
               # path= was provided and was valid.
               # cas is not a valid path/filename on it's own  - try with path added
               wStr      <- paste0(RV$basePath,cas)
               #cat("path&cas=",wStr,"\n")
               
               if (file.exists(wStr)) {
                  # path and cas point to a good existing file
                  RV$casFN   <- wStr      # save full path/file
                  RV$v_Cas   <- TRUE
               } else {
                  # path= and cas= don't point to an existing file - throw error and ignore.
                  xmsg <- paste0("The path= and cas= parameters combined are not an existing cas file.  The cas= parameter ignored.")
                  warning(xmsg,call.=FALSE)
                  xmsg <- paste0("  path=",RV$basePath,"  and cas=",cas)
                  warning(xmsg,call.=FALSE)
                  RV$casFN   <- casFN_def
                  RV$v_Cas   <- FALSE
               }
            } 
         }
      } else {
         # not character - throw error and ignore.
         xmsg <- paste0("The cas parameter is not a character vector. The cas parameter is ignored.")
         warning(xmsg, call.=TRUE)
         RV$casFN  <- casFN_def
         RV$v_Cas  <- FALSE
      }
   }
   
   #
   #  # 4 - pop - population file parameter (can also be used for cas file location)
   #
   popFN_def   <- ""
   
   # was the pop= call parameter specified by caller?
   if (is.null(pop) || is.na(pop)) {
      # parameter not provided, set to "" (no name)
      RV$popFN   <- popFN_def    # set defaults - no pop file
      RV$v_Pop   <- FALSE
   } else {
      #  caller provided pop= parameter
      pop <- pop[[1]][1]
   
      if (is.character(pop)){    # is it a character vector?
         #
         # test to see if cas points to an existing file. 
         #    That means the filename provide works without the use of the "path"
         #
         if (file.exists(pop)) {
            # pop contains a valid filename - it can stand on it's own.
            RV$popFN  <- pop
            RV$v_Pop  <- TRUE
            
         } else {
            # not valid filename - check to see if it works with path added?
            # does path exist?
            if (!RV$v_Path) {
               # no path provided - throw error
               xmsg <- paste0("The pop= filename does not exist. The pop= parameter is ignored.")
               warning(xmsg, call.=FALSE)
               
               RV$popFN  <- popFN_def
               RV$v_Pop  <- FALSE
            } else {
               # pop is not a valid path/filename on it's own  - try with path added
               wStr <- paste0(RV$basePath,pop)
               #cat("path&pop=",wStr,"\n")
               
               if (file.exists(wStr)) {
                  # path and pop are a file reference - good
                  RV$popFN   <- wStr      # save full path/file
                  RV$v_Pop   <- TRUE
               } else {
                  xmsg <- paste0("The path= and pop= parameters combined are not an existing cas file.  The pop= parameter ignored.")
                  warning(xmsg,call.=FALSE)
                  xmsg <- paste0("  path=",RV$basePath,"  and pop=",pop)
                  warning(xmsg,call.=FALSE)
                  RV$popFN   <- popFN_def
                  RV$v_Pop   <- FALSE
               }
            } 
         }
      } else {
         # not character
         xmsg <- paste0("The pop parameter is not a character vector. The pop parameter is ignored.")
         warning(xmsg, call.=TRUE)
         RV$popFN     <- popFN_def
         RV$v_Pop     <- FALSE
      }
   }
   
   #
   #  # 5 - geo parameter
   #
   geoFN_def   <- ""
   RV$geoWrite <- FALSE
   RV$v_Geo_def<- FALSE
   
   # did caller provide a geo= parameter?
   if (is.null(geo) || is.na(geo)) {
      # parameter not provided, set to ""
      RV$geoFN   <- geoFN_def     # set defaults = no GEO
      RV$v_Geo   <- RV$v_Geo_def
      #  notify caller, no geo= parameter, so function will use the cas (then pop) filename
      #  with the extension of .geo for the generated .geo file.
      
      #xmsg       <- paste0("No output geo= filename has been specified. ")
      #warning(xmsg,call.=FALSE)
      
      #WhoHelped <- ""
      #if (RV$v_Cas)  {
      #   pFN       <- ParsePFE(RV$casFN)
      #   RV$geoFN  <- paste0(pFN[1],pFN[2],pFN[3],".geo")
      #   RV$v_Geo  <- TRUE
      #   WhoHelped <- "cas"
      #} else {
      #   # no cas filename - see if pop exists
      #   if (RV$v_Pop) {
      #      pFN      <- ParsePFE(RV$popFN)
      #      RV$geoFN <- paste0(pFN[1],pFN[2],pFN[3],".geo")
      #      RV$v_Geo <- TRUE
      #      WhoHelped <- "pop"
      #   }
      #}
      #if (!RV$v_Geo) {
         ## a geo filename was created from the cas or pop filename. (not permitted,  If no geo file then return variable.)

         xmsg      <- paste0("No geo= specified. Coordinates data are returned as the value of this call.")
         warning(xmsg,call.=FALSE)

         #
         #  still need to check for possible overwrite condition
      #}
   } else {
      # caller did provide geo= parameter
      geo <- geo[[1]][1]
      
      if (is.character(geo)) {   # is it a character vector?

         # is a character string

         pFN  <- ParsePFE(geo)     # parse the filename into components.

         # test to see if the geo file already exists?
         
         if (file.exists(geo)) {
         
            # geo contains a valid existing filename, save info, check for overwrite later.
            RV$geoFN   <- geo
            RV$v_Geo   <- TRUE
         
         } else {
            #  files does not exist - good
            #  is the path part of the filename good?
            geoPath          <- paste0(pFN[1],pFN[2])   # get drive and path
            if (geoPath != "") {
               # geoPath is not empty, is it good?
               if (dir.exists(geoPath)) {
                  # path in geo string exists.  So assume file is ok and will be created
                  RV$geoFN   <- geo             # save path/filename
                  RV$v_Geo   <- TRUE
               } else {
                  # geo path part is not valid
                  xmsg <- paste0("The geo= call parameter does not contain a valid directory/path. ")
                  warning(xmsg, call.=FALSE)
                  xmsg <- paste0("The geo output file cannot be accessed: ", geoPath, " The geo= parameter is ignored.")
                  warning(xmsg, call.=FALSE)
                  RV$geoFN   <- geoFN_def
                  RV$v_Geo   <- FALSE
               }
            } else {
               # no path part in the geo= parameter.  geo is only the filename..
               if (!RV$v_Path) {
                  # no extra path= available and geo does not exist.
                  #  Must assume the file can be created when it is written.
                  RV$geoFN     <- geo    # only geo provided to create output path/filename.
                  RV$v_Geo     <- TRUE
               
               } else {
                  # have path= (basePath) - attach it to the geo filename. 
                  wStr <- paste0(RV$basePath,geo)
                  # check to see if the path= and geo= combinate exist as a file.
                  RV$geoFN   <- wStr
		  RV$v_Geo   <- TRUE
               }
            }
         }
      } else {
         # not character
         xmsg <- paste0("The geo parameter is not a character vector. The geo parameter is ignored.")
         warning(xmsg, call.=TRUE)
         RV$geoFN    <- geoFN_def
         RV$v_Geo    <- FALSE
      }
   }
   #
   #  almost done with geo parameter - if filename set, we must check to see if 
   #  it exists and if it can be overwritten.
   #
   if (RV$v_Geo) {
      # we have or created an GEO output filename.  does it exist? and can we overwrite?
      if (file.exists(RV$geoFN)) {
         xmsg <- paste0("The geo= path/filename of (",RV$geoFN,") already exists.")
         warning(xmsg, call.=FALSE)
         # yes = output filename existx.
         if (RV$overwrite) {
            xmsg <- paste0("The overwrite= call paramter is set to TRUE, so the file will be replaced/overwritten.")
            warning (xmsg, call.=FALSE)
            RV$geoFN        <- geo
            RV$v_Geo        <- TRUE
            RV$geoWrite     <- TRUE
         } else {
            xmsg <- paste0("Overwriting of the existing geo file is not permitted. No output file will be created.")
            warning(xmsg, call.=FALSE)
            RV$geoFN        <- geoFN_def
            RV$v_Geo        <- FALSE
            RV$geoWrite     <- FALSE
         }
      } else {
         # the file does not exist - Write it.
         RV$geoWrite        <- TRUE
      }
   }
   
   #
   #  # 6 - censusYear
   #
   censusYear_def <-  "2000"
   
   if (is.null(censusYear) || is.na(censusYear)) {
      # parameter not provided, set to ""
      RV$censusYear  <- censusYear_def
   } else {
      censusYear <- censusYear[[1]][1]
      
      censusYear <- as.character(censusYear)   # convert to character string
      # we actually do accept numeric and integer input. HA!
      
      if (censusYear != "2000" && censusYear != "2010") {
         # census year is not set to 2000 or 2010.  Invalid use default
         xmsg <- paste0("The censusYear parameter is not set to '2000' or '2010'.  The default census year of '",censusYear_def,"' is used.\n")
         RV$censusYear <- censusYear_def
      } else {
         RV$censusYear <- censusYear
      }
   }
   
   #
   #  # 7 - header on cas or pop files
   #
   header_def <- TRUE
   
   if (is.null(header) || is.na(header) ) {
      header <- header_def
   } else {
      if (!is.logical(header)) {
         # wrong type
         xmsg <- paste0("The header= parameter is not a logical variable.  The default of ",header_def," will be used.\n")
         warning(xmsg,call.=FALSE)
         header <- header_def
      }
   }
   RV$header <- header
   
   #
   #  Time to cross check the input parameters.
   # 
   
   #
   #  at least one cas or pop must be accessible.
   #
   if (!(RV$v_Cas || RV$v_Pop)) {
      xmsg <- paste0("At least one valid input SaTScan cas or pop file must exist.  Can not create geo coordinates file.")
      stop(xmsg,call.=TRUE)
   }
   
   #
   #   check for other errors that keep us from running.
   #
   if (ErrFnd) {
      xmsg <- paste0("Errors were found during the validation of the call parameters.  Please correct.")
      stop(xmsg,call.=TRUE)
   }
  
   #
   #  if no Geo filename, the data.frame is returned at the end of the call.
   #
  
   #
   #  Read cas and/or pop file.
   #
   
   #
   #  Input format:  
   #    CAS File:  <location ID>  <case-count>  (date/time> <attributes> <censored> <weight> <covariates>
   #      We are only interested in the <location ID> field. The rest of the fields are ignored. 
   #      According to SaTScan manual, the location ID is any numerical or string of characters.
   #      Empty spaces may not form part of the id. (example: New York is not permitted.
   #      Quoted strings are not supported. 
   #
   #    POP File:  <location ID> <date/time> <population> <covariates> 
   #      We are only interested in the <location ID> field.  The rest of the fields are ignored.
   #
   #    In both cases, it's only the first field that is of interest.
   #    The delimiter can be a comma, space or tab. The field could be in quotes.
   #    
   #    Validate tests:  numeric (0-9) and proper number of digits.
   #       state              = 1 or 2 characters
   #       state/county       = 4 or 5 characters
   #       state/county/tract = 8, 9, 10 or 11 characters
   #      
   #       if 1, 4, 8, 10 characters = add leading "0"
   #       if 8 or 9 characters - add trailing two "0".
   #
   #       Invalid IDs = 0, 3, 6, 7 and > 11 characters.
   #
   #       After edit corrections, 2 = state, 5 = county, 11 = tract.
   #
   
   #
   #  Write geo file.
   #
   #  Output format:
   #    GEO File:  <location ID> <x> <y>
   #
   #    The cartisian coordinates provided in the GEO file are a translation from the Lat/Long 
   #    coordinates based on the proj4 parameter: "+proj=longlat +datum=NAD83" to an Albers
   #    equal area projection based on the proj4 parameters: 
   #    "+proj=eqdc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=96w +units=m" with the unit of measure
   #    of meters.
   #
   #    Any other coordinates used for the Grid File, 
   #
   # The files must be ASCII formatted files with space, tab or comma delimitered fields.
   # 
   #   OrigCRS      <- CRS("+proj=longlat +datum=NAD83")
   #
   #  Transform the State, State/County, State/County/Census Tract
   #  boundary polygons from long/lat to Equidistance Conic projection.
   #
   #   Projection = Equidistance-Conical => simpleconic
   #   Lat Parallel 1   = 33
   #   Lat Parallel 2   = 49   (should be 45.
   #   Origin of Lat    = 39
   #   central Meridian = -96   (96W)
   #
   #    lat_1=33,   lat_2=49  (dif = 16, center 41)  lat_0=39 not 41
   #    lat_1=33,   lat_2=45  (dif = 12, center 39)
   #    lat_1=29.5, lat_2=45.5(dif = 16, center 37.5)  lat_0=37.5
   #
   #  ProjCRS <- CRS("+proj=eqdc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=96w +units=m")
   #
   # ESRI -> +proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs 
   #
   
   casTable <- NULL
   popTable <- NULL
   lnames   <- c("LOC_ID","P1","P2","P3","P4","P5","P6","P7","P8","P9")
   
   #
   #  Read cas and pop files.
   #
   if (RV$v_Cas) {
      # read .cas file
      #cat("Reading cat:",RV$casFN,"\n")
      casTable  <- read.table(RV$casFN, header=RV$header, colClasses="character",
                              blank.lines.skip=TRUE,
                              )
      nR        <- dim(casTable)[2]     # number of rows
      names(casTable) <- lnames[1:nR]
   }
   if (RV$v_Pop) {
      # read .pop file
      #cat("Reading pop:",RV$popFN,"\n")
      popTable  <- read.table(RV$popFN, header=RV$header, colClasses="character",
                              blank.lines.skip=TRUE
                              )
      nR        <- dim(popTable)[2]
      names(popTable) <- lnames[1:nR]
   }
   #
   #  Collect, combine and fine unique location IDs
   #
   
   LocList  <- c(casTable$LOC_ID,popTable$LOC_ID)
   ULocList <- str_trim(unique(LocList))
      
   #cat("List of unique LocList:\n")
   #print(ULocList)
 
   x <- ULocList == ""
   y <- is.na(x)
   x <- x | y
   
   if (any(x)) {
      xmsg <- paste0("Empty or NA value locations IDs found in the data.  Entries ignored.")
      warning(xmsg,call.=FALSE)
      ULocList <- ULocList[!x]
   }
   
   xR <- grepl("^[0-9]{1,11}$",ULocList)
   
   if (!all(xR)) {
      # we found a non-numericf
      BadList <- ULocList[!xR]
      xmsg <- paste0("At least of the location IDs is non-numeric. is empty, or greater than 11 characters.")
      warning(xmsg, call.=FALSE)
      xmsg <- paste0("  Invalid Location IDs: ",paste0(BadList,collapse=" "))
      warning(xmsg, call.=FALSE)
      xmsg <- paste0("Entry is ignored and run will continue, but please investigate and correct.")
      warning(xmsg, call.=FALSE)
      ULocList <- ULocList[xR]
      ##  add quots for empty entries... or something.
   }
   
   # classify the location ID type
   nLocList     <- nchar(ULocList)
   nLocRange    <- range(nLocList)
   
   #  Classify the location IDs   (State, County, Census Tracts)  1-2, 4-5, 8-11
   #
   #    S              1
   #    SS           2
   #
   #    SCCC           4
   #    SSCCC        5
   #
   #    SCCCTTTT       8
   #    SSCCCTTTT    9
   #    SCCCTTTTTT    10
   #    SSCCCTTTTTT 11
   #
   #    invalid = <0, 0, 3, 6, 7, and > 11
   #
   #  Fix FIPS codes if needed.
   #
   ErrFnd <- FALSE
   tVL <- nLocRange[1]
   tVH <- nLocRange[2]
   
   TestLocID <- function(x) {
          res <- TRUE
          lx  <- nchar(x)
          if (lx <= 0 ) res <- FALSE
          if (lx == 3 ) res <- FALSE
          if (lx == 6 ) res <- FALSE
          if (lx == 7 ) res <- FALSE
          if (lx > 11 ) res <- FALSE
          return (res)
   }
   TestLocList <- function(y) {
      sapply(y, function(z) TestLocID(z))
   }
   
   
   xR <- switch(tVH,
                ifelse (tVL == 1, FALSE, TRUE),
                ifelse ((tVL == 1 || tVL == 2), FALSE, TRUE),
                TRUE,
                ifelse (tVL == 4, FALSE, TRUE),
                ifelse ((tVL == 4 || tVL == 5), FALSE, TRUE),
                TRUE,
                TRUE,
                ifelse (tVL == 8 ,FALSE, TRUE),
                ifelse ((tVL >= 8 && tVL <= tVH), FALSE, TRUE),
                ifelse ((tVL >= 8 && tVL <= tVH), FALSE, TRUE),
                ifelse ((tVL >= 8 && tVL <= tVH), FALSE, TRUE),
               TRUE
              )
              
   gLocList <- TestLocList(ULocList)
   
              
   #cat("xR:",xR,"  tVH:",tVH,"  tVL:",tVL,"\n")
   
   if (xR) {
      # the location ID are the wrong format.
      BadList <- ULocList[!gLocList]
      xmsg <- paste0("The location IDs range from ",tVL, " to ",tVH," digits and are not valid FIPS codes. Bad IDs ignored.")
      warning(xmsg,call.=FALSE)
      xmsg <- paste0("  Bad location IDs:",paste0(BadList,collapse=" "))
      warning(xmsg,call.=FALSE)
      
      ##  need to quote the error list incase a "" is present.
      
      ULocList <- ULocList[gLocList] 
   }
   
   nLocList <- nchar(ULocList)
   
   #
   #  fix up the location IDs
   #
   #  Leading zero for the 1 digit state codes.
   sAdj <- c(1,4,8,10)
   xM   <- !is.na(match(nLocList,sAdj))
   if (any(xM)) {
     ULocList[xM]  <- paste0("0",ULocList[xM])
   } 
   
   #  Trailing 2 digits for the short tract codes.  missing 1 trailing digit should not happen.
   tAdj       <- c(8,9)
   xM         <- !is.na(match(nLocList,tAdj))
   
   if (any(xM)) {
      ULocList[xM]  <- paste0(ULocList[xM],"00")
   }
   
   nLocList  <- nchar(ULocList)
   nLocRange <- range(nLocList)
   
   #
   #  at this point everyone should pass the test
   #
  
   xR <- TestLocList(ULocList)
   
   #cat("xR after edits and final test:",all(xR),"\n")

   if (nLocRange[1] != nLocRange[2]) {
      #  Not all of the Location IDs are at the same level.
      xmsg <- paste0("Not all of the Loc_IDs are at the same FIPS code level.  Check the Loc_ID and re-run.\n")
      warning(xmsg,call.=FALSE)
      
      x <- unique(nLocList)
      xmsg <- paste0("The function has found Loc_ID codes  with ",x," number of digits.\n")
      stop(xmsg, call.=FALSE)
   }   
   
   idMode <- nLocRange[2]   # get number of digits and use as idMode.
   
   idName <- switch(idMode,
                        "INVALID",    # 1
   			"STATE",      # 2
   		        "INVALID",    # 3 
   		        "INVALID",    # 4
   			"COUNTY",     # 5
   			"INVALID",    # 6
   			"INVALID",    # 7
   			"INVALID",    # 8
   			"INVALID",    # 9
   			"INVALID",    # 10
                        "CENSUS TRACT",# 11
   			"INVALID"
   	           )
   	           
   #cat("idMode:",idMode,"  idName:",idName,"\n")
   
   #cat("Centroids will be collected for the Location IDs at the",idName,"level.\n")   
   
   #cat("Evaluate Location ID, Select Level, and get centroid values.\n")
   
   #cat("Always load state tables.\n")
   data(st99_data,  envir = environment())
   data(st99_M_data,envir = environment())
   st99_data$ID      <- row.names(st99_data)
   st99_M_data$ID    <- row.names(st99_M_data)
   
   #cat("Update to TRUE coordinates\n")
   stIDs             <- st99_M_data$ID
   st99_data[stIDs,] <- st99_M_data[stIDs,]  # replace original with scewed centroid.
   
   #cat("Validate state list from location IDs\n")
   
   stList     <- sort(unique(str_sub(ULocList,1,2)))
   # we will return the centroids for all states referenced in the data location ids
   
   xM         <- match(stList,st99_data$ID)  # is if all stIDs are good.
   # need to validate stList.
   
   xMN        <- is.na(xM)    
   
   if (any(xMN)) {
      # have one or more bad state ids.
      BadList <- stList[xMN]
      xmsg <- paste0("Found one or more invalid state IDs in the location IDs: ",paste0(BadList,collapse=" "))
      stop(xmsg,call.=FALSE)
   }
   
   #cat("Get centroids based on Location ID type (State, county, tract\n")
   if (idMode == 2) {
      #cat("state\n")
      xM          <- match(st99_data$ID,stList)
      xFnd        <- !is.na(xM)
      
      LocCentroid <- st99_data[xFnd,c("c_X","c_Y")]
      #print(LocCentroid)
    }
   
   
   if (idMode == 5) {
      #cat("county - load all county data\n")
      
      data(co99_data,envir = environment())
      data(co99_M_data,envir = environment())
      co99_data$ID      <- row.names(co99_data)
      co99_M_data$ID    <- row.names(co99_M_data)
      
      #cat("Update county to TRUE coordinates.\n")
      coIDs             <- co99_M_data$ID
      co99_data[coIDs,] <- co99_M_data[coIDs,]  # replace original with screwd centroids.
      
      co99_data$stID    <- as.character(co99_data$stID)
      co99_M_data$stID  <- as.character(co99_M_data$stID)
      
      xM                <-   match(co99_data$stID, stList)
      xFnd              <- !is.na(xM)
      
      if (RV$censusYear == "2000") {
         LocCentroid    <- co99_data[xFnd,c("c_X_00","c_Y_00","y")]
         xMM            <- LocCentroid$y == 1 | LocCentroid$y == 3
      } else {
         LocCentroid    <- co99_data[xFnd,c("c_X_10","c_Y_10","y")]
         xMM            <- LocCentroid$y == 2 | LocCentroid$y == 3
      }
      
      LocCentroid       <- LocCentroid[xMM,]
      colnames(LocCentroid) <- c("c_X","c_Y","y")
      LocCentroid$y     <- NULL
   }
   if (idMode == 11) {
      #cat("tract - load tract data\n")
      data(tr99_data,envir = environment())
      data(tr99_M_data,envir = environment())
      tr99_data$ID     <- row.names(tr99_data)
      tr99_M_data$ID   <- row.names(tr99_M_data)
   
      #cat("Update to TRUE centroid coordinates.\n")
      trIDs            <- tr99_M_data$ID
      tr99_data[trIDs,]<- tr99_M_data[trIDs,]
      
      tr99_data$stcoID <- as.character(tr99_data$ID)
      tr99_data$stID   <- str_sub(tr99_data$ID,1,2)
      
      xM               <- match(tr99_data$stID,stList)
      xFnd             <- !is.na(xM)
      
      if (RV$censusYear == "2000") {
         LocCentroid   <- tr99_data[xFnd,c("c_X_00","c_Y_00","y")]
         xMM           <- LocCentroid$y == 1 | LocCentroid$y == 3
      } else {
         LocCentroid   <- tr99_data[xFnd,c("c_X_10","c_Y_10","y")]
         xMM           <- LocCentroid$y == 2 | LocCentroid$y == 3
      }
   
      LocCentroid     <- LocCentroid[xMM,]
      colnames(LocCentroid) <- c("c_X","c_Y","y")
      LocCentroid$y         <- NULL
   }
   
   #head(LocCentroid,20)
   
   LocCentroid$LocID  <- row.names(LocCentroid)
   IDList             <- row.names(LocCentroid)
   LocCentroid        <- LocCentroid[,c("LocID","c_X","c_Y")]
   #print(str(LocCentroid))
   #print(IDList)
  
   xM   <- match(ULocList,IDList)    # see if the location IDs match the boundary tables
   xMN  <- is.na(xM)
   if (any(xMN)) {
     xmsg <- paste0("The following Locations IDs in the cas and/or pop data do not match any FIPS codes for the census year:")
     warning(xmsg,call.=FALSE)
     BadList <- ULocList[xMN]
     xmsg <- paste0("   Bad Location IDs:",paste0(BadList,collapse=" "))
     warning(xmsg,call.=FALSE)
     xmsg <- paste0("Recommend research the cause and correcting.")
     warning(xmsg,call.=FALSE)
   }
   #
   #  Match up the full Location ID list to the propriate ID table.
   # 

   #
   #  Convert Centroids to equal area projection
   #
   #OrigCRS   <- CRS("+proj=longlat +datum=NAD83")
   #
   #  All centroid coordinate values (X, Y) are in the following 
   #
   #ProjCRS   <- CRS("+proj=eqdc +lat_1=33 +lat_2=49 +lat_0=39 +lon_0=96w +units=m")

   if (RV$v_Geo) {
       #cat("Writing Geo File:",RV$geoFN,"\n")
       write.table(LocCentroid, file=RV$geoFN, sep="  ",
            row.names=FALSE,
            col.names=RV$header)
   }
   
   return(LocCentroid)  # return data.frame to caller.   
}  


