#-------------------------------------------------------------------------------
#       fleets.om(fleets, biols, fleets.ctrl, year, season,... )
# 
#   Output: Updated FLFleets 
# 
#   Projects the bfleets from: - [year-1, ns] to [year,1] 
#                          or, - [year, season-1] to [year,season]
# 
#       1. Calculate effort. [seasonal].
#       2. Update landings and discards (total and at age).  [seasonal].
#       3. Calculate price. [seasonal].
#       4. Calculate capacity and catchability [annual].
#
#
# Dorleta Garcia
# Created: 01/12/2010 16:24:56
# Changed: 09/12/2010 08:36:59
#-------------------------------------------------------------------------------

fleets.om <- function(fleets, biols, GDGTs, SRs, BDs, covars, advice, biols.ctrl, fleets.ctrl, advice.ctrl, year, season){
   
    flnms <- names(fleets)
    
    fleets         <- unclass(fleets) # unclass fleets to speed up the algorithm
    fleets.ctrl.aux <- fleets.ctrl
    
    # 1. Calculate effort. 
    print('~~~~~~~~~~ EFFORT ~~~~~~~~')
    
    for(fl in flnms){
    
    print(fl)
        
        dyn.model <- fleets.ctrl[[fl]]$effort.model
        
        # # For selective fisheries:
        # if (dyn.model=='MaxProfitSeq' & is.null(fleets.ctrl[[fl]]$q2zero)) 
        #   stop(paste("fleets.ctrl[['",fl,"']]$q2zero is missing.",sep=""))
        
        if (!is.null(fleets.ctrl[[fl]]$q2zero)) 
          fleets <- catchability2zero(fleets = fleets, flnm = fl, advice = advice, fleets.ctrl = fleets.ctrl, year = year)
             
        res <- eval(call(dyn.model, biols = biols, fleets = fleets, BDs = BDs, flnm = fl, advice = advice,
                    year = year, season = season, biols.ctrl=biols.ctrl, fleets.ctrl = fleets.ctrl, covars = covars, advice.ctrl = advice.ctrl)) 
         
        fleets[[fl]]         <- res$fleets[[fl]]
        fleets.ctrl.aux[[fl]] <- res$fleets.ctrl[[fl]]
        remove(res)
    }
    
    fleets.ctrl <- fleets.ctrl.aux
    
    # 2. Update landings and discards (total and at age).  [seasonal].

    print('~~~~~~ UPDATE CATCH ~~~~~~')

    catchRet <- updateCatch(fleets, biols, GDGTs, SRs, BDs, advice, biols.ctrl, fleets.ctrl, advice.ctrl, year = year, season = season, covars = covars)
    
    ## For gadget, we need also to update biols
    fleets <- catchRet$fleets
    biols <- catchRet$biols
    SRs   <- catchRet$SRs
    # For gadget
    GDGTs <- catchRet$GDGTs

    #  3. Calculate price. [seasonal].
    print('~~~~~~~~~~ PRICE ~~~~~~~~~~')

    fleets <- unclass(fleets)
    for(fl in flnms){
    print(fl)
        
        sts <- catchNames(fleets[[fl]])
        
        for(st in sts){
            dyn.model <- fleets.ctrl[[fl]][[st]]$price.model  
        
            res <- eval(call(dyn.model, fleets = fleets, flnm = fl, stnm = st, year = year, season = season, fleets.ctrl = fleets.ctrl, covars = covars)) 
         
            fleets[[fl]]         <- res[[fl]]
        }
    }
    fleets <- FLFleetsExt(fleets)
    
    # 4. CAPITAL DYNAMICS: Calculate capacity and catchability [annual].
    if(season == dim(biols[[1]]@n)[4]){
       print('****************************** CAPITAL ******************************')

        for(fl in flnms){
            print(fl)
        
            dyn.model <- fleets.ctrl[[fl]]$capital.model  
        
            res <- eval(call(dyn.model,  fleets = fleets, advice = advice,
                flnm = fl, year = year, season = season, fleets.ctrl = fleets.ctrl, covars = covars)) 
         
            fleets <- res$fleets
            covars <- res$covars
        }
    }
    
    fleets <- FLFleetsExt(fleets)
    
    return(list(fleets = fleets, biols=biols, SRs = SRs, GDGTs = GDGTs, covars = covars, fleets.ctrl = fleets.ctrl))

}










