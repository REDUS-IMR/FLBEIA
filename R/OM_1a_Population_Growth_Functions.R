#-------------------------------------------------------------------------------
#                   POPULATION GROWTH FUNCTIONS             
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMARK: '...' in the arguments of the functions are neccesary in order to be
#   generalistic inside 'biols.om' function ('eval(call(...)').
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#------------------------------------------------
# ** ASPG(biol, SR, fleets, biol.control) ~  Age Structured Population Growth.
# Projects in one season an age structured population using
#    << N[s] = (N[s-1]*exp(-M[s-1]/2) - Catch[s-1])*exp(-M[s-1]/2)  >>
# The assumption is that natural mortality is constant and continious 
# in each season but the catch (fishing mortality) occurs instantaneously 
# in the middle of 
# the season. 
#------------------------------------------------
#
#------------------------------------------------
# # ** BDPG(biol, BD, fleets, biol.control, ...) ~  Biomass Dynamic Population Growth.
# Projects in one season a biomass dynamic population using
#    << B[s] = B[s-1] - C[s-1] + Production[s-1]  >>
# Production is determined using BDsim. 
#------------------------------------------------
#
# Dorleta Garcia
# Created: 21/10/2010 10:08:16
# Changed: 27/10/2010 08:34:17
#-------------------------------------------------------------------------------

fixedPopulation <- function(biols, SRs, fleets, year, season, stknm, ...)  return(list(biol = biols[[stknm]]))

#-------------------------------------------------------------------------------
# ASPG(biol, SR, fleets, biol.control)
# - OUTPUT: list(biol = biol, SR = SR) - Upadated FLBiol and FLSRsim objects.
#-------------------------------------------------------------------------------

ASPG <- function(biols, SRs, fleets, year, season, stknm, ...){

    cat('-----------------ASPG-----------\n')

    if(length(year) > 1 | length(season) > 1)
        stop('Only one year and season is allowed' )
    
    biol <- biols[[stknm]]     
    SR   <- SRs[[stknm]]
        
    dimnms <- dimnames(biol@n)
    
    # If year/season/iter numerics => indicate position 
    # else names => get positions.
    
    if(length(year) > 1 | length(season) > 1)
        stop('Only one year and season is allowed' )
    
    # 'year' dimension.
    yr <- year
    if(is.character(year)) yr <- which(dimnms[[2]] %in% year)
    if(length(yr) == 0) stop('The year is outside object time range')  
    
    # 'season' dimension.
    ss <- season
    if(is.character(season)) ss <- which(dimnms[[4]] %in% season)
    if(length(ss) == 0) stop('The season is outside object season range')  
    
    na <- dim(biol@n)[1]
    ns <- dim(biol@n)[4]
    stock <- biol@name
    
    # IF season = 1 THEN age groups move to the next. 
    if(ss == 1){
        # total catch in year [y-1] season [ns].
        catch.n <- catchStock(fleets,stock)[,yr-1,,ns]
        # middle ages
        biol@n[-c(1,na),yr,,ss] <- (biol@n[-c(na-1,na),yr-1,,ns]*exp(-biol@m[-c(na-1,na),yr-1,,ns]/2) - catch.n[-c(na-1,na),])*
                                            exp(-biol@m[-c(na-1,na),yr-1,,ns]/2) 
        # plusgroup
        biol@n[na,yr,,ss]       <- (biol@n[na-1,yr-1,,ns]*exp(-biol@m[na-1,yr-1,,ns]/2) - catch.n[na-1,])*exp(-biol@m[na-1,yr-1,,ns]/2) + 
                                   (biol@n[na,yr-1,,ns]*exp(-biol@m[na,yr-1,,ns]/2) - catch.n[na,])*exp(-biol@m[na,yr-1,,ns]/2)

    }
    else{
        # total catch in year [yr] season [ss-1].
        catch.n <- catchStock(fleets,stock)[,yr,,ss-1]
        # middle ages      # for unit == ss  and age = 1, it will be equal NA but be updated after with SRsim.
        biol@n[,yr,,ss] <- (biol@n[,yr,,ss-1]*exp(-biol@m[,yr,,ss-1]/2) - catch.n)*exp(-biol@m[,yr,,ss-1]/2) 
    }

    # Update SSB.
    SR@ssb[,yr,,ss] <- unitSums(quantSums(n(biol) * wt(biol) * fec(biol)*mat(biol) * 
                                            exp(-biol@m*spwn(biol)), na.rm=TRUE))[,yr,,ss]

    # RECRUITMENT
    if(dim(biol@n)[3] > 1 & dim(biol@n)[3] == dim(biol@n)[4]){
        # 'number_of_units > 1' => Recruitment occurs in all seasons, the 'unit' correspond with the recruitment 'season'.
        SR <- SRsim(SR, year = yr, season = ss, iter = 'all') 
        biol@n[1,yr,ss,ss] <- SR@rec[,yr,,ss]
        biol@n[1,yr,-(1:ss),ss] <- 0  # The recruitment is 0 in [units > ss].
    }
    else{  # dim(biol@n)[3] = 1, The recruitment only occurs in 1 season every year. 
        if(SR@proportion[,yr,,ss,,1] == 1){ # If the recruitment season is 'ss' generate it otherwise
            SR <- SRsim(SR, year = yr, season = ss, iter = 'all') 
            biol@n[1,yr,,ss] <- SR@rec[,yr,,ss]
       } else if (ss==1) { # If the recruitment season is NOT the first one, set numbers at first age group to 0 in this season.
      biol@n[1,yr,,ss] <- 0
    } # If the recruitment season is NOT the first one, do nothing, the population in first age group is just the survivors of previous season.
    
    
    }
    
    if(any(biol@n[,yr,,ss]<0)){
        biol <- correct.biomass.ASPG(biol, yr, ss)
    }
    
    return(list(biol = biol, SR = SR))

} 


#-------------------------------------------------------------------------------
# BDPG(biol, BD, fleets, biol.control)
# - OUTPUT: list(biol = biol, BD = BD) - Upadated FLBiol and FLSRsim objects.
#-------------------------------------------------------------------------------

BDPG <- function(biols, BDs, fleets, year, season, stknm, ...){

    cat('-----------------BDPG-----------\n')

    if(length(year) > 1 | length(season) > 1)
        stop('Only one year and season is allowed' )
    
    biol <- biols[[stknm]]    
    BD   <- BDs[[stknm]]
    
    dimnms <- dimnames(biol@n)
    
    # If year/season/iter numerics => indicate position 
    # else names => get positions.
    
    if(length(year) > 1 | length(season) > 1)
        stop('Only one year and season is allowed' )
    
    # 'year' dimension.
    yr <- year
    if(is.character(year)) yr <- which(dimnms[[2]] %in% year)
    if(length(yr) == 0) stop('The year is outside object time range')  
    
    # 'season' dimension.
    ss <- season
    if(is.character(season)) ss <- which(dimnms[[4]] %in% season)
    if(length(ss) == 0) stop('The season is outside object season range')  
    
    ns <- dim(biol@n)[4]
    stock <- biol@name
    
    # IF season = 1 
    if(ss == 1){
        yr0 <- yr - 1
        ss0 <- ns

    }
    else{
        yr0 <- yr
        ss0 <- ss - 1 
    }
    
    # total catch in year [yr0] season [ss0].
    BD@catch[,yr0,,ss0] <- catchStock(fleets,stock)[,yr0,,ss0]
        
    # Update FLBDsim object.
    BD <- BDsim(BD,yr,ss)
        
    # Update FLBiol Object
    biol@n[,yr,,ss] <- BD@biomass[,yr,,ss]/biol@wt[,yr,,ss]
    
    return(list(biol = biol, BD = BD))

} 

# correct.biomass.ASPG <- function(biol)
correct.biomass.ASPG <- function(biol, year, season){

    for(i in 1:dim(biol@n)[6]){
        N  <- c(biol@n[,year,, season,,i])
        wt <- c(biol@wt[,year,,season,,i])
    
        # identify the ages that are < 0, this should be only happen in the 
        # first year of simulation, when catch has not been calculated using CobbDoug.
        negs <- which(N<0)
        #biomass proportion in the positive ages, the negative biomass is discounted proportionally 
        # depending on the abundace of each age group. The calculation is done for all cohorts (units) at the same time.
        p <- N[-negs]*wt[-negs]/sum(N[-negs]*wt[-negs])
        # negative biomass
        NegB <- -sum(N[negs]*wt[negs])
        # New Biomass
        Nnew <- N                                       
        Nnew[negs]  <- 0
        Nnew[-negs] <- N[-negs] - (NegB*p/wt[-negs]) 
        
        biol@n[,year,, season,,i] <- matrix(Nnew, dim(biol@n)[1])
        }
        return(biol)
        
}


#-------------------------------------------------------------------------------
# The same as ASPG but the catch is produced using Baranov catch equation, i.e, 
#  the catch is done all along the year using an instantaneous constant 
#  fishing mortality rate. 
#-------------------------------------------------------------------------------

ASPG_Baranov <- function(biols, SRs, fleets, year, season, stknm, ...){
  
  cat('-----------------ASPG-----------\n')
  
  if(length(year) > 1 | length(season) > 1)
    stop('Only one year and season is allowed' )
  
  biol <- biols[[stknm]]     
  SR   <- SRs[[stknm]]
  
  dimnms <- dimnames(biol@n)
  
  # If year/season/iter numerics => indicate position 
  # else names => get positions.
  
  if(length(year) > 1 | length(season) > 1)
    stop('Only one year and season is allowed' )
  
  # 'year' dimension.
  yr <- year
  if(is.character(year)) yr <- which(dimnms[[2]] %in% year)
  if(length(yr) == 0) stop('The year is outside object time range')  
  
  # 'season' dimension.
  ss <- season
  if(is.character(season)) ss <- which(dimnms[[4]] %in% season)
  if(length(ss) == 0) stop('The season is outside object season range')  
  
  na <- dim(biol@n)[1]
  ns <- dim(biol@n)[4]
  stock <- biol@name
  
  findF <- function(E, qa, Ca, Ma, Na){
    Ca. <- (Fa/(Fa+Ma))*(1-exp(-Ma-Fa))*Na
    res <- sum((Ca.-Ca)^2)
    return(res)
  }
  
  # IF season = 1 THEN age groups move to the next. 
  if(ss == 1){
    # total catch in year [y-1] season [ns].
    catch.n <- catchStock(fleets,stock)[,yr-1,,ns]
    
    
    # middle ages
    biol@n[-c(1,na),yr,,ss] <- (biol@n[-c(na-1,na),yr-1,,ns]*exp(-biol@m[-c(na-1,na),yr-1,,ns]/2) - catch.n[-c(na-1,na),])*
      exp(-biol@m[-c(na-1,na),yr-1,,ns]/2) 
    # plusgroup
    biol@n[na,yr,,ss]       <- (biol@n[na-1,yr-1,,ns]*exp(-biol@m[na-1,yr-1,,ns]/2) - catch.n[na-1,])*exp(-biol@m[na-1,yr-1,,ns]/2) + 
      (biol@n[na,yr-1,,ns]*exp(-biol@m[na,yr-1,,ns]/2) - catch.n[na,])*exp(-biol@m[na,yr-1,,ns]/2)
    
  }
  else{
    # total catch in year [yr] season [ss-1].
    catch.n <- catchStock(fleets,stock)[,yr,,ss-1]
    # middle ages      # for unit == ss  and age = 1, it will be equal NA but be updated after with SRsim.
    biol@n[,yr,,ss] <- (biol@n[,yr,,ss-1]*exp(-biol@m[,yr,,ss-1]/2) - catch.n)*exp(-biol@m[,yr,,ss-1]/2) 
  }
  
  # Update SSB.
  SR@ssb[,yr,,ss] <- unitSums(quantSums(n(biol) * wt(biol) * fec(biol)*mat(biol) * 
                                          exp(-biol@m*spwn(biol)), na.rm=TRUE))[,yr,,ss]
  
  # RECRUITMENT
  if(dim(biol@n)[3] > 1 & dim(biol@n)[3] == dim(biol@n)[4]){
    # 'number_of_units > 1' => Recruitment occurs in all seasons, the 'unit' correspond with the recruitment 'season'.
    SR <- SRsim(SR, year = yr, season = ss, iter = 'all') 
    biol@n[1,yr,ss,ss] <- SR@rec[,yr,,ss]
    biol@n[1,yr,-(1:ss),ss] <- 0  # The recruitment is 0 in [units > ss].
  }
  else{  # dim(biol@n)[3] = 1, The recruitment only occurs in 1 season every year. 
    if(SR@proportion[,yr,,ss,,1] == 1){ # If the recruitment season is 'ss' generate it otherwise
      SR <- SRsim(SR, year = yr, season = ss, iter = 'all') 
      biol@n[1,yr,,ss] <- SR@rec[,yr,,ss]
    } else if (ss==1) { # If the recruitment season is NOT the first one, set numbers at first age group to 0 in this season.
      biol@n[1,yr,,ss] <- 0
    } # If the recruitment season is NOT the first one, do nothing, the population in first age group is just the survivors of previous season.
    
    
  }
  
  if(any(biol@n[,yr,,ss]<0)){
    biol <- correct.biomass.ASPG(biol, yr, ss)
  }
  
  return(list(biol = biol, SR = SR))
  
}

    
#-------------------------------------------------------------------------------
# gadgetGrowthReal(biol, GDGT, fleets, biol.control)
# - OUTPUT: list(biol = biol, SR = SR, GDGT = GDGT) - Upadated FLBiol, SR and GDGT objects.
#-------------------------------------------------------------------------------

gadgetGrowth <- function(biols, GDGTs, SRs, fleets, year, season, stknm, ...){

	biol <- biols[[stknm]]
	SR   <- SRs[[stknm]]

	GDGT <- GDGTs

        if(is.null(GDGT$runNow) || GDGT$runNow == FALSE)
		return(list(biol = biol, SR = SR, GDGT = GDGT))

	print(paste("Runing Gadget", season, year, GDGT$gadget.inputDir))

	# Name converter
	curGadgetStockName <- convertStockName[[stknm]]

	# Load params
	stockParams <- eval(parse(text=paste0(curGadgetStockName, ".params")))

	# Check if this is gadget first run?
	if(GDGT$firstRun){

		startYear <- GDGT$startYear
		curYear <- startYear

		print(paste("Start year is", curYear, (curYear - startYear + 1)))

		# Run until the current FLBEIA year
		gadgetSimOut <- runUntil(projYear-1)

		# Expand FLStock and FLIndex
		gadgetSimOut[[curGadgetStockName]]$stk <- expand(gadgetSimOut[[curGadgetStockName]]$stk, year=firstYear:finalYear)
		gadgetSimOut[[curGadgetStockName]]$idx <- expand(gadgetSimOut[[curGadgetStockName]]$idx, year=firstYear:finalYear)

		# Collect hindcast stats
		SR@ssb <- ssb(gadgetSimOut[[curGadgetStockName]]$stk)
		SR@rec <- rec(gadgetSimOut[[curGadgetStockName]]$stk)
		biol@wt<- stock.wt(gadgetSimOut[[curGadgetStockName]]$stk)
		biol@n <- stock.n(gadgetSimOut[[curGadgetStockName]]$stk)

		# Run for this year and collect the stats
		simInfo <- getEcosystemInfo()
		curYear <- simInfo[["time"]][["currentYear"]]
		stats <- runYear()

		# Create the list for the stats
		GDGT[["currentStats"]] <- list()

		# Put the statistics information for this year
		GDGT[["currentStats"]][[as.character(year)]] <- stats

		# Put the FLstock and FLindex information
		GDGT[["gadgetSimOut"]] <- gadgetSimOut

		GDGT$firstRun <- FALSE
	}else{
		# Check whether this is another species or a start of the year
		simInfo <- getEcosystemInfo()
		curYear <- simInfo[["time"]][["currentYear"]]
		startYear <- GDGT$startYear

		print(paste("Year now is", curYear, (curYear - startYear + 1)))

		if((curYear - startYear + 1) == year){
			# Run subsequent steps
			stats <- runYear()
			# Put the statistics information for this year
			GDGT[["currentStats"]][[as.character(year)]] <- stats
			# Start from the first stock
		}else{
			# Get stats for the current year
			stats <- GDGT[["currentStats"]][[as.character(year)]]
			curYear <- curYear - 1
		}
	}

	iter <- 1

	# Update FLStock(and FLIdx) (using the latest simout)
	gadgetSimOut <- GDGT[["gadgetSimOut"]]
	gadgetSimOut[[curGadgetStockName]] <- updateFLStock(curGadgetStockName, stats, as.character(curYear), gadgetSimOut[[curGadgetStockName]]$stk, gadgetSimOut[[curGadgetStockName]]$idx, stockParams[["stockStep"]])

	# Put the updated FLstock and FLindex information
	GDGT[["gadgetSimOut"]] <- gadgetSimOut

	SR@ssb[, as.character(curYear)]  <- ssb(gadgetSimOut[[curGadgetStockName]]$stk)[, as.character(curYear)]
	SR@rec[, as.character(curYear)]  <- rec(gadgetSimOut[[curGadgetStockName]]$stk)[, as.character(curYear)]
	biol@wt[, as.character(curYear)] <- stock.wt(gadgetSimOut[[curGadgetStockName]]$stk)[, as.character(curYear)]
	biol@n[, as.character(curYear)] <- stock.n(gadgetSimOut[[curGadgetStockName]]$stk)[, as.character(curYear)]

	# If reaches the end of gadget simulation
	print("Current gadget information:\n")
	simInfo <- getEcosystemInfo()
	print(simInfo)

	if(simInfo[["time"]][["finished"]] == 1){
		# Sim cleanup
		finalizeSim()

		# Get the output
		out <- finalize()
	}

	print("Biol@wt")
	print(biol@wt)
	print("Biol@n")
	print(biol@n)

	print("SR@ssb")
	print(SR@ssb)
	print("SR@rec")
	print(SR@rec)

	print("Gadget ends")

	return(list(biol = biol, SR = SR, GDGT = GDGT))
}
