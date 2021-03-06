#-------------------------------------------------------------------------------
#                    ** create.biols.ctrl **
#
# The existing functions to simulate population growth:
#  (ASPG, BDPG and 'fixedPopulation')
# do not have extra arguments, thus we don't call additional create.XYZ.ctrl  
# within create.biols.ctrl
#
#' biols.ctrl object creator
#' 
#' It creates the biols.ctrl object to be used in the call to the main function FLBEIA.
#' 
#   :: ARGUMENTS ::
#
#' @param stksnames A vector with the name of the stocks in the OM.
#' @param growth.models A character vector of the same length as stksnames with the name of the model used to project the stock populations in the simulation.
#' @param immediate logical, indicating if the warnings should be output immediately.
#' @param ... any extra arguments necessary in the model specific creators. '...' are extracted using 'list(...)', this generates a named list with the extra arguments.
#'        To assure the correct functioning the extra arguments must have a name.
#' 
#' @return A list of lists with the basic structure of the biols.ctrl object.

# Dorleta Garc?a - Azti Tecnalia
# 28/05/2013 10:42:07
#-------------------------------------------------------------------------------

create.biols.ctrl <- function(stksnames, growth.models = NULL, immediate = FALSE,...){
    
    growth.models.available <- c('fixedPopulation', 'ASPG', 'BDPG', 'gadgetGrowth')
    
    nstk <- length(stksnames) 
    res  <- vector('list', nstk)
    names(res) <- stksnames
    
    extra.args <- list(...)
    
    if(is.null(growth.models)) growth.models <- rep('fixedPopulation',nstk)
    else{ 
      if(length(growth.models) != nstk) stop("'growth.models' must be NULL or must have the same length as stknames'")
      if(!all(growth.models %in% growth.models.available)){ 
        wmod <- unique(growth.models[which(!(growth.models %in% growth.models.available))])  
        warning(paste(unique(wmod), collapse = "-")," in 'growth.models' is not an internal FLBEIA growth model. If you want to use create.biols.ctrl you must create, ", paste('create', paste(unique(wmod), collapse = ', ') ,'ctrl', sep = ".")," function.", immediate. = immediate)
    }}
   
    for(stk in 1:nstk){
        res[[stk]]                   <- vector('list',1)
        names(res[[stk]])            <- 'growth.model'
        res[[stk]][['growth.model']] <- growth.models[stk]
    } 
    
    for(st in 1:nstk){
        
        growthmodcreator <- paste('create', growth.models[st],  'ctrl', sep = '.')
        res[[st]] <- eval(call(growthmodcreator, res = res[[st]], stkname = st, largs = extra.args))
    }
    
    return(res) 
} 


#-------------------------------------------------------------------------------
#                       ** create.fixedPopulation.ctrl **
# No extra arguments needed
#-------------------------------------------------------------------------------
create.fixedPopulation.ctrl <- function(resst,stkname,largs) return(resst)

#-------------------------------------------------------------------------------
#                       ** create.ASPG.ctrl **
# No extra arguments needed
#-------------------------------------------------------------------------------
create.ASPG.ctrl <- function(resst,stkname,largs) return(resst)

#-------------------------------------------------------------------------------
#                       ** create.BDPG.ctrl **
# No extra arguments needed
#-------------------------------------------------------------------------------
create.BDPG.ctrl <- function(resst,stkname,largs) return(resst)

#-------------------------------------------------------------------------------
#                       ** create.gadgetGrowth.ctrl **
# No extra arguments needed
#-------------------------------------------------------------------------------
create.gadgetGrowth.ctrl <- function(resst,stkname,largs) return(resst)
