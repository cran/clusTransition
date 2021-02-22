#==============================================================================#
#' Clustering.
#'
#' @name Clusters
#' @rdname Clusters-methods
#'
setGeneric("Clusters",
           def = function(object, x, k)
           {
               standardGeneric("Clusters")
           }
)

#==============================================================================#
#' Overlap
#'
#' @name Overlap
#' @rdname Overlap-methods
#'
#' @exportMethod Overlap
#'
setGeneric("Overlap",
           def = function(object, e1, e2)
           {
               standardGeneric("Overlap")
           }
)

#==============================================================================#
#' External Transition Count
#'
#' @name extTransitionCount
#' @rdname extTransitionCount-methods
#'
setGeneric("extTransitionCount",
           def = function(object)
           {
               standardGeneric("extTransitionCount")
           }
)

#==============================================================================#
#' External Transition Candidate.
#'
#' @name extTransitionCan
#' @rdname extTransitionCan-methods
#'
setGeneric("extTransitionCan",
           def = function(object)
           {
               standardGeneric("extTransitionCan")
           }
)

#==============================================================================#
#' Internal Transition Candidates.
#'
#' @description This method identify internal transition of the survived clusters,
#' obtained from 'extTransitionCan()' method.
#'
#' @name internalTransition
#' @rdname internalTransition-methods
#'
setGeneric("internalTransition",
           def = function(object)
           {
               standardGeneric("internalTransition")
           }
)
