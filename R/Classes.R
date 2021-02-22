#' Class Clustering
#'
#' Partition data into clusters
#'
#' @name Clustering-class
#' @rdname Clustering-class
#' @docType class
#'
#' @slot Cluster List of matrices, where each element of the list include data items belonging to the corresponding cluster.
#' @slot Centers Matrix of cluster centers.
#' @slot k Number of centers.
#' @slot clusterMem Numeric vector of cluster membership.
#'
#' @import flexclust
#'
#' @details Object of class \code{Clustering} containing clustering solution of cumulative dataset D_i. The object of class
#' \code{Clustering} comprise of four slots. Slot \code{Clusters} contain data items of each cluster, slot \code{Centers} contain
#' cluster centers, slot \code{k} contain the number of centers, while slot \code{clusterMem} contain cluster memberships vector.
#'
Clustering <- setClass("Clustering",
                       slots = c(
                           Cluster = "list",
                           k = "numeric",
                           Centers = "matrix",
                           clusterMem = "numeric"
                       )
)

#==============================================================================#
#' Overlap between clusters
#'
#' @description Contains matrix of similarity indices between clusters, after clustering dynamic datasets at consecutive time points.
#'
#' @name OverLap-class
#' @rdname OverLap-class
#' @docType class
#'
#' @slot Overlap A numeric matrix containing the similarity index between clusters extracted at time point \code{t_1} and \code{t_2}.
#' The rows of the matrix illustrate clusters extracted from first clustering \eqn{\xi_1(time point t_1)},whereas columns represent
#' clusters extracted from second clustering \eqn{\xi_2(time point t_2)}.
#' @slot rx A numeric vector containg radius of each cluster from first clustering \eqn{\xi_1}.
#' @slot ry A numeric vector containg radius of each cluster from second clustering \eqn{\xi_2}.
#' @slot Centersx A numeric vector containing centers of clusters from first clustering \eqn{\xi_1}.
#' @slot Centersy A numeric vector containing centers of clusters from second clustering \eqn{\xi_2}.
#' @slot avgDisx A numeric vector containing average distance between points in a cluster from its center in first clustering \eqn{\xi_1}.
#' @slot avgDisy A numeric vector containing average distance between points in a cluster from its center in second clustering \eqn{\xi_2}.
#' @slot clusterMem A vector of integers containing cluster membership from second clustering \eqn{\xi_2}.
#'
#' @exportClass OverLap
#' @export OverLap
#'
OverLap <- setClass("OverLap",
                    slots = c(
                        Overlap = "matrix",
                        Centersx = "matrix",
                        Centersy = "matrix",
                        clusterMem = "numeric",
                        avgDisx = "numeric",
                        avgDisy = "numeric",
                        rx = "numeric",
                        ry = "numeric"
                    )
)

#==============================================================================#
#' External Transition Count
#'
#' @description Trace cluster solutions of dynamic datasets at consecutive
#' time points and counts the clusters that experiences external transition.
#' External transition includes Survive, Split, Merge, newly emerged, and Died candidates.
#'
#' @name TransitionCount-class
#' @rdname TransitionCount-class
#' @docType class
#'
#' @slot Survive Number of candidates survive from first clustering \eqn{\xi_1}.
#' @slot Split Number of candidates from first clustering \eqn{\xi_1} that split into several daughter clusters at second clustering \eqn{\xi_2}.
#' @slot Merge Number of candidates from first clustering \eqn{\xi_1} that merge toghter at second clustering \eqn{\xi_2}.
#' @slot Died Number of candidates from first clusterin \eqn{\xi_1} that disapeared at second clustering \eqn{\xi_2}.
#' @slot SurvivalRatio Ratio of survive clusters to total number of clusters from first clusering \eqn{\xi_1}.
#' @slot AbsorptionRatio Ratio of Merged clusters to total number of clusters from first clusering \eqn{\xi_1}.
#' @slot passforwardRatio Sum of SurvivalRatio and AbsorptionRatio.
#' @slot Survival_thrHold Threshold for survival of clusters.
#' @slot Split_thrHold Threhold for split of clusters.
#' @slot Cluster_Tracex Vector containing each cluster result from first clustering \eqn{\xi_1}.
#'
TransitionCount <- setClass("TransitionCount",
                            contains = "OverLap",
                            slots = c(
                                Survive = "numeric",
                                Split = "numeric",
                                Merge = "numeric",
                                Died = "numeric",
                                SurvivalRatio = "numeric",
                                AbsorptionRatio = "numeric",
                                passforwardRatio = "numeric",
                                Survival_thrHold = "numeric",
                                Split_thrHold = "numeric",
                                Cluster_Tracex = "vector"
                            ),

                            prototype = prototype(
                                Survive = 0,
                                Split = 0,
                                Died = 0,
                                Merge = 0,
                                SurvivalRatio = 0,
                                AbsorptionRatio = 0,
                                passforwardRatio = 0
                            ),

                            # Function for testing data consistency
                            validity = function(object)
                            {
                                if(object@Survival_thrHold < 0.0 | object@Survival_thrHold > 1.0)
                                    stop("Error: The survival threshold is out of bounds")
                                else if(object@Split_thrHold < 0.0 | object@Split_thrHold > 1.0)
                                    stop("Error: The Split threshold is out of bounds")
                                else
                                    return(TRUE)
                            }
)

#==============================================================================#
#' External Transition Candidates
#'
#' @description Class containing candidates that adopted external transition from first clustering \eqn{\xi_1},
#' and emerged as new clusters at second clustering \eqn{\xi_2}.
#'
#' @name TransitionCan-class
#' @rdname TransitionCan-class
#' @docType class
#'
#' @slot SurvivalCanx Vector of integers comprising Candidates that Survive from first clustering \eqn{\xi_1}.
#' @slot SurvivalCany Vector of integers comprising Candidates that Survive to second clustering \eqn{\xi_2}.
#' @slot SplitCanx Vector of integers comprising Candidates that Sliced into Various daughter Clusters from
#' first clustering \eqn{\xi_1}.
#' @slot SplitCany List of integer vectors comprising Candidates that emerged as daughter clusters in second clustering \eqn{\xi_2}
#' because of Split from first clustering \eqn{\xi_1}.
#' @slot MergeCanx List of integer vectors comprising Candidates from first clustering \eqn{\xi_1} that are merged.
#' Each slot of list indicates the clusters that merge together from first clustering.
#' @slot MergeCany Vector of integers comprising Candidates that emerged in second clustering \eqn{\xi_2} because
#' of merging various clusters from first clustering \eqn{\xi_1}.
#' @slot EmergCan Newley emerged candidates which are not a result of any external transition from first clustering \eqn{\xi_1}.
#' @slot Cluster_Tracey Vector of Cluster Trace from second clustering \eqn{\xi_2}.
#'
TransitionCan <- setClass("TransitionCan",
                          contains = "TransitionCount",
                          slots = c(
                              SurvivalCanx = "vector",
                              SurvivalCany = "vector",
                              SplitCanx = "vector",
                              SplitCany = "vector",
                              MergeCanx = "vector",
                              MergeCany = "vector",
                              EmergCan = "vector",
                              Cluster_Tracey = "vector"
                          ),

                          prototype = list(
                              SplitCany = vector(mode = "list", length = 10),
                              MergeCanx = vector(mode = "list", length = 20)
                          )
)

#==============================================================================#
#' Internal Transition Candidates
#'
#' @description Class containing results of Internal Transition of survived clusters from first clustering \eqn{\xi_1}.
#'
#' @name intTransitionCan-class
#' @rdname intTransitionCan-class
#' @param object An object of class Transitioncan
#'
#' @slot Location.diff Vector of integers containing difference in location (= Distance bw cluster centers/min(rx,ry)).
#' @slot Compactness.diff Vector of integers containing Change in density of survived clusters (d(rx, ry)).
#' @slot Location_thrHold Minimum value of threshold for shift in location.
#' @slot Density_thrHold Minimum value of threshold for change in density.
#' @slot ShiftLocCan Vector of integers containing Survived candidates with shift in their location.
#' @slot NoShiftLocCan Vector of integers containing Survived candidates with no Shift in their Location.
#' @slot MoreCompactCan Vector of integers containing Survived Candidates Which becomes more compact.
#' @slot MoreDiffuseCan Vector of integers containing Survived Candidates Which becomes more diffuse.
#' @slot NoChangeCompactCan Vector of integers containing Survived candidates with no change in compactness.
#'
intTransitionCan <- setClass("intTransitionCan",
                             contains = "TransitionCan",
                             slots = c(
                                 Location.diff = "numeric",
                                 Compactness.diff = "numeric",
                                 ShiftLocCan = "numeric",
                                 NoShiftLocCan = "numeric",
                                 MoreCompactCan = "numeric",
                                 MoreDiffuseCan = "numeric",
                                 NoChangeCompactCan = "numeric",
                                 Location_thrHold = "numeric",
                                 Density_thrHold = "numeric"
                             ),

                             prototype = prototype(
                                 Location.diff = 0,
                                 Compactness.diff = 0,
                                 ShiftLocCan = numeric(0),
                                 NoShiftLocCan = numeric(0),
                                 MoreCompactCan = numeric(0),
                                 MoreDiffuseCan = numeric(0),
                                 NoChangeCompactCan = numeric(0)
                             ),

                             # Function for testing data consistency
                             validity = function(object)
                             {
                                 if(object@Location_thrHold < 0.0 | object@Location_thrHold > 1.0)
                                     stop("Error: The location threshold is out of bounds")
                                 else if(object@Density_thrHold < 0.0 | object@Density_thrHold > 1.0)
                                     stop("Error: The density threshold is out of bounds")
                                 else
                                     return(TRUE)
                             }
)

#==============================================================================#
#' An S4 class that contain time steps
#'
#' @name Monic
#' @rdname Monic-class
#' @param object An object of class \code{Transitioncan}
#'
#' @slot TimeStep Time Steps
#'
setClass("Monic",
         slots = c(
             TimeStep = "list"
         )
)

