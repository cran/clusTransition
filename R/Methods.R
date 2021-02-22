#==============================================================================#
#' Partition data into clusters
#'
#' @description Initialize slots of class \code{Clustering} by partitioning the dataset into \code{k} clusters.
#'
#' @rdname Clusters-methods
#' @aliases Clusters,Clustering-method
#' @docType methods
#'
#' @param object An object of class \code{Clustering}.
#' @param x Numeric matrix of data.
#' @param k Number of centers.
#'
#' @details Runs \code{cclust} function from "flexclust" package with default settings i.e.
#' method = "kmeans", dist = "euclidean", and partition the dataset. Returns object of class \code{Clustering}.
#'
#' @import flexclust
#'
#' @return An object of class \code{Clustering}
#'
setMethod("Clusters",
          signature(object = "Clustering", x = "matrix", k = "numeric"),
          definition = function(object, x, k)
          {
              fit <- cclust(x, k, dist = "euclidean", method = "kmeans")
              for(i in 1:k)
                  object@Cluster[[i]] <- x[fit@cluster==i,]
              object@Centers <- fit@centers
              object@clusterMem <- fit@cluster
              object@k <- k
              return(object)
          }
)

#==============================================================================#
#' Similarity index between clusters
#'
#' @description Initialize slots of class \code{OverLap} by importing clustering solutions of dynamic
#' datasets at two consecutive time points. Clusters at each time point should be provided as a list of
#' matrices, where each matrix contains dataset belongs to the corresponding cluster.
#'
#' @rdname Overlap-methods
#' @aliases OverLap, Clustering, Clustering-method
#' @docType methods
#'
#' @param object An object of class \code{OverLap}
#' @param e1 An object of class \code{Clustering}, or any object that can be coerced, such as list of matrices
#' or data frames that contain clusters from first clustering.
#' @param e2 An object of class \code{Clustering}, or any object that can be coerced, such as list of matrices
#' or data frames that contain clusters from second clustering.
#' @return Return an object of class \code{OverLap}.
#'
setMethod("Overlap",
          signature(object = "OverLap", e1 = "Clustering", e2 = "Clustering"),
          definition = function(object, e1, e2)
          {
              object@Overlap <- matrix(object@Overlap, nrow = e1@k, ncol = e2@k)
              object@Centersx <- e1@Centers
              object@Centersy <- e2@Centers
              object@clusterMem <- e2@clusterMem

              for(i in 1:e1@k)
              {
                  for(j in 1:e2@k)
                      object@Overlap[i,j] <- nrow(
                          merge(e1@Cluster[[i]], e2@Cluster[[j]], all.x = F, all.y = F))/nrow(e1@Cluster[[i]])
              }

              for(i in 1:e1@k)
              {
                  object@rx[i] <- max(apply(e1@Cluster[[i]], 1, function(x)sqrt(sum(x - e1@Centers[i,])^2)))
                  object@avgDisx[i] <- sum(apply(e1@Cluster[[i]], 1, function(x)sqrt(sum(x - e1@Centers[i,])^2)))/length(e1@Cluster[[i]])
              }

              for(i in 1:e2@k)
              {
                  object@ry[i] <- max(apply(e2@Cluster[[i]], 1, function(x)sqrt(sum(x - e2@Centers[i,])^2)))
                  object@avgDisy[i] <- sum(apply(e2@Cluster[[i]], 1, function(x)sqrt(sum(x - e2@Centers[i,])^2)))/length(e2@Cluster[[i]])
              }

              rownames(object@Overlap) <- paste("Cluster", 1:e1@k, sep = ".")
              colnames(object@Overlap) <- paste("Cluster", 1:e2@k, sep = ".")
              return(object)
          }
)

#==============================================================================#
#' @rdname Overlap-methods
#' @aliases OverLap,ANY,ANY-methods
#' @docType methods
#'
#' @exportMethod Overlap
#' @export Overlap
#'
setMethod("Overlap",
          signature(object = "OverLap", e1 = "ANY", e2 = "ANY"),
          definition = function(object, e1, e2)
          {
              object@Overlap <- matrix(object@Overlap, nrow = length(e1), ncol = length(e2))
              object@Centersx <- matrix(object@Centersx, nrow = length(e1), ncol = ncol(e1[[1]]))
              object@Centersy <- matrix(object@Centersy, nrow = length(e2), ncol = ncol(e2[[1]]))

              for(i in 1:length(e1))
              {
                  for(j in 1:length(e2))
                      object@Overlap[i,j] <- (nrow(
                          merge(e1[[i]], e2[[j]], all.x = F, all.y = F)))/nrow(e1[[i]])
              }

              for(i in 1:length(e1))
              {
                  object@Centersx[i,] <- colMeans(e1[[i]])
                  object@rx[i] <- max(apply(e1[[i]], 1, function(x)sqrt(sum(x - object@Centersx[i,])^2)))
                  object@avgDisx[i] <- sum(apply(e1[[i]], 1, function(x)sqrt(sum(x - object@Centersx[i,])^2)))/length(e1[[i]])
              }

              for(i in 1:length(e2))
              {
                  object@Centersy[i,] <- colMeans(e2[[i]])
                  object@ry[i] <- max(apply(e2[[i]], 1, function(x)sqrt(sum(x - object@Centersy[i,])^2)))
                  object@avgDisy[i] <- sum(apply(e2[[i]], 1, function(x)sqrt(sum(x - object@Centersy[i,])^2)))/length(e2[[i]])
              }
              colnames(object@Centersx) <- colnames(e1[[1]])
              colnames(object@Centersy) <- colnames(e2[[1]])
              rownames(object@Overlap) <- paste("Cluster", 1:length(e1), sep = ".")
              colnames(object@Overlap) <- paste("Cluster", 1:length(e2), sep = ".")
              return(object)
          }
)

#==============================================================================#
#' External Transiton Counts
#'
#' @description Trace cluster solutions of dynamic datasets and count the number of clusters that experiences
#' external transition from first clustering. The external transition includes survived, split into various daughters,
#' spliced into one, disappeared, and newly emerged candidates.
#'
#' @rdname extTransitionCount-methods
#' @aliases extTransitionCount-method
#' @docType methods
#'
#' @param object An object of class \code{Transitioncount}
#' @return Return an object of class \code{TransitionCount}
#'
setMethod("extTransitionCount",
          signature(object ="TransitionCount"),
          definition = function(object)
          {
              R <- vector()
              SpO <- vector()
              AbC <- vector()

              R <- max.col(object@Overlap)
              p <- apply(object@Overlap, 1, max)
              l <- apply(object@Overlap, 2, max)
              AbC <- apply(object@Overlap, 2,
                           function(c)sum(c >=  object@Survival_thrHold))
              SpO <- apply(object@Overlap, 1,
                           function(x)sum(x[(x <= object@Survival_thrHold) &
                                                (x >= object@Split_thrHold) & (l<object@Survival_thrHold)]))


              for(i in 1:nrow(object@Overlap))
              {
                  #----- Count for Survival
                  if((p[i] >= object@Survival_thrHold) & (AbC[R[i]] == 1))
                  {
                      object@Survive <- object@Survive + 1
                      object@Cluster_Tracex[i] <- "Survive"
                  }
                  #----- Count for Absorptions
                  else if((p[i] >= object@Survival_thrHold) & (AbC[R[i]] > 1))
                  {
                      object@Merge <- object@Merge + 1
                      object@Cluster_Tracex[i] <- "Merged"
                  }
                  #----- Count for Split
                  else if((p[i] < object@Survival_thrHold) & (SpO[i] >= object@Survival_thrHold))
                  {
                      object@Split <- object@Split + 1
                      object@Cluster_Tracex[i] <- "Splited"
                  }
                  #---- Count for Disapearence
                  else
                  {
                      object@Died <- object@Died + 1
                      object@Cluster_Tracex[i] <- "Died"
                  }
              }
              object@SurvivalRatio = object@Survive / nrow(object@Overlap)
              object@AbsorptionRatio = object@Merge / nrow(object@Overlap)
              object@passforwardRatio = object@SurvivalRatio + object@AbsorptionRatio
              return(object)
          }
)

#==============================================================================#
#' External Transiton Candidates
#'
#' @description This S4 method trace cluster solutions of dynamic dataset, and identify the candidates
#' that experience external transition from first clustering and emerged at second clustering.
#'
#' @rdname extTransitionCan-methods
#' @aliases extTransitionCan-method
#' @docType methods
#'
#' @param object An object of class Transitioncan
#' @return Return an object of class TransitionCan
#'
setMethod("extTransitionCan",
          signature(object ="TransitionCan"),
          definition = function(object)
          {

              length(object@Cluster_Tracey) <- ncol(object@Overlap)

              R <- vector()
              SpO <- vector()
              AbC <- vector()

              R <- max.col(object@Overlap)
              p <- apply(object@Overlap, 1, max)
              l <- apply(object@Overlap, 2, max)
              AbC <- apply(object@Overlap, 2,
                           function(c)sum(c >=  object@Survival_thrHold))
              SpO <- apply(object@Overlap, 1,
                           function(x)sum(x[(x <= object@Survival_thrHold) &
                                                (x >= object@Split_thrHold) & (l<object@Survival_thrHold)]))

              for(i in 1:nrow(object@Overlap))
              {
                  #---- Survived Candidates threshold
                  if((p[i] >= object@Survival_thrHold) & (AbC[R[i]] <= 1))
                  {
                      object@SurvivalCanx <- c(object@SurvivalCanx, i)
                      object@SurvivalCany <- c(object@SurvivalCany, R[i])
                      object@Cluster_Tracey[R[i]] <- "Survived"
                  }
                  #---- Split Candidates threshold
                  else if((p[i] < object@Survival_thrHold) & (SpO[i] >= object@Survival_thrHold))
                  {
                      object@SplitCanx <- c(object@SplitCanx, i)
                      for(j in 1:ncol(object@Overlap))
                      {
                          if(object@Overlap[i,j]>=object@Split_thrHold)
                          {
                              object@SplitCany[[i]] <- c(object@SplitCany[[i]], j)
                              object@Cluster_Tracey[j] <- "Splited"
                          }
                      }
                  }
              }
              object@SplitCany <- object@SplitCany[-which(sapply(object@SplitCany, is.null))]

              #---- Absorb Candidates threshold
              for(i in 1:ncol(object@Overlap))
              {
                  if(AbC[i] > 1)
                  {
                      object@MergeCany <- c(object@MergeCany, i)
                      for(j in 1:nrow(object@Overlap))
                      {
                          if(object@Overlap[j,i] >= object@Survival_thrHold)
                          {
                              object@MergeCanx[[i]] <- c(object@MergeCanx[[i]], j)
                              object@Cluster_Tracey[i] <- "Absorbed"
                          }
                      }
                  }
              }
              object@MergeCanx <- object@MergeCanx[-which(sapply(object@MergeCanx, is.null))]

              #---- Emerging Candidates
              object@Cluster_Tracey[which(is.na(object@Cluster_Tracey))] <- "Emerged"
              object@EmergCan <- as.numeric(which(object@Cluster_Tracey == "Emerged"))
              return(object)
          }
)

#==============================================================================#
#' Internal Transition Candidates
#'
#' @description Trace clustering solutions of cumulative datasets and identify the survived clusters
#' experiencing Internal transitions. Internal transition includes the change in location and density
#' of the survived candidates.
#'
#' @rdname internalTransition-methods
#' @aliases internalTransition-method
#' @docType methods
#'
#' @param object An object of class \code{intTransitionCan}
#' @return Return an object of class \code{intTransitionCan}
#'
setMethod("internalTransition",
          signature(object ="intTransitionCan"),
          definition = function(object)
          {
              if(object@Survive > 0)
              {
                  for(i in 1:object@Survive)
                  {
                      distance <- sqrt(sum(
                          (object@Centersx[object@SurvivalCanx[i],] - object@Centersy[object@SurvivalCany[i],])^2))

                      object@Location.diff[i] <- distance/min(object@rx[object@SurvivalCanx[i]], object@ry[object@SurvivalCany[i]])
                      object@Compactness.diff[i] <- object@avgDisx[object@SurvivalCanx[i]]- object@avgDisy[object@SurvivalCany[i]]

                      #---- Check for location transition
                      if(object@Location.diff[i] >= object@Location_thrHold)
                          object@ShiftLocCan <- c(object@ShiftLocCan, object@SurvivalCanx[i])
                      else
                          object@NoShiftLocCan <- c(object@NoShiftLocCan, object@SurvivalCanx[i])

                      #---- Check for compactness transition
                      if(abs(object@Compactness.diff[i]) >= object@Density_thrHold)
                      {
                          if(object@Compactness.diff[i] > 0)
                              object@MoreCompactCan <- c(object@MoreCompactCan, object@SurvivalCanx[i])
                          else
                              object@MoreDiffuseCan <- c(object@MoreDiffuseCan, object@SurvivalCanx[i])
                      }
                      else
                          object@NoChangeCompactCan <- c(object@NoChangeCompactCan, object@SurvivalCanx[i])
                  }
              }
              return(object)
          }
)

#==============================================================================#
#' Show Method for output
#'
#' @rdname show-methods
#' @aliases show-method
#' @docType methods
#'
#' @param object An object of class Monic
#'
setMethod("show",
          "Monic",
          definition = function(object)
          {
              cat("\n Survival Threshold = ", object@TimeStep[[1]]@Survival_thrHold,
                  ", Split Threshold = ", object@TimeStep[[1]]@Split_thrHold)
              cat("\n\n External Cluster Transition Count \n\n")
              Table1 <-data.frame(sapply(object@TimeStep,attr, 'Survive'),
                                  sapply(object@TimeStep,attr, 'Merge'),
                                  sapply(object@TimeStep,attr, 'Split'),
                                  sapply(object@TimeStep,attr, 'Died'),
                                  lengths(lapply(object@TimeStep,attr, 'EmergCan')))

              Table2 <-data.frame(lengths(lapply(object@TimeStep,attr, 'ShiftLocCan')),
                                  lengths(lapply(object@TimeStep,attr, 'NoShiftLocCan')),
                                  lengths(lapply(object@TimeStep,attr, 'MoreCompactCan')),
                                  lengths(lapply(object@TimeStep,attr, 'MoreDiffuseCan')),
                                  lengths(lapply(object@TimeStep,attr, 'NoChangeCompactCan')))

              colnames(Table1)<-c("Survive", "Merged", "Split", "Died", "New.Emerged")
              row.names(Table1)<-paste("Time Step", 1:(length(object@TimeStep)), sep = ".")
              print(Table1)

              cat("\n Location Threshold = ", object@TimeStep[[1]]@Location_thrHold,
                  ", Density Threshold = ", object@TimeStep[[1]]@Density_thrHold)
              cat("\n\n Internal Cluster Transition Count \n\n")

              colnames(Table2)<-c("Shif", "No.Shift", "Compact", "Diffuse", "No.change.Compactness")
              row.names(Table2)<-paste("Time Step", 1:(length(object@TimeStep)), sep = ".")
              print(Table2)
          }
)

#==============================================================================#
#' plot Method for output
#'
#' @name moplot
#' @rdname moplot-methods
#' @exportMethod moplot
#'
setGeneric("moplot",
           def = function(object)
           {
               standardGeneric("moplot")
           }
)

#==============================================================================#
#' plot.Monic
#'
#' @description This method plot 3 barplot and 1 line graph. The first stack barplot shows
#' SurvivalRatio and AbsorptionRatio, second barplot shows number of newly emerged clusters
#' at each time stamp, third barplot shows number of disapeared clusters at each time stamp. The line
#' graph shows passforward Ratio and Survival Ratio.
#'
#' @rdname moplot-methods
#' @aliases moplot-method
#' @docType methods
#'
#' @param object An object of class Monic
#'
#' @exportMethod moplot
#'
#' @importFrom graphics par legend lines
#'
setMethod("moplot",
          signature = "Monic",
          definition = function(object)
          {
              SurRatio <- sapply(object@TimeStep,attr, 'SurvivalRatio')
              AbsRatio <- sapply(object@TimeStep,attr, 'AbsorptionRatio')
              Birth <- lengths(lapply(object@TimeStep,attr, 'EmergCan'))
              Died <- sapply(object@TimeStep,attr, 'Died')
              passforward <- sapply(object@TimeStep,attr, 'passforwardRatio')
              height <- rbind(SurRatio, AbsRatio)
              c <- seq(1:length(object@TimeStep))
              xlab <- "Time Steps"

              oldpar <- par(mfrow = c(2,2), no.readonly = TRUE)
              on.exit(par(oldpar))

              barplot(height, beside = TRUE,
                      names.arg = c,
                      ylim = c(0, 1.6),
                      col = c("lightpink","skyblue"),
                      xlab = xlab
              )

              legend("topright",
                     legend = c("SurRatio","AbsRatio"),
                     fill = c("lightpink","skyblue")
              )

              barplot(Birth,
                      ylab = "Number of emergence",
                      xlab = xlab,
                      names.arg = c,
                      col = "skyblue"
              )

              barplot(Died,
                      ylab = "number of disapearence",
                      xlab = xlab,
                      names.arg = c,
                      col = "skyblue"
              )

              plot(passforward,
                   type = "b",
                   xlab = xlab,
                   ylab = "Pass Farward Ratio",
                   col = "skyblue",
                   ylim = c(0,1.20))
              lines(SurRatio,
                    type = "b",
                    col = "light pink")
          }
)
