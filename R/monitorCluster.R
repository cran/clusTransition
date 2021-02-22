#' Monitor Transitions in Cluster Solutions.
#'
#' @description Model and trace the evolution of clusters evolving over time in cumulative
#' datasets. A typical call to \code{Transition()} function involves three essential pieces:
#' the data input \code{(listdata, listclus, overlap)}, choice of window \code{swSize},
#' and the threshold parameters. The function either receive a list of datasets arriving at
#' time points \code{t_1, t_2, t_3, ..., t_n} respectively, list of clustering solutions
#' extracted from cumulative datasets at successive time points, or list of objects of class
#' \code{OverLap} (see \strong{Details}).
#'
#' @param listdata List of numeric matrices containing datasets \code{d_1, d_2, ..., d_n},
#' or a list of objects that can be coerced to such matrices, for instance, data frames.
#' Each element of the list contain dataset \code{d_i} evolving at corresponding time point
#' \code{t_i}. The number of clusters in each accumulative data matrix is specified by the
#' argument \code{k}.
#'
#' @param swSize Integer value (1, length(listdata)) indicating size of the sliding window. As time goes
#' by, each window consist only objects that fall in the interval [t-swSize+1, t], while older objects
#' are discarded. The default value of \code{swSize = 1} indicate landmark window model, where objects over
#' the entire history are included i.e. [1, t]. Size of sliding window can only be provided if \code{listdata}
#' arguments is choosen. If there are total \code{n} time stamps and a window of size \code{swSize} is
#' selected then entire history would be devided into \code{n-swSize+2} window panes.
#'
#' @param Overlap A list of objects as produced by the \code{Overlap()} method. The object contains a matrix of similarity
#' indices between clusters, and the summaries of clusters extracted at first and second clustering.
#'
#' @param listclus \code{listclus} is a list of nested lists containing clustering solutions {\eqn{\xi_1, \xi_2, ..., \xi_n}} at
#' time points \code{{t1, t2,···, tn}} respectively, and having the same length as the number of time points. The \code{i^th}
#' element of \code{listclus} is a nested list that contain set of clusters as matrices at corresponding time point \code{t_i}
#' i.e. \eqn{\xi_i = {X1, X2,···, Xki}}. For more details, \emph{see} \strong{Examples}.
#'
#' @param typeind Type indicator. \code{typeind = 1} indicates that the raw data is provided in
#' \code{listdata} argument, \code{typeind = 2} indicates that the \code{OverLap} objects are provided,
#' whereas \code{typeind = 3} indicates that list of clusters are provided using \code{listclus} argument.
#'
#' @param Survival_thrHold A numeric value (0,1) indicating minimum threshold value for survival of clusters.
#'
#' @param Split_thrHold A numeric value (0,1) indicating minimum threshold value for split of clusters.
#'
#' @param location_thrHold A numeric value (0,1) indicating minimum threshold value for shift in location of survived clusters.
#'
#' @param density_thrHold A numeric value (0,1) indicating minimum threshold value for changes in density of Survived clusters.
#'
#' @param k Numeric Vector of length \code{vector("numeric", length = n-swSize+2)}. In the case of landmark window, its length
#' is \code{n}, whereas in case of sliding window model its length is \code{n-swSize+2}, where \code{n} is the number of time points
#' and \code{swSize} is the size of the sliding window. This argument should only be provided if \code{listdata} argument is chosen.
#'
#' @import methods
#'
#' @return Returns A list of class \code{Monic}.
#'
#' \item{Survive}{Number of clusters survived.}
#' \item{Merged}{Number of clusters merged.}
#' \item{Split}{Number of clusters split.}
#' \item{Died}{Number of clusters disappeared.}
#' \item{new.Emerged}{Number of newly emerged clusters, which are not upshot of any external transition.}
#' \item{SurvivalCanx}{A vector of integers indicating candidates from the first clustering
#' survived to the latter time stamp}
#' \item{SurvivalCany}{A vector of integers indicating candidates of second clustering, that
#' clinch the survival candidates from first clustering.}
#' \item{SplitCanx}{A vector of integers indicating candidate(s) that split into various daughter clusters from
#'  first clustering.}
#' \item{SplitCany}{List of integer vector(s) designating candidates appeared, as a result of splits
#'  from first clustering.}
#' \item{MergeCanx}{List of integer vector(s) designating Candidates that spliced together to form
#' new clusters. Each element of the list gives candidates that merge together to form one.}
#' \item{MergeCany}{Vector of integers designating candidates that emerged, as a result of merger
#' of different candidates from first clustering.}
#' \item{EmergCan}{Vector of integers contain Newly emerged candidates, which are not result of
#' any external transition.}
#' \item{SurvivalRatio}{The Ratio of survived clusters at second clustering to the total number
#' of clusters at first clustering.}
#' \item{AbsorptionRatio}{Ratio of number of merged clusters to total number of clusters at first
#' clustering.}
#' \item{passforwardRatio}{Sum of SurvivalRatio and AbsorptionRatio. This gives the ratio of
#' clusters that is also present at second clustering either in the form of survival or absorption.}
#' \item{Overlap}{A numeric matrix containing overlap of the two clustering. The rows of matrix
#' indicate first clustering, while columns indicate second clustering.}
#' \item{Centersx}{A matrix of cluster centers from first clustering.}
#' \item{Centersx}{A matrix of cluster centers from second clustering.}
#' \item{rx}{A numeric vector containing radius of each cluster from first clustering.}
#' \item{ry}{A numeric vector containing radius of each cluster from second clustering.}
#' \item{avgDisx}{A numeric vector containing average distance of points in a cluster from its center in first clustering.}
#' \item{avgDisy}{A numeric vector containing average distance of points in a cluster from its center in second clustering.}
#' \item{ShiftLocCan}{A vector of integers comprises of Survived candidates with shift in location.}
#' \item{NoShiftLocCan}{A vector of integers comprises of Survived candidates with no shift in location.}
#' \item{MoreCompactCan}{A Vector of integers comprises of Survived candidates, which becomes more compact.}
#' \item{MoreDiffuseCan}{A Vector of integers comprises of Survived candidates, which becomes more diffuse.}
#' \item{NoChangeCompactCan}{A Vector of integers comprises of Survived candidates, with no changes in compactness.}
#' \item{Location.diff}{A numeric vector containing Distance between the centers of survived clusters.}
#' \item{Compactness.diff}{A numeric vector containing Difference between compactness of survived clusters.}
#' \item{Cluster_Tracex}{A vector containing result of each cluster from first clustering.}
#' \item{Cluster_Tracey}{A Vector representing result of each cluster from second clustering.}
#' \item{clusterMem}{A vector of integers (from 1 to k) indicating the point to which cluster it is allocated from second clusterig.}
#'
#' @details The \code{Transition()} function apply 'MONIC' algorithm presented by Spiliopoulou et.al (2006) to trace
#' changes in cluster solutions of dynamic data sets. The changes includes two types of transition i.e. External transition
#' and Internal transition. External Transition consist of 'Survive', 'Split', 'Merge', 'Disappeared' and 'newly emerged' candidates,
#' while Internal transition consist of changes in location and cohesion of the survived clusters. The \code{listdata} argument
#' allow user to import dynamic datasets as a list of matrices or data frames, where each element of the list is a matrix containing
#' data set at a single time point. Each dataset are clustered by 'kmeans' algorithm using default settings of \code{cclust()} function
#' from \code{flexclust} package. The number of clusters at each time stamp can be import by \code{k} argument of the function,
#' which is a vector of integers encompassing number of partitions in corresponding datasets of \code{listdata} argument. Once the datasets are
#' clustered, the 'Overlap' matrices in clustering at consecutive time stamps are calculated. The Overlap matrix is
#' calculated by using algorithm presented by Ntoutsi, I., et.al (2012). These 'Overlap' matrices are used to trace the
#' transitions occurred in cluster solutions.
#' Alternatively, the user can directly import list of 'Overlap' matrices between consecutive clustering. The Overlap
#' matrix can be calculated using \code{Overlap(obj, e1, e2)} method of the package, where 'obj' is the object of class
#' \code{OverLap} and e1, e2 are any clustering at time stamp i and j respectively.
#' As a third option user can provide list of clusters at each data point utilizing \code{listclus} argument. Each element
#' of the \code{listclus} is a nested list, which holds clusters at a single time stamp.
#'
#' @export Transition
#'
#' @examples
#'
#' ### Example 1: typeind = 1 (listdata Argument)
#'
#' d1 <- Data2D[[1]][c("X1", "X2")]
#' d2 <- Data2D[[2]][c("X1", "X2")]
#' d3 <- Data2D[[3]][c("X1", "X2")]
#'
#' listdata <- list(d1, d2, d3)
#'
#' p <- Transition(listdata = listdata, swSize = 1, typeind = 1, Survival_thrHold = 0.8,
#'                 Split_thrHold = 0.3, density_thrHold = 0.3, location_thrHold = 0.3, k = c(3,3,2))
#'
#' ### Example 2: typeind = 3 (listclus Argument)
#'
#' D1 <- d1
#' D2 <- merge(d1, d2, all.x = TRUE, all.y = TRUE)
#' D3 <- merge(D2, d3, all.x = TRUE, all.y = TRUE)
#'
#' set.seed(10)
#' f1 <- kmeans(D1, 3)
#' C1 <- list()
#' for(i in 1:3)C1[[i]] <- D1[f1$cluster == i, ]
#' f2 <- kmeans(D2, 3)
#' C2 <- list()
#' for(i in 1:3)C2[[i]] <- D2[f2$cluster == i, ]
#' f3 <- kmeans(D3, 2)
#' C3 <- list()
#' for(i in 1:2)C3[[i]] <- D3[f3$cluster == i, ]
#'
#' listclus <- list(C1, C2, C3)
#'
#' p <- Transition(listclus = listclus, typeind = 3, Survival_thrHold = 0.8,
#'                 Split_thrHold = 0.3, density_thrHold = 0.3, location_thrHold = 0.3)
#'
#' ### Example 3: typeind = 3 (Overlap Argument)
#'
#' obj <- new("OverLap")
#' Overlap1 <- Overlap(obj, e1 = C1, e2 = C2)
#' Overlap2 <- Overlap(obj, e1 = C2, e2 = C3)
#'
#' Overlap <- list(Overlap1, Overlap2)
#' p <- Transition(Overlap = Overlap, typeind = 2, Survival_thrHold = 0.8,
#'                 Split_thrHold = 0.3, density_thrHold = 0.3, location_thrHold = 0.3)
#'
#'
#' @references
#'
#' Spiliopoulou, M., Ntoutsi, I., Theodoridis, Y., Schult, R. MONIC: modeling and monitoring cluster
#' transitions. In: Eliassi-Rad, T., Ungar, L. H., Craven, M., Gunopulos, D. (eds.) ACM SIGKDD 2006, pp. 706-711. ACM, Philadelphia (2006).
#'
#====================================================================================

Transition <- function(listdata, swSize = 1, Overlap = NULL, listclus = NULL, typeind = 1, Survival_thrHold = 0.7,
                       Split_thrHold = 0.3, location_thrHold = 0.3, density_thrHold = 0.3, k = NULL)
{
    R<-new("Monic")

    #... if 'listdata' argument is provided
    if(typeind == 1)
    {
        if(!is.list(listdata))
            stop("Provide 'data frames' or 'matrices' of datasets as list")
        if(is.null(k))
            stop("Provide 'k' as vector of integers, The No of Centers at each accumulative data point")
        if((length(k) != (length(listdata)-swSize+2)) & swSize != 1)
            stop("Number of accumulative datasets and length of 'k' does not match")
        listdata <- lapply(listdata, function(x) {if(any(class(x)=="data.frame")) as.matrix(x) else x})
        y <- listdata[[1]]
        l1 <- Clusters(new("Clustering"), x = y, k = k[1])
        ifelse(swSize == 1, (itr = length(listdata)-1), (itr = length(listdata)-swSize+1))
        for(j in 1:itr)
        {
            listd <- listdata[c(j:(j+swSize-1))]
            ifelse(swSize == 1,
                   (y <- as.matrix(merge(y, listdata[[j+1]], all.x = TRUE, all.y = TRUE))),
                   (y <- as.matrix(Reduce(function(x,y) (merge(x,y, all = TRUE)),listd))))

            l2 <- Clusters(
                new("Clustering"),
                x = y, k = k[j+1]
            )

            l3 <- Overlap(
                new("OverLap"),
                e1 = l1, e2 = l2
            )

            l4 <- extTransitionCount(
                new("TransitionCount",
                    Survival_thrHold = Survival_thrHold,
                    Split_thrHold = Split_thrHold,
                    l3
                )
            )

            l5 <- extTransitionCan(
                new("TransitionCan",
                    l4
                )
            )

            l6 <- internalTransition(
                new("intTransitionCan",
                    Location_thrHold = location_thrHold,
                    Density_thrHold = density_thrHold,
                    l5
                )
            )
            R@TimeStep[[j]] <- l6
            l1 <- l2
        }
        return(R)
    }
    #... if list of 'Overlap, matrices is provided
    else if(typeind == 2)
    {
        if((is.null(Overlap)) & (!is.list(Overlap)))
            stop("Provide 'Overlap' matrices as list")
        for(i in 1:length(Overlap))
        {
            l1 <- extTransitionCount(
                new("TransitionCount",
                    Survival_thrHold = Survival_thrHold,
                    Split_thrHold = Split_thrHold,
                    Overlap[[i]]
                )
            )

            l2 <- extTransitionCan(
                new("TransitionCan",
                    l1
                )
            )

            l3 <- internalTransition(
                new("intTransitionCan",
                    Location_thrHold = location_thrHold,
                    Density_thrHold = density_thrHold,
                    l2
                )
            )
            R@TimeStep[[i]] <- l3
        }
        return(R)
    }
    #... if 'listclus' argument is selected
    else if(typeind == 3)
    {
        if((is.null(listclus)) & (!is.list(listclus)))
            stop("Provide 'Cluster solutions' as list")
        for(j in 1:(length(listclus)-1))
        {
            l1 <- Overlap(
                new("OverLap"),
                e1 = listclus[[j]], e2 = listclus[[j+1]]
            )

            l2 <- extTransitionCount(
                new("TransitionCount",
                    Survival_thrHold = Survival_thrHold,
                    Split_thrHold = Split_thrHold,
                    l1
                )
            )

            l3 <- extTransitionCan(
                new("TransitionCan",
                    l2
                )
            )

            l4 <- internalTransition(
                new("intTransitionCan",
                    Location_thrHold = location_thrHold,
                    Density_thrHold = density_thrHold,
                    l3
                )
            )
            R@TimeStep[[j]] <- l4
        }
        return(R)
    }
}
#==============================================================================#

