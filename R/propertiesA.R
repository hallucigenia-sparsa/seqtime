#' Print properties of the interaction matrix
#'
#' The interaction matrix is converted into a graph and a number
#' of graph properties (modularity, average clustering coefficient, average path length) are computed using igraph functions.
#' @param A interaction matrix
#' @param plot.degree plot the degree distribution
#' @export

propertiesA<-function(A, plot.degree=FALSE){
  # has to be 0 or 1
  A[A!=0]=1
  g=graph.adjacency(A, mode="directed")
  fc <- fastgreedy.community(as.undirected(g),weights=NULL)
  print("Fastgreedy modularity of final network: ")
  print(max(fc$modularity))
  if(plot.degree==TRUE){
    plot(degree.distribution(g,cum=T),log="xy")
  }
  print(paste("Average clustering coefficient of final network",transitivity(g, type="average")))
  print(paste("Average path length of final network",average.path.length(g)))
}
