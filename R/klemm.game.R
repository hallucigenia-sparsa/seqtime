
################# Klemm-Eguiluz #########################
# Code provided by Stefanie Widder
# Reference: Growing Scale-Free Networks with Small World Behavior, Klemm & Eguiluz
# URL: http://arxiv.org/pdf/cond-mat/0107607v1.pdf

# Functions:
# klemm.game
# make.full.clique
# edges.relink
# sparse_IAs
# sparse_democratic
# fill.graph.to.N
# reset.activity.lists
# reduce_deactive_by_identical

# Generate a graph with N nodes
# if verb is true, plot the graph
# clique is including clique.size nodes,
# threshold for link between non-clique nodes is 0.8,
# probability of re-wiring is 0.01
klemm.game<-function(N,verb, clique.size){
  c<-make.full.clique(clique.size)
  c<-fill.graph.to.N(N,0.8,c)
  g<-edges.relink(c,0.01)
  if (verb==TRUE){plot(g)}
  return(g)
}

# make a graph of m nodes where each node
# is connected with each other node
# (fully connected clique)
make.full.clique<-function(m){
  source<-vector()
  target<-vector()
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      source<-c(source,i)
      target<-c(target,j)
      source<-c(source,j)
      target<-c(target,i)
    }
  }
  active<-seq(1,m)
  nw<-graph(rbind(source,target))
  list(g=nw,active_nodes=active,deactive_nodes=vector(),m=m)
}

edges.relink<-function(c,d){
  #work on edgelist
  edges<-get.edgelist(c$g)
  all<-seq(1,c$n)
  #foreach node consider re-linking to a non-adjacent node
  for(i in 1:c$n){
    adj<-as.numeric(V(c$g)[get.adjlist(c$g,mode="out")[[i]]])
    adji<-c(adj,i)
    theothers<-all[!all%in%adji]
    if(length(theothers)>0){
      for(j in 1:(length(adj))){
        x<-runif(1,0.0,1.0)
        if(d>x){
          #remove old edge
          aux=which(edges[,1]==i)
          aux = aux[which(edges[aux,2] == adj[j])]
          edges<-edges[-aux,]
          #add new edge auf namen ebene
          h<-sample(1:length(theothers),1,replace=TRUE)
          edges<-rbind(edges,c(i,theothers[h]))
        }
      }}
  }
  #remove identical edges
  return(graph.edgelist(unique(edges)))

}

# randomly remove edges
sparse_IAs<-function(g,hm){
  edg<-get.edgelist(g)
  # select a random edge subset
  tt<- sample(1:nrow(edg),nrow(edg)*hm)
  # remove selected edges
  edgn<-edg[-tt,]
  ng<-graph.edgelist(edgn)
  #keep nodenumber integrity
  zzz<-length(V(g))-length(V(ng))
  if(zzz>0){
    ng<-add.vertices(ng,zzz)
  }
  return(ng)
}

# randomly remove incoming edges
sparse_democratic<-function(g,hm){
  kedg<-vector()
  edg<-get.edgelist(g)
  adjedg<-get.adjedgelist(g,mode="in")
  tt<- sample(1:length(adjedg),nrow(edg)*hm,replace=TRUE)
  nt<-table(tt)
  cand<-as.numeric(names(nt))#indices, zahlen
  for(i in 1:length(cand)){
    hmy<-as.numeric(nt[i])
    if(hmy>length(adjedg[[cand[i]]])  ){hmy<-length(adjedg[[cand[i]]])}
    kedg<-c(kedg,sample(adjedg[[cand[i]]],hmy))
  }
  edgn<-edg[-kedg,]
  return(graph.edgelist(edgn))
}

# N = node number
# mu = probability of rich get richer
# g = graph
fill.graph.to.N<-function(N,mu,g){
  #names in active
  active<-g$active_nodes
  #deactive<-g$m
  deactive<-vector()
  m<-g$m
  edges<-get.edgelist(g$g)
  # for non-clique nodes, place interactions when their randomly selected strength is above mu
  for(i in (m+1):N){
    for(j in 1:length(active)){
      x<-runif(1,0.0,1.0)	# 1 uniform random number between 0 and 1
      if((mu>x)||(length(deactive)==0)){
        edges<-rbind(edges,c(i,active[j]),c(active[j],i))
      }
      else{
        connected<-"F"
        rd<-reduce_deactive_by_identical(deactive,i,edges)
        while(connected=="F"){
          if(length(rd>0)){
            jj<-sample(1:length(rd),1)
            x<-runif(1,0.0,1.0)
            E<-sum(degree(g$g)[deactive])
            if(((degree(g$g)[rd[jj]])/E)>x){
              edges<-rbind(edges,c(i,rd[jj]),c(rd[jj],i))
              connected<-"T"
            }
          }
          else{
            edges<-rbind(edges,c(i,active[j]),c(active[j],i))
            connected<-"T"
          }
        }
      }
    }#for active
    g$g<-graph.edgelist(edges)
    g<-reset.activity.lists(g,i)
    active<-g$active
    deactive<-g$deactive
  }#for rest nodes
  nw<-graph.edgelist(edges)

  list(g=nw,active_nodes=active,deactive_nodes=deactive,m=g$m,n=N,mu=mu)
}

reset.activity.lists<-function(c,who){
  active<-c(c$active_nodes,who)
  deactive<-c$deactive_nodes
  chosen<-"F"
  while(chosen=="F"){
    j<-sample(1:length(active),1)
    p<-(1/(degree(c$g)[active[j]]))*(1/sum(degree(c$g)[active]))
    x<-runif(1,0.0,1.0)
    if(p>x){
      chosen<-"T"
      deactive<-c(deactive,active[j])
      active<-active[-j]
    }
  }
  c$active_nodes<-active
  c$deactive_nodes<-deactive
  return(c)
}

reduce_deactive_by_identical<-function(deactive,who,edges){
  reduced<-vector()
  aux = which(edges[,1]==who)
  for(i in 1:length(deactive)){
    out<- aux[which(edges[aux,2] == deactive[i])]
    if(length(out)==0){reduced<-c(reduced,deactive[i])}
  }
  return(reduced)
}
