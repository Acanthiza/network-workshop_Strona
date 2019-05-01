#set seed to ensure reproducibility
set.seed(123456)

#load the igraph package
library (igraph)

#create your first network
g <- graph( edges=c(1,2, 2,3, 3, 1))


#or load network from file
a<-read.csv('example_network.csv',header=F)

g<-graph.edgelist(as.matrix(a),directed=TRUE)


#visualize the network
plot(g)

#change node color
V(g)$color<-'blue'

plot(g)


#rename nodes
V(g)$name<-c('node1','node2','node3')


#change node size
V(g)$size<-100


#change label color
V(g)$label.color<-'white'


#apply individual changes to nodes
V(g)$size<-c(100,200,300)


#more at https://igraph.org/r/doc/plot.common.html
#e.g.
E(g)$lty<-2


plot(g)


#generate a more interesting random network 
g<-erdos.renyi.game(10,0.5)

plot(g)


#note that layout changes

for (i in 1:10){
  plot(g)
  Sys.sleep(0.5)}



#fix the layout

lay = layout.auto(g)

plot(g,layout = lay)

print (lay)

V(g)$size<-50
for (i in 1:100){
  plot(g,layout = lay)
  V(g)$color<-rainbow(100)[i] #make the example more interesting adding color variation
  Sys.sleep(0.1)}


#generate custom layout
lay_x <- c(1:10)
lay_y <- c(1:10)

lay<-cbind(lay_x,lay_y)

plot(g,layout = lay)


#modify edge curvature to make edges visible

E(g)$curved<-2
plot(g,layout = lay)


#randomize edge curvature

E(g)$curved<-runif(length(E(g)),-2,2)

plot(g,layout = lay)

#randomize edge width
E(g)$width<-runif(length(E(g)),1,3)

plot(g,layout = lay)


#edge direction
E(g)

is.directed(g)


#create bidirectional links
gdm<-as.directed(g, mode = c("mutual"))

E(gdm)

plot(gdm)


#create unidirectional links, with arbitrary direction
gda<-as.directed(g, mode = c("arbitrary"))

E(gda)

plot(gda)


#explore node degree
g<-erdos.renyi.game(1000,0.5)

deg<-degree(g)
h<-hist(deg)
plot(log(h$mids),log(h$density),xlab="log(degree)",ylab="log(frequency)")


g<-sample_smallworld(1,100, 5, 0.05)
deg<-degree(g)
h<-hist(deg)
plot(log(h$mids),log(h$density),xlab="log(degree)",ylab="log(frequency)")


#compare to degree distribution of power law networks
#first play a bit with power law network model parameter
g<-barabasi.game(10,2)
plot(g)


for (node_n in 5:20){
    g<-barabasi.game(node_n,0.1)
  plot(g)
  Sys.sleep(0.5)  
}



g<-barabasi.game(1000,0.2)

deg<-degree(g)
h<-hist(deg)
x<-h$mids[h$density>0]
y<-h$density[h$density>0]
plot(log(x),log(y),xlab="log(degree)",ylab="log(counts)")




###
###Explore other node properties


###Network paths and diameter
diam<-diameter(g)


####shortest paths
g<-graph(c(1,2,2,3,3,4,4,1),directed=TRUE)
plot(g)
shortest.paths(g,4,1,mode='out')
shortest.paths(g,4,1,mode='in')

all_shp<-shortest.paths(g,mode='in')
max(all_shp[all_shp<Inf])

shortest.paths(g,c(1,3),c(2,4,5))


##compute trophic level
###assume nodes with out-degree = 0 (not consuming other nodes) are basal
basal<-which(degree(g,mode='in')==0)

dist_mat<-shortest.paths(g,to=basal,mode='in')

tl<-apply(dist_mat, 1, FUN=min)

tl


#############co-extinction experiments
##random removal
g<-barabasi.game(1000,0.2)
#g<-erdos.renyi.game(1000,0.5)

plot(0,0,col='white',xlim=c(0,1),ylim=c(0,1),xlab='species removed',ylab='diversity')


g1<-g
N<-length(V(g))
x<-c(0)
y<-c(1)
for (i in 1:length(V(g))){
  nodes<-V(g1)
  to_del<-sample(nodes,1)
  g1<-delete_vertices(g1,to_del)
  count<-sum(degree(g1)>0)
  x<-c(x,i/N)
  y<-c(y,count/N)  
}

lines(x,y,lwd=2,col='black')



################removal according to node degree
g1<-g
N<-length(V(g))
x<-c(0)
y<-c(1)
for (i in 1:length(V(g))){
  to_del<-which(degree(g1)==max(degree(g1)))
  g1<-delete_vertices(g1,to_del)
  count<-sum(degree(g1)>0)
  x<-c(x,i/N)
  y<-c(y,count/N)  
}

lines(x,y,lwd=2,col='red')

###########removal according to pagerank
#see page.rank output structure:
page.rank(g)



g1<-g
N<-length(V(g))
x<-c(0)
y<-c(1)
for (i in 1:length(V(g))){
  pr<-page.rank(g1,directed=F)$vector
  to_del<-which(pr==max(pr))
  g1<-delete_vertices(g1,to_del)
  count<-sum(degree(g1)>0)
  x<-c(x,i/N)
  y<-c(y,count/N)  
}

lines(x,y,lwd=2,col='blue')




####Empirical bipartite network
###Go to http://www.web-of-life.es/map.php


pol_mat<-as.matrix(read.csv("M_PL_029.csv",header=TRUE,row.names=1))
n_plant<- dim(pol_mat)[1]
n_pol<- dim(pol_mat)[2]


g<-graph_from_incidence_matrix(pol_mat, directed = FALSE)

plot(g,layout=layout.bipartite)

V(g)$name<-''
E(g)$arrow.size=0.1

V(g)$color <-1+V(g)$type*1

colnames(pol_mat)
row.names(pol_mat)

#rows are plants; node type == FALSE 

###Simulate co-extinction by removing plants in random order
plot(0,0,xlim=c(0,1),ylim=c(0,1),xlab = 'plant diversity',ylab='pollinator diversity')

g1<-g
x<-c(0)
y<-c(1)
for (i in 1:n_plant){
    to_del<-sample(list(which(V(g1)$type == FALSE))[[1]],1)
    g1<-delete_vertices(g1,to_del)

    to_del_co_ex<-which(degree(g1)==0)
    g1<-delete_vertices(g1,to_del_co_ex)
    
  count<-sum(V(g1)$type)
  x<-c(x,i/n_plant)
  y<-c(y,count/n_pol)
}

lines(x,y,lwd=2,col='black')


###from te most connected to the least connected plant
g1<-g
x<-c(0)
y<-c(1)
for (i in 1:n_plant){
  deg_plant<-degree(g1)[which(V(g1)$type == FALSE)]
  to_del<-sample(list(which(V(g1)$type == FALSE & degree(g1)==max(deg_plant)))[[1]],1)

  g1<-delete_vertices(g1,to_del)
  
  to_del_co_ex<-which(degree(g1)==0)
  g1<-delete_vertices(g1,to_del_co_ex)
  
  count<-sum(V(g1)$type)
  x<-c(x,i/n_plant)
  y<-c(y,count/n_pol)
  
}


lines(x,y,lwd=2,col='red')


###from the least connected to the most connected
g1<-g
x<-c(0)
y<-c(1)
for (i in 1:n_plant){
  deg_plant<-degree(g1)[which(V(g1)$type == FALSE)]
  to_del<-sample(list(which(V(g1)$type == FALSE & degree(g1)==min(deg_plant)))[[1]],1)
  
  g1<-delete_vertices(g1,to_del)
  
  to_del_co_ex<-which(degree(g1)==0)
  g1<-delete_vertices(g1,to_del_co_ex)
  
  count<-sum(V(g1)$type)
  x<-c(x,i/n_plant)
  y<-c(y,count/n_pol)
  
}


lines(x,y,lwd=2,col='blue')









#######################
###Empirical food web

#obtain data
library(rglobi)

#create a bounding box around an area of interest
serengeti <- get_interactions_in_area( bbox=c(34.20, -2.62, 35.31, -1.48))

table(serengeti$interaction_type)

#limit the interactions to prey-predator
pred <- serengeti[(serengeti$interaction_type=="preysOn")|(serengeti$interaction_type=="eats"),]
colnames(pred)

#identify which species are plants
plants<-pred[grepl('Plantae', pred[,8]),7]


#generate subtable with prey-->predator; col 2 is predator; col 7 is prey
fweb<-pred[,c(7,2)]

#convert to igraph directed network
g<-graph.edgelist(as.matrix(fweb),directed=TRUE)


#explore trophic structure

basal<-which(degree(g,mode='in')==0) #no incoming food/resource

#find and remove which basal are not plants
not_rooted<-setdiff(V(g)$name[basal],plants)
length(not_rooted)
g<-delete_vertices(g,not_rooted)

#check

basal<-which(degree(g,mode='in')==0) #no incoming food/resource
not_rooted<-setdiff(V(g)$name[basal],plants)
length(not_rooted)

while(length(not_rooted)>0){
  g<-delete_vertices(g,not_rooted)
  basal<-which(degree(g,mode='in')==0) #no incoming food/resource
  not_rooted<-setdiff(V(g)$name[basal],plants)
}


length(not_rooted)

##########################################
#reminder of how shortest.path works...
g1<-graph(c(1,2,2,3,3,4,4,1),directed=TRUE)
plot(g1)
out_1_to_4<-as.numeric(shortest.paths(g1,1,4,mode='out'))
in_1_to_4<-as.numeric(shortest.paths(g1,1,4,mode='in'))
out_1_to_4
in_1_to_4
###########################################

dist_mat<-shortest.paths(g,to=basal,mode='in') #move from focal node towards basal resources
tl<-apply(dist_mat, 1, FUN=min)

tl

####Visualization
node_n<-length(V(g))
lay_y<-tl
max_freq<-max(table(tl))

#lay_x<-runif(node_n,0,max_freq)
lay_x<-rev(order(degree(g)))



lay<-cbind(lay_x,lay_y)

plot(g,layout = lay)

V(g)$name<-''
V(g)$size<-10*degree(g)/max(degree(g))
V(g)$color<-tl
E(g)$arrow.size<-0.1

plot(g,layout = lay)


#g<-simplify(g)
E(g)$curved<-runif(length(E(g)),-0.5,0.5)
pdf("network.pdf")
plot(g,layout = lay)
dev.off()


########################
#Co-extinctions in food webs
#identify basal species and others
basal<-which(degree(g,mode='in')==0) #no incoming food/resource
dist_mat<-shortest.paths(g,to=basal,mode='in') #move from focal node towards basal resources
tl<-apply(dist_mat, 1, FUN=min)
V(g)$is_basal<-'no'
V(g)[basal]$is_basal<-'yes'
V(g)$trophic_level<-(tl+1/degree(g))


plot(0,0,col='white',xlim=c(0,1),ylim=c(0,1),xlab='species removed',ylab='diversity')


g1<-g

N<-length(V(g))
x<-c(0)
y<-c(1)
for (i in 1:N){
  if (length(V(g1))>0){
    nodes<-V(g1)
    to_del<-sample(1:length(nodes),1)
    while (length(to_del)>0){
      g1<-delete_vertices(g1,to_del)
      to_del<-setdiff(which(degree(g1,mode='in')==0),which(V(g1)$is_basal=='yes'))
    }
  }
  count<-length(V(g1))
  x<-c(x,i/N)
  y<-c(y,count/N)  
}

lines(x,y,lwd=2,col='black')

###Remove also species remaining with no paths to basal resources

g1<-g

N<-length(V(g))
x<-c(0)
y<-c(1)
for (i in 1:N){
  nodes<-V(g1)
  if (length(nodes)>0){
      to_del<-sample(1:length(nodes),1)
      while (length(to_del)>0){
          g1<-delete_vertices(g1,to_del)
          rem_bas<-which(V(g1)$is_basal=='yes')
          dist_mat<-shortest.paths(g1,to=rem_bas,mode='in')
          tl_rem<-apply(dist_mat, 1, FUN=min,na.rm=TRUE)
          to_del_co_ex<-which(tl_rem == Inf)  
          g1<-delete_vertices(g1,to_del_co_ex)
          to_del<-setdiff(which(degree(g1,mode="in")==0),which(V(g1)$is_basal=='yes'))
          }
      }
  count<-length(V(g1))#sum(degree(g1)>0)
  x<-c(x,i/N)
  y<-c(y,count/N)  
}


lines(x,y,lwd=2,col='black',lty=2)


################trophic level


###from highest to lowest - best case scenario
g1<-g
N<-length(V(g))
x<-c(0)
y<-c(1)
for (i in 1:N){
    nodes<-V(g1)
    if (length(nodes)>0){
      to_del<-which(V(g1)$trophic_level==max(V(g1)$trophic_level))[1]
      while (length(to_del)>0){
        g1<-delete_vertices(g1,to_del)
        rem_bas<-which(V(g1)$is_basal=='yes')
        dist_mat<-shortest.paths(g1,to=rem_bas,mode='in')
        tl_rem<-apply(dist_mat, 1, FUN=min, na.rm=TRUE)
        to_del_co_ex<-which(tl_rem == Inf)  
        g1<-delete_vertices(g1,to_del_co_ex)
        to_del<-setdiff(which(degree(g1,mode='in')==0),which(V(g1)$is_basal=='yes'))
      }
    }
    count<-length(V(g1))#sum(degree(g1)>0)
    x<-c(x,i/N)
    y<-c(y,count/N)  
  }
  
lines(x,y,lwd=2,col='blue')



#from lowest to highest - worst case scenario
g1<-g
N<-length(V(g))
x<-c(0)
y<-c(1)
for (i in 1:N){
  nodes<-V(g1)
  if (length(nodes)>0){
    to_del<-which(V(g1)$trophic_level==min(V(g1)$trophic_level))[1]
    while (length(to_del)>0){
      g1<-delete_vertices(g1,to_del)
      rem_bas<-which(V(g1)$is_basal=='yes')
      dist_mat<-shortest.paths(g1,to=rem_bas,mode='in')
      tl_rem<-apply(dist_mat, 1, FUN=min, na.rm=TRUE)
      to_del_co_ex<-which(tl_rem == Inf)  
      g1<-delete_vertices(g1,to_del_co_ex)
      to_del<-setdiff(which(degree(g1,mode='in')==0),which(V(g1)$is_basal=='yes'))
    }
  }
  count<-length(V(g1))#sum(degree(g1)>0)
  x<-c(x,i/N)
  y<-c(y,count/N)  
}

lines(x,y,lwd=2,col='red')

