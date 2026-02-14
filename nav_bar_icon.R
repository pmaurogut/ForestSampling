library(lidR)
library(purrr)
trees <- list.files("/media/TLM_PROJECT2/1_TLS_DRON/data/LiDAR/7_MODELING_REPO/cif_12/in_12_m/",pattern = ".las",full.names = T)

ordered<-map_dfr(trees,function(x){
  
  h<-readLASheader(x)
  
  res <- data.frame(file=x,zmin=h$`Min X`,zmax=h$`Max Z`)
  res$h <- res$zmax-res$zmin
  res
  
})
ordered<-ordered[order(ordered$h),]
ordered$offset <-rank(ordered$h)

make_figure <- function(ordered,where){
  
  list_of_las<- list()
  for(i in 1:dim(ordered)[1]){
    las <- readLAS(ordered[i,"file"])
    las <- las@data
    las$X <- las$X-mean(las$X) -1.5*i
    las$Y <- las$Y - mean(las$Y)
    las$Z <- las$Z - min(las$Z)
    las<-LAS(las)
    list_of_las[[i]]<-las
    png(gsub(".png",paste(i,".png",sep="_"),where),width=100,height=160,bg = "transparent" )
    par(mar=rep(0,4),oma=rep(0,4))
    plot(las$X,las$Z,pch=10,col=rgb(las$R/max(las$R),las$G/max(las$G),las$B/max(las$B),0.2),bg=rgb(0,0,0,0.2),cex=0.1,
         xaxt='n', axes=FALSE,ann=FALSE)
    dev.off()
  }
  
  big_las <- do.call(rbind,list_of_las)
  big_las <- decimate_points(big_las,random(500))
  
  
  png(where,width=1000,height=160,bg = "transparent" )
  par(mar=rep(0,4),oma=rep(0,4))
  plot(big_las$X,big_las$Z,pch=10,col=rgb(0,0,0,0.2),bg=rgb(0,0,0,0.2),cex=0.1,
       xaxt='n', axes=FALSE,ann=FALSE)
  dev.off()
}




