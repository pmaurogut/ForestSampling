
make_population <-function(N,L){
  
  A<-(L*L)/10000
  Npop<-round(N*A)
  res <-data.frame(
    id = 1:Npop,
    x=runif(Npop,0,L),
    y=runif(Npop,0,L),
    diam = round(runif(Npop,5,40),1)
  )
  
  res$ht <- round(5+(res$diam-5)*0.5,1)
  
  res$r_fijo <- 15
  res$area_fijo <-  pi*(res$r_fijo^2)/10000
  res$fac_exp_fijo<- 1/res$area_fijo
  
  res$r_variable <- ifelse(res$diam<15,10,20)
  res$area_variable <-  pi*(res$r_variable^2)/10000
  res$fac_exp_variable <- 1/res$area_variable
  
  res$r_relascopio <- res$diam/2
  res$area_relascopio <-  pi*(res$r_relascopio^2)/10000
  res$fac_exp_relascopio <- 1/res$area_relascopio
  return(res)
}

sampling_points <-function(n,L){
  data.frame(Parc=1:n,x=runif(n,0,L),y=runif(n,0,L))
}

get_trees <- function(population,point,type){
  
  pick <- sqrt((population$x-point$x)^2+(population$y-point$y)^2)<=population[[type]]
  if(all(!pick)){
    return(data.frame(
      type=gsub("radio","_",type),x0=NA,y0=NA,x=NA,y=NA,Parc=point$Parc,
      diam=NA,ht=NA,gi_m2=NA,A_parc_ha=NA,EXP_FAC=NA))
  }else{
    res<-population[pick, ]
    res$Parc<- point$Parc
    res$type <- gsub("radio","r_",type)
    res$x0 <- point$x
    res$y0 <- point$y
    res$gi_m2 <- (pi*(1/4)*res$diam^2)/10000
    res<-res[,c("type","x0","y0","x","y","Parc","diam","ht","gi_m2",type)]
    res$area_parc_ha<-(pi*res[,type]^2)/10000
    res$EXP_FAC <- 1/res$area_parc_ha
    res<-res[order(res$diam,decreasing = TRUE),]
    colnames(res)<-gsub("area_","A_",colnames(res))
    return(res)
  }
  
}

get_n_points <- function(population,points,type){
  
  res<-list()
  for(i in 1:dim(points)[1]){
    res[[i]]<-get_trees(population,points[i,],type)
  }
  map_dfr(res,function(x){x})
}

parametros_interes <- function(poblacion, lado,rotate=TRUE){
  A<-lado*lado/10000
  res<-data.frame(
    Area_ha = A,
    Total_N = length(poblacion[[1]]),
    Total_G = (1/10000)*sum(pi*poblacion$diam^2)/4,
    Total_h = sum(poblacion$ht)
  )
  res$N <- res$Total_N/A
  res$G <- res$Total_G/A
  res$h_media <- mean(poblacion$ht)
  res$dg <- sqrt(mean(poblacion$diam^2))
  poblacion <- poblacion[order(poblacion$diam,decreasing=TRUE),]
  
  if(res$N<100){
    res$Ho<-mean(poblacion$ht)
  }else{
    k <- round(dim(poblacion)[1]*(100/res$N))
    res$Ho<-mean(poblacion[1:k,"ht"])
  }
  if(rotate){
    names<- colnames(res)
    res <- as.data.frame(t(res))
    colnames(res)<-"Valor"
    res$parametro<-names
    res
  }else{
    res
  }
  
}

estimacion <- function(sample,lado,rotate=TRUE){
  print(sample)
  if(is.na(sample$diam[1])){
    res <- data.frame(Parc=sample$Parc[1],Total_N=0,Total_G=0,
                      Total_h=0,N=0,G=0,hmedia=0,dg=0,Ho=0)
    if(rotate){
      names <- colnames(res)
      res <- data.frame(t(res))
      colnames(res)[1]<-"Estimacion"
      res$Variable <- names
      return(res[,c("Variable","Estimacion")])
    }else{
      return(res)
    }
    
    
    
  }
  sample <- sample[order(sample$diam,decreasing = TRUE),]
  sample$cum_sum<-cumsum(sample$EXP_FAC)
  
  A <- (lado*lado)/10000
  res<-data.frame(Parc=sample$Parc[1])
  res$Total_N <- sum(sample$EXP_FAC)*A
  res$Total_G <- sum(sample$EXP_FAC*sample$gi_m2)*A
  res$Total_h <- sum(sample$EXP_FAC*sample$ht)*A
  res$N <- res$Total_N/A
  res$G<- res$Total_G/A
  res$hmedia<- res$Total_h/res$Total_N
  res$dg<-sqrt((res$G/res$N)*(4/pi))
  if(sum(sample$EXP_FAC)>=100){
    last <- which(sample$cum_sum>100)[1]
    s2 <- sample[1:last,]
    s2[last,"EXP_FAC"]<-s2[last,]$cum_sum-100
    res$Ho<-sum(s2$EXP_FAC*s2$ht)/sum(s2$EXP_FAC)
    
  }else{
    res$Ho <- sum(sample$EXP_FAC*sample$ht)/sum(sample$EXP_FAC)
  }
  if(rotate){
    names <- colnames(res)
    res <- data.frame(t(res))
    colnames(res)[1]<-"Estimacion"
    res$Variable <- names
    return(res[,c("Variable","Estimacion")])
  }else{
    return(res)
  }
  
}

n_estimaciones<-function(sample,lado,rotate=FALSE){
  map_dfr(group_split(sample,Parc),estimacion,lado=lado,rotate=rotate)
}

pop_plot <- function(forest_data,lado){
  rect <- data.frame(x=c(0,0,lado,lado,0),y=c(0,lado,lado,0,0))
  ggplot(forest_data) +
    geom_polygon(data=rect,aes(x=x,y=y),col="red",fill="darkgreen",alpha=0.1)+
    geom_circle(aes(x0=x,y0=y,r=diam/20),col="black",fill="burlywood4")+
    xlim(c(-20,lado+20)) + ylim(c(-20,lado+20))+
    coord_fixed(ratio = 1) +
    labs(x="x(m)",y="y(m)") +
    theme_bw(base_size = 16) +
    theme(axis.title = element_blank())
}

plot_selection <- function(p,selected,samp_points,type,tree_center=TRUE,all=FALSE,add_hd=FALSE){
  
  title <- switch(type,
                  r_fijo = "Radio fijo 15m",
                  r_variable = "R anidados d<15cm 10m, d>=15cm 20m",
                  r_relascopio = "Relascopio BAF=1"
  )
  
  if(add_hd){
    p <- p +
      geom_label(aes(x=x,y=y-3,label=paste("d: ",diam)),size=4,fill="darkgreen",alpha=0.3)+
      geom_label(aes(x=x,y=y-8,label=paste("h: ",diam)),size=4,fill="blue",alpha=0.3)
  }
  
  if(all){
    if(type == "r_relascopio"){
      p <- p + geom_circle(aes(x0=x,y0=y,r=.data[[type]],fill=.data[[type]]),alpha=0.4)
    }else{
      p <- p + geom_circle(aes(x0=x,y0=y,r=.data[[type]],fill=factor(.data[[type]])),alpha=0.4)
    }
    
  }
  
  if(!is.na(selected$diam[1])){
    if(tree_center){
      p <- p  + geom_circle(data=selected,aes(x0=x,y0=y,r=.data[[type]]),fill="purple",alpha=0.2)
    }else{
      selected2 <- selected |> group_by(!! sym(type)) |> filter(row_number()==1) |> ungroup()
      p <- p  + geom_circle(data=selected2,aes(x0=x0,y0=y0,r=.data[[type]]),fill="purple",alpha=0.2)  
    }
    p <- p + geom_circle(data=selected,aes(x0=x,y0=y,r=diam/20),col="green",fill="darkgreen",lwd=0.5)
  }
  p <- p + geom_point(data=samp_points[1,],aes(x=x,y=y),shape=13,col="red",size=4)
  p <- p + guides(fill=FALSE)+ggtitle(title)
  p
}

plot_n_selections <- function(p,selected,samp_points,type,tree_center=TRUE,all=FALSE){
  
  selected$Parc <- as.factor(selected$Parc)
  samp_points$Parc <- as.factor(samp_points$Parc)
  title <- switch(type,
                  r_fijo = "Radio fijo 15 m",
                  r_variable = "R anidados d<15 cm 10, d>=15 cm 20m",
                  r_relascopio = "Relascopio BAF=1"
  )
  if(all){
    p <- p + geom_circle(aes(x0=x,y0=y,r=.data[[type]]/20), fill="grey50",alpha=0.2)
  }

  if(!all(is.na(selected$diam))){
    if(tree_center){
      p <- p  + geom_circle(data=selected,aes(x0=x,y0=y,r=.data[[type]],fill=Parc),alpha=0.2)
    }else{
      selected2 <- selected |> group_by(Parc, !!sym(type)) |> filter(row_number()==1) |> ungroup()
      print("Hola")
      print(selected2)
      p <- p  + geom_circle(data=selected2,aes(x0=x0,y0=y0,r=.data[[type]],fill=Parc),alpha=0.2)  
    }
    p <- p + geom_circle(data=selected,aes(x0=x,y0=y,r=diam/20),fill="green")
  }
  p <- p + geom_point(data=samp_points,aes(x=x,y=y,col=Parc),shape=13,size=4)
  p <- p + guides(fill=FALSE,color=FALSE)+ggtitle(title)
  p
}

add_samples_plot<-function(p_int,first,parametro="G"){
  
  names<-p_int$parametro
  p_int <- data.frame(t(p_int[,1,drop=FALSE]))
  colnames(p_int)<-names
  
  print(p_int)

  max <- max(first$G)
  print(max)
  print("Hola")
  
  p <- ggplot(p_int)+
    geom_vline(aes(xintercept=G),col="red")+ylim(c(0,1.5))+xlim(c(-0.1*max,2.1*max))
  
  variation <- data.frame(
    mean=mean(first$G,na.rm=TRUE),
    sd = sd(first$G,na.rm=TRUE)
  )
  variation$xmin <- variation$mean + variation$sd*2
  variation$xmax <- variation$mean-variation$sd*2
  variation$xmin2 <- min(first$G)
  variation$xmax2 <- max(first$G)
  p <- p + geom_point(data=first,aes(x=G,y=0.5),col="red",shape=20,alpha=0.5,size=3)+
    geom_linerange(data=variation,aes(y=1,xmin=xmin2,xmax=xmax2),col="blue")+
    geom_point(data=variation,aes(x=mean,y=1),col="blue",size=5)
  p
  
}


add_samples_n_plots<-function(p_int,all,parametro="G"){
  
  names<-p_int$parametro
  p_int <- data.frame(t(p_int[,1,drop=FALSE]))
  colnames(p_int)<-names
  
  means <- all|> group_by(Rep)|> summarise_all(mean)
  
  means$type_est <- "n-parcelas"
  all$type_est <- "1 parcela"
  
  variation <- data.frame(
    mean=mean(all$G,na.rm=TRUE),
    sd = sd(all$G,na.rm=TRUE)
  )
  variation$xmin <- variation$mean + variation$sd*2
  variation$xmax <- variation$mean-variation$sd*2
  # variation$xmin2 <- min(all$G)
  # variation$xmax2 <- max(all$G)
  variation$type_est <- "1 parcela"
  
  variation2 <- data.frame(
    mean=mean(means$G,na.rm=TRUE),
    sd = sd(means$G,na.rm=TRUE)
  )
  variation2$type_est <- "n-parcelas"
  variation2$xmin <- variation2$mean + variation2$sd*2
  variation2$xmax <- variation2$mean-variation2$sd*2

  
  all <- rbind(means,all)
  all$type_est <- factor(all$type_est,levels=c("1 parcela","n-parcelas"),ordered=TRUE)
  print(all)
  
  variation<- rbind(variation,variation2)
  variation$type_est <- factor(variation$type_est,levels=c("1 parcela","n-parcelas"),ordered=TRUE)
  
  ggplot(all) +
    facet_grid(rows="type_est",scales="fixed")+
    geom_point(aes(x=G,y=0.25,col=type_est,fill=type_est),shape=20,size=4)+
    geom_density(aes(x=G,fill=type_est,col=type_est),alpha=0.4) +
    geom_linerange(data=variation,aes(y=0.5,xmin=xmin,xmax=xmax,col=type_est))+
    geom_vline(data=p_int,aes(xintercept=G),col="red")+
    scale_fill_manual(values=c("1 parcela"="red","n-parcelas"="blue"))+
    scale_color_manual(values=c("1 parcela"="red","n-parcelas"="blue"))

  
}



