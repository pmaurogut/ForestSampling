ff <- function(clust){
  clust <- clust
  function(x,y){
    0.25*(2+sin(2*pi*(x/100))+cos(2*pi*(y/100)))
  }
  
}
make_population <-function(N,L){
  
  A<-(L*L)/10000
  Npop<-round(N*A)
  x<-runif(Npop,0,L)
  y<-runif(Npop,0,L)
  cluster <- runif(1,min=0.25,max=1)
  
  # if(cluster==0){
  #   mean <- sin(2*pi*L/100)
  #   
  #   dn_cm <- round(
  #     ifelse(type,pmax(5,rnorm(Npop,10,10)),
  #            pmin(80,rnorm(Npop,50,20))),1)
  # }else{
  #   
  #   
  #   alpha<-ff(clust)(x,y)
  #   d1 <- pmax(5,rnorm(Npop,10,cluster))
  #   d2 <- max(0,pmin(90,rnorm(Npop,70,cluster)))
  #   dn_cm <- (alpha)*d1+(1-alpha)*d2
  # }
  alpha<-ff(clust)(x,y)
  d1 <- pmax(5,rnorm(Npop,10,1.5))
  d2 <- max(0,pmin(90,rnorm(Npop,60,1.5)))
  dn_cm <-cluster*((alpha)*d1+(1-alpha)*d2)+ runif(Npop,5,90)*(1-cluster)
  
  res <-data.frame(
    id = 1:Npop,
    x=x,
    y=y,
    dn_cm = dn_cm
  )
  
  res$ht_m <- round(5+(res$dn_cm-5)*0.5,1)
  res$gi_m2 <- pi*(res$dn_cm/200)^2
  res$vcc_m3 <- res$ht_m*res$gi_m2*0.7
  
  res$r_fijo <- 15
  res$area_fijo <-  pi*(res$r_fijo^2)/10000
  res$fac_exp_fijo<- 1/res$area_fijo
  
  res$r_variable <- ifelse(res$dn_cm<15,10,20)
  res$area_variable <-  pi*(res$r_variable^2)/10000
  res$fac_exp_variable <- 1/res$area_variable
  
  res$r_relascopio <- res$dn_cm/2
  res$area_relascopio <-  pi*(res$r_relascopio^2)/10000
  res$fac_exp_relascopio <- 1/res$area_relascopio
  return(res)
}

make_population_old<-function(N,L){
  
  A<-(L*L)/10000
  Npop<-round(N*A)
  res <-data.frame(
    id = 1:Npop,
    x=runif(Npop,0,L),
    y=runif(Npop,0,L),
    dn_cm = round(runif(Npop,5,40),1)
  )
  
  res$ht_m <- round(5+(res$dn_cm-5)*0.5,1)
  res$gi_m2 <- pi*(res$dn_cm/200)^2
  res$vcc_m3 <- res$ht_m*res$gi_m2*0.7
  
  res$r_fijo <- 15
  res$area_fijo <-  pi*(res$r_fijo^2)/10000
  res$fac_exp_fijo<- 1/res$area_fijo
  
  res$r_variable <- ifelse(res$dn_cm<15,10,20)
  res$area_variable <-  pi*(res$r_variable^2)/10000
  res$fac_exp_variable <- 1/res$area_variable
  
  res$r_relascopio <- res$dn_cm/2
  res$area_relascopio <-  pi*(res$r_relascopio^2)/10000
  res$fac_exp_relascopio <- 1/res$area_relascopio
  return(res)
}

parametros_interes <- function(poblacion, lado,rotate=TRUE){
  A<-lado*lado/10000
  res<-data.frame(
    Area_ha = A,
    Total_N = length(poblacion[[1]]),
    Total_G = (1/10000)*sum(pi*poblacion$dn_cm^2)/4,
    Total_h = sum(poblacion$ht_m),
    Total_V = sum(poblacion$vcc_m3)
  )
  res$N <- res$Total_N/A
  res$G <- res$Total_G/A
  res$V <-res$Total_V/A
  
  res$h_media <- mean(poblacion$ht_m)
   
  res$dg <- sqrt(mean(poblacion$dn_cm^2))
  poblacion <- poblacion[order(poblacion$dn_cm,decreasing=TRUE),]
  
  if(res$N<=100){
    res$ho<-mean(poblacion$ht_m)
  }else{
    k <-round(100*A)
    res$ho<-mean(poblacion[1:k,"ht_m"])
  }
  if(rotate){
    names<- colnames(res)
    res <- as.data.frame(t(res))
    colnames(res)<-"Valor"
    res$parametro<-names
    res[,c(2,1)]
  }else{
    res
  }
  
}

sampling_points <-function(n,L){
  data.frame(Parc=1:n,xp=runif(n,0,L),yp=runif(n,0,L))
}

get_trees <- function(population,point,type){
  
  pick <- sqrt((population$x-point$xp)^2+(population$y-point$yp)^2)<=population[[type]]
  if(all(!pick)){
    return(data.frame(
      Type=type,
      Parc=point$Parc,xp=point$xp,yp=point$yp,x=NA,y=NA,
      dn_cm=NA,ht_m=NA,gi_m2=NA,vcc_m3=NA,radio_sel_m=NA,A_parc_ha=NA,EXP_FAC=NA))
  }else{
    
    res<-population[pick, ]
    res$Parc<- point$Parc
    res$Type <- type
    res$xp <- point$xp
    res$yp <- point$yp
    
    res$gi_m2 <- (pi/4)*(res$dn_cm/100)^2
    res$radio_sel_m <- res[,type]
    res$A_parc_ha<-(pi*res[,type]^2)/10000
    res$EXP_FAC <- 1/res$A_parc_ha
    res<-res[order(res$dn_cm,decreasing = TRUE),]
    res<-res[,c("Type","Parc","xp","yp",
                "x","y","dn_cm","ht_m","gi_m2","vcc_m3",
                "radio_sel_m","A_parc_ha","EXP_FAC")]
    colnames(res)<-gsub("area_","A_",colnames(res))
    return(res)
  }
  
}

get_all_trees <- function(population,points){
  
  points_list <- group_split(points,Parc)
  map_dfr(c("r_fijo","r_variable","r_relascopio"),
          function(r){
            
            map_dfr(points_list,function(x,population,type){
              get_trees(population,x,type)
            },population=population,type=r)

          })

}


estimacion <- function(sample,lado,rotate=TRUE){
  
  A <- (lado*lado)/10000
  res <- data.frame(Type=sample$Type[1],Parc=sample$Parc[1],
                    xp=sample$xp[1],yp=sample$yp[1],
                    Total_N=0,Total_G=0,Total_V=0,Total_h=0,N=0,G=0,V=0,h_media=NA,dg=NA,ho=NA)
  if(!is.na(sample$dn_cm[1])){

    sample <- sample[order(sample$dn_cm,decreasing = TRUE),]
    sample$cum_sum<-cumsum(sample$EXP_FAC)
    res$Total_N <- sum(sample$EXP_FAC)*A
    res$Total_G <- sum(sample$EXP_FAC*sample$gi_m2)*A
    res$Total_h <- sum(sample$EXP_FAC*sample$ht_m)*A
    res$Total_V <- sum(sample$EXP_FAC*sample$vcc_m3)*A
    res$N <- res$Total_N/A
    res$G<- res$Total_G/A
    res$V <- res$Total_V/A
    res$h_media<- res$Total_h/res$Total_N
    res$dg<-sqrt((res$G/res$N)*(4/pi))*100
    
    if(sum(sample$EXP_FAC)>100){
      last <- which(sample$cum_sum>100)[1]
      s2 <- sample[1:last,]
      s2[last,"EXP_FAC"]<-s2[last,]$cum_sum-100
      res$ho<-sum(s2$EXP_FAC*s2$ht_m)/sum(s2$EXP_FAC)
      
    }else{
      res$ho <- sum(sample$EXP_FAC*sample$ht_m)/sum(sample$EXP_FAC)
    }
  }
  
  
  
  if(rotate){
    names <- colnames(res)
    res <- data.frame(t(res[,-c(1:5)]))
    colnames(res)[1]<-"Estimacion"
    res$Variable <- names
    res$Parc <-sample$Parc[1]
    res$xp<- sample$xp[1]
    res$yp <- sample$yp[1]
    return(res[,c("Type","Parc","xp","yp","Variable","Estimacion")])
  }else{
    return(res)
  }
  
}

n_estimaciones<-function(sample,lado,rotate=FALSE){
  map_dfr(group_split(sample,Parc,Type),estimacion,lado=lado,rotate=rotate)
}

pop_plot <- function(forest_data,lado){
  rect <- data.frame(x=c(0,0,lado,lado,0),y=c(0,lado,lado,0,0))
  ggplot(forest_data) +
    geom_polygon(data=rect,aes(x=x,y=y),col="red",fill="darkgreen",alpha=0.1)+
    geom_circle(aes(x0=x,y0=y,r=dn_cm/20),col="black",fill="burlywood4")+
    xlim(c(-20,lado+20)) + ylim(c(-20,lado+20))+
    coord_fixed(ratio = 1) +
    labs(x="x(m)",y="y(m)") +
    theme_bw(base_size = 16) +
    theme(axis.title = element_blank())
}



plot_selection <- function(p,trees,tree_center=TRUE,all=FALSE,add_hd=FALSE){
  
  
  type <- trees$Type[1]
  title <- switch(type,
                  r_fijo = "Radio fijo 15m",
                  r_variable = "R anidados d<15cm 10m, d>=15cm 20m",
                  r_relascopio = "Relascopio BAF=1"
  )
  
  if(add_hd){
    p <- p +
      geom_label(aes(x=x,y=y-3,label=paste("d: ",dn_cm)),size=4,fill="darkgreen",alpha=0.3)+
      geom_label(aes(x=x,y=y-8,label=paste("h: ",dn_cm)),size=4,fill="blue",alpha=0.3)
  }
  
  if(all){
    if(type == "r_relascopio"){
      p <- p + geom_circle(aes(x0=x,y0=y,r=.data[[type]],fill=.data[[type]]),alpha=0.4)
    }else{
      p <- p + geom_circle(aes(x0=x,y0=y,r=.data[[type]],fill=factor(.data[[type]])),alpha=0.4)
    }
    
  }
  
  if(!is.na(trees$dn_cm[1])){
    if(tree_center){
      p <- p  + geom_circle(data=trees,aes(x0=x,y0=y,r=radio_sel_m),fill="purple",alpha=0.2)
    }else{
      trees2 <- trees |> group_by(radio_sel_m) |> filter(row_number()==1) |> ungroup()
      p <- p  + geom_circle(data=trees2,aes(x0=xp,y0=yp,r=radio_sel_m),fill="purple",alpha=0.2)  
    }
    p <- p + geom_circle(data=trees,aes(x0=x,y0=y,r=dn_cm/20),col="green",fill="darkgreen",lwd=0.5)
  }
  p <- p + geom_point(data=trees[1,],aes(x=xp,y=yp),shape=13,col="red",size=4)
  p <- p + guides(fill=FALSE)+ggtitle(title)
  p
}

plot_n_selections <- function(p,trees,tree_center=TRUE,all=FALSE){
  
  trees$Parc <- as.factor(trees$Parc)

  type <- trees$Type[1]
  points <- trees |> group_by(Parc, Rep) |> filter(row_number()==1) |> ungroup()
  title <- switch(type,
                  r_fijo = "Radio fijo 15m",
                  r_variable = "R anidados d<15cm 10m, d>=15cm 20m",
                  r_relascopio = "Relascopio BAF=1"
  )
  if(all){
    p <- p + geom_circle(aes(x0=x,y0=y,r=radio_sel_m), fill="grey50",alpha=0.2)
  }

  if(!all(is.na(trees$dn_cm))){
    if(tree_center){
      p <- p  + geom_circle(data=trees,aes(x0=x,y0=y,r=radio_sel_m,fill=Parc),alpha=0.2)
    }else{
      trees2 <- trees |> group_by(Parc, radio_sel_m) |> filter(row_number()==1) |> ungroup()
      p <- p  + geom_circle(data=trees2,aes(x0=xp,y0=yp,r=radio_sel_m,fill=Parc),alpha=0.2)  
    }
    p <- p + geom_circle(data=trees,aes(x0=x,y0=y,r=dn_cm/20,fill=Parc))
  }
  p <- p + geom_point(data=points,aes(x=xp,y=yp,col=Parc),shape=13,size=4)
  p <- p + guides(fill=FALSE,color=FALSE)+ggtitle(title)
  p
}

prepare_long1 <- function(data){

  data_long <- pivot_longer(data[,c("Parc","N","G","V","h_media","dg","ho")],
                            cols = c("N","G","V","h_media","dg","ho"),
                            names_to = "parametro",values_to = "estimacion")
  means <- data_long|> group_by(parametro)|> summarise_all(mean,na.rm=TRUE)
  
  means$type_est <- "n-parcelas"
  data_long$type_est <- "1-parcela"
  
  all <- rbind(means,data_long)
  all$type_est <- factor(all$type_est,levels=c("1-parcela","n-parcelas"),ordered=TRUE)
  
  variation <- data_long|> group_by(type_est,parametro)|> 
    summarise(mean=mean(estimacion,na.rm=TRUE),sd=sd(estimacion,na.rm=TRUE),.groups = "keep") |>
    transmute(xmin = mean - 2*sd,xmax=mean + 2*sd) |> ungroup()
  
  variation2 <- variation
  variation<- rbind(variation,variation2)
  variation$type_est <- factor(variation$type_est,levels=c("1-parcela","n-parcelas"),ordered=TRUE)
  return(list(all=all,variation=variation))
}


add_samples_plot<-function(p_int,first,variation){
  
  p_int <- p_int[p_int$parametro%in%c("N","G","V","h_media","dg","ho"),]
  p_int2 <- p_int
  p_int$type_est <- "1-parcela"
  p_int2$type_est <- "n-parcelas"
  p_int <- rbind(p_int,p_int2)
  p_int$type_est <- factor(p_int$type_est,levels=c("1-parcela","n-parcelas"),ordered=TRUE)

  # print(first)
  to_plot <- prepare_long1(first)
  
  variation <- variation[variation$Type==first$Type[1],]
  variation$x_min <- variation$mean-3*variation$sd
  variation$x_max <- variation$mean+3*variation$sd
  variation$type_est <- "1-parcela"
  
  variation2 <- variation
  variation2$type_est <- "n-parcelas"
  
  variation <- rbind(variation,variation2)
  variation$type_est <- factor(variation$type_est,levels=c("1-parcela","n-parcelas"),ordered=TRUE)

  to_plot$all$y <- ifelse(to_plot$all$type_est=="1-parcela",1,0)
  
  print(to_plot)
  p <- ggplot(variation) +
    facet_wrap(.~parametro,scales="free")+
    geom_point(data=to_plot$all,aes(x=estimacion,y=y,col=type_est,fill=type_est),shape=20,size=4)
  
  if(max(to_plot$all$Parc)>1){
    densities <- to_plot$all[to_plot$all$type_est=="1-parcela",]
    densities <- group_split(densities,parametro)
    print(densities)
    
    densities <- map_dfr(densities,function(x){
      
      x<-x[!is.na(x$estimacion),]
      if(dim(x)[1]<2){
        return(NULL)
      }else{
        dens <- density(x$estimacion)
        dens <- data.frame(x=c(dens$x,dens$x[1]),y=c(dens$y,dens$y[1]),
                           parametro=x$parametro[1])
        dens$y <- 0.5+0.4*dens$y/max(dens$y)
        return(dens)
      }
    })
    
    if(dim(densities)[1]>0){
      print("dens")
      print(densities)
      p <- p + geom_polygon(data=densities,aes(x=x,y=y),col="black",fill="grey",linewidth=0.5,alpha=0.5)
    }
    
  }
  

  
  p + geom_linerange(data=variation[variation$type_est=="1-parcela",],aes(y=0,xmin=x_min,xmax=x_max),col="red")+
    geom_vline(data=p_int,aes(xintercept=Valor),col="black")+
    scale_fill_manual(values=c("1-parcela"="red","n-parcelas"="blue"))+
    scale_color_manual(values=c("1-parcela"="red","n-parcelas"="blue")) +
    guides(fill=NULL,color=NULL)+
    theme(legend.position = "bottom")+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
  
}


prepare_long_n <- function(data){

  data_long <- pivot_longer(data[,c("Rep","Parc","N","G","V","h_media","dg","ho")],
                            cols = c("N","G","V","h_media","dg","ho"),
                            names_to = "parametro",values_to = "estimacion")
  
  means <- data_long|> group_by(Rep,parametro)|> summarise_all(mean,na.rm=TRUE)
  means$y <- 0.5
  means$type_est <- "n-parcelas"
  data_long$type_est <- "1-parcela"
  data_long$y <- 0.75
  
  variation <- data_long|> group_by(type_est,parametro)|> 
    summarise(mean=mean(estimacion,na.rm=TRUE),sd=sd(estimacion,na.rm=TRUE),.groups = "keep") |>
    transmute(xmin = mean - 2*sd,xmax=mean + 2*sd) |> ungroup()
  variation$y <- 1
  variation2 <- means|> group_by(type_est,parametro)|> 
    summarise(mean=mean(estimacion),sd=sd(estimacion),.groups = "keep") |>
    transmute(xmin = mean - 2*sd,xmax=mean + 2*sd) |> ungroup()
  variation2$y<-0.25
  
  all <- rbind(means,data_long)
  all$type_est <- factor(all$type_est,levels=c("1-parcela","n-parcelas"),ordered=TRUE)
  
  variation<- rbind(variation,variation2)
  variation$type_est <- factor(variation$type_est,levels=c("1-parcela","n-parcelas"),ordered=TRUE)

  return(list(all=all,variation=variation))
}

add_samples_n_plots<-function(p_int,all,variation,n){
  
  p_int <- p_int[p_int$parametro%in%c("N","G","V","h_media","dg","ho"),]
  p_int2 <- p_int
  p_int$type_est <- "1-parcela"
  p_int2$type_est <- "n-parcelas"
  p_int <- rbind(p_int,p_int2)
  p_int$type_est <- factor(p_int$type_est,levels=c("1-parcela","n-parcelas"),ordered=TRUE)
  

  
  variation <- variation[variation$Type==all$Type[1],]
  variation$x_min <- variation$mean-2*variation$sd
  variation$x_max <- variation$mean+2*variation$sd
  variation$type_est <- "1-parcela"
  
  variation2 <- variation
  variation2$type_est <- "n-parcelas"
  variation2$x_min <- variation$mean-2*variation$sd/sqrt(n)
  variation2$x_max <- variation$mean+2*variation$sd/sqrt(n)
  
  variation <- rbind(variation,variation2)
  variation$type_est <- factor(variation$type_est,levels=c("1-parcela","n-parcelas"),ordered=TRUE)
  
  to_plot <- prepare_long_n(all)
  
  to_plot$all <- ungroup(to_plot$all)
  print("to_plot")
  print(to_plot)
  
  ggplot(variation) +
    facet_grid(type_est~parametro,scales="free_x")+
    geom_point(data=to_plot$all,aes(x=estimacion,y=0.25,col=type_est,fill=type_est),shape=20,size=4)+
    geom_density(data=to_plot$all,aes(x=estimacion,fill=type_est,col=type_est),alpha=0.4) +
    geom_linerange(data=variation,aes(y=0.75,xmin=x_min,xmax=x_max,col=type_est))+
    geom_vline(data=p_int,aes(xintercept=Valor),col="black")+
    scale_fill_manual(values=c("1-parcela"="red","n-parcelas"="blue"))+
    scale_color_manual(values=c("1-parcela"="red","n-parcelas"="blue")) +
    guides(fill=NULL,color=NULL)+
    theme(legend.position = "bottom")+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())

  
}

normal_approx <- function(estimates,p_int,n,type,variation,K){
  
  reps <- 100
  last <- reps*n
  colnames(p_int)[2]<-"target"
  p_int <- p_int[p_int$parametro%in%variation$parametro,]
  print(p_int)
  variation <- variation[variation$Type==type,]
  variation$x_min <- variation$mean-2*variation$sd
  variation$x_max <- variation$mean+2*variation$sd
  variation$x_min2 <- variation$mean-2*variation$sd
  variation$x_max2 <- variation$mean+2*variation$sd
  
  estimates<-estimates[estimates$Type==type,]
  
  positions<-sample(1:dim(estimates)[1],last,replace=TRUE)
  estimates <- estimates[positions,]
  estimates$Parc <- rep(1:n,times=reps)
  estimates$Rep <- rep(1:reps,each=n)
  
  limits <- pivot_longer(variation[,c("x_min","x_max","parametro")],
                                       cols = c("x_min","x_max"),
                                       names_to = "type",values_to = "value")
  
  
  estimates <-  pivot_longer(estimates[,c("Rep","Parc","N","G","V","h_media","dg","ho")],
                             cols = c("N","G","V","h_media","dg","ho"),
                             names_to = "parametro",values_to = "estimacion")
  estimates <- group_by(estimates,parametro,Rep)|>summarize(estimacion=mean(estimacion,na.rm=TRUE))|>ungroup()

 
  estimates<-merge(estimates,variation,by="parametro")
  estimates$sd_n <- estimates$sd/sqrt(n)
  estimates <- merge(estimates,p_int,by="parametro")

  print("hola")
  print(estimates)
  norm <- estimates %>% 
    group_by(parametro) %>% 
    reframe(x=seq(min(target,na.rm=TRUE)-3*mean(sd),max(target,na.rm=TRUE)+3*mean(sd),length.out=200),
            y = dnorm(x, mean = mean(target,na.rm=TRUE), sd = mean(sd_n,na.rm=TRUE) )) 
  print(norm)

  ggplot(data=norm) +
    facet_wrap(.~ parametro,scales="free") +
    geom_line(data=norm,aes(x=x,y = y),col="blue") + 
    geom_point(data=estimates,aes(x=estimacion,y=0),shape=20,col="red")+
    geom_density(data=estimates,aes(x=estimacion),fill="red",colour = "red",alpha=0.2)+
    geom_vline(data=p_int,aes(xintercept=target),colour = "red")+
    geom_point(data=limits,aes(x=value,y=0),colour = "red",alpha=0)+
    ggtitle("Aproximación a una distribución normal al aumentar n (100 repeticiones)")+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
  
}

standard_dev<- function(var,n){
  a<-data.frame(n=1:50,id=1)
  var <- merge(var,a)
  var$sd_n <- var$sd/sqrt(var$n)
  var$color <- ifelse(var$n==n,"red","black")
  
  red <- var[var$color=="red",]
  
  ggplot(var,aes(x=n,y=sd_n))+facet_wrap(.~parametro,scales="free_y")+
    geom_point(aes(color=color))+geom_path() + 
    geom_point(data=red,aes(color=color),pch=20,size=3)+
    xlab("Desiviación tipica del estimador final")+
    scale_color_manual(values=c("red"="red","black"="black"))+
    guides(color="none")+
    ggtitle("Cambio en la desviación típica al aumentar n")
}

standard_dev2<- function(var,n,samples=NULL){
  
  a<-data.frame(n=1:50,id=1)
  var$id <- 1
  var <- merge(var,a)
  
  var$sd_n <- var$sd/sqrt(var$n)
  var$color <- ifelse(var$n==n,"red","black")
  
  red <- var[var$color=="red",]

  
  p <- ggplot(var,aes(x=n,y=sd_n))+facet_wrap(.~parametro,scales="free_y")+
    geom_point(aes(y=1.5*sd_n),alpha=0)+
    geom_point(aes(color=color))+geom_path() + 
    geom_point(data=red,aes(color=color),pch=20,size=3)
  if(!is.null(samples)){
    samples <- samples |> group_by(Rep,parametro)|> summarise(sd_n=sd(estimacion)/sqrt(n())) |> ungroup()
    samples$n <- n
    p <- p + geom_point(data=samples,color="blue",alpha=0.5,pch=20,size=3)
  }
  p <- p + 
    xlab("Desiviación tipica del estimador final")+
    scale_color_manual(values=c("red"="red","black"="black"))+
    guides(color="none")+
    ggtitle("Cambio en la desviación típica al aumentar n")
  p
}

get_estimatesIC <- function(estimates,type,n,K,reps){
  
  estimates <- estimates |> filter(Type==type)
  print(estimates)
  estimates <- estimates[sample(1:K,n*reps,replace=TRUE),]

  estimates$Parc <- rep(1:n,reps)
  estimates$Rep <- rep(1:reps,each=n)
  
  estimates <- pivot_longer(estimates[,c("Rep","Parc","N","G","V","h_media","dg","ho")],
                            cols = c("N","G","V","h_media","dg","ho"),
                            names_to = "parametro",values_to = "estimacion")
  
  return(estimates)
}

get_pilot <- function(estimates,type,n_pilot,K,wide=FALSE){
  
  estimates <- estimates |> filter(Type==type)
  estimates <- estimates[sample(1:K,n_pilot,replace=TRUE),]
  estimates$Parc <- 1:n_pilot
  estimates$Rep <- 1
  
  if(wide){
    return(estimates)
  }else{
    estimates <- pivot_longer(estimates[,c("Rep","Parc","N","G","V","h_media","dg","ho")],
                              cols = c("N","G","V","h_media","dg","ho"),
                              names_to = "parametro",values_to = "estimacion")
    
    return(estimates)
  }
  
}




confint_plot<-function(estimates, var, par_int,conf){
  
  par_int <- merge(par_int, var, by="parametro")
  
  var$x_min <- var$mean - 5 * var$sd
  var$x_max <- var$mean + 5 * var$sd
  
  reps <- max(estimates$Rep)
  
  # print(estimates)
  estimates <- estimates[!is.na(estimates$estimacion),]
  means <- estimates |> group_by(Rep,parametro) |> summarize(mean=mean(estimacion,na.rm=TRUE),
                                                   sd=sd(estimacion,na.rm=TRUE)/sqrt(n()),n=n()) |> ungroup()
  conf <- -qnorm((1-conf)/2)
  means$x_min <- means$mean - conf*means$sd
  means$x_max <- means$mean + conf*means$sd
  
  means$rel_error <- round(100*conf*means$sd/means$mean,1)
  means$label_sd1 <- paste("hat( 'sd' )( hat( mu )[p] )==",round(means$sd*sqrt(means$n),1),sep="")
  means$label_sdn<- paste("hat( 'sd' )( hat( mu )[final] )==",round(means$sd,1),sep="")
  means$label_rel_error<- paste("hat( epsilon )[(1-alpha)]==",round(means$rel_error,1),"~'%'",sep="")
  
  labels <- merge(par_int,means,by="parametro")
  labels$pos1 <- labels$mean.x - 2*labels$sd.x
  labels$pos2 <- labels$mean.x + 2*labels$sd.x
    # par_int <- filter(par_int, parametro%in%var$parametro)

  line_range <- merge(var,expand.grid(parametro=unique(var$parametro),Rep=1:reps),by="parametro")
  
  # print(par_int)
  # print(labels)
  print(means)
  p <- ggplot(var) + facet_wrap(.~parametro,scales="free_x") + 
    geom_point(aes(x=x_min,y=reps),alpha=0)+
    geom_point(aes(x=x_max,y=1),alpha=0)+
    geom_linerange(data=line_range,aes(xmin=x_min,xmax=x_max,y=Rep-0.6),col="grey30",linetype=2)+
    
    geom_point(data=estimates,aes(x=estimacion,y=Rep),shape=20,alpha=0.5,col="red")+
    
    geom_point(data=means,aes(x=mean,y=Rep-0.4),shape=20,alpha=0.5,col="blue",size=3)+
    geom_linerange(data=means,aes(y=Rep-0.4,xmin=x_min,xmax=x_max),col="blue")+
    geom_vline(data=par_int,aes(xintercept=Valor),col="black")
    
  if(max(estimates$Parc)>1){
    p <- p + geom_text(data=labels,aes(x=pos2,y=Rep-0.4,label=label_rel_error),parse=TRUE,hjust=0,size=4.5,col="blue") +
             geom_text(data=labels,aes(x=pos1,y=Rep-0.4,label=label_sdn),parse=TRUE,hjust=1,size=4.5,col="blue") +
             geom_text(data=labels,aes(x=pos1,y=Rep,label=label_sd1),parse=TRUE,hjust=1,size=4.5,col = "red")
      
  }
    
    p <- p+ guides(fill=NULL,color=NULL)+theme(legend.position = "bottom") + xlab("Estimación")
            theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                  element_text(family="mono"))
    p
  
  
}


prepare_error_pol <- function(means_sd,conf_level=0.95){
  
  n <- max(means_sd$n)
  
  variation_n <- expand.grid(parametro=unique(means_sd$parametro),n_samp=n:50,side=1:2)
  variation_n <- merge(means_sd,variation_n,by="parametro")

  variation_n$q <- qt((1-conf_level)/2,n-1)
  
  variation_n$bound <- ifelse(variation_n$side==1,variation_n$q*variation_n$sd,
                              -variation_n$q*variation_n$sd)
  variation_n$bound <- variation_n$mean + variation_n$bound/sqrt(variation_n$n_samp)
  
  v1 <- filter(variation_n,side==1) |> group_by(parametro) |> 
    arrange(desc(n_samp)) |> mutate(order=n:50)|> ungroup()
  
  
  v2 <- filter(variation_n,side==2)|> group_by(parametro) |> 
    arrange(n_samp) |> mutate(order=(51:(100-(n-1))))|> ungroup()
  
  v3 <- v1[v1$order==1,]
  v3$order <- 100-(n-2)
  res <- rbind(v1,v2,v3)
  
  return(res[order(res$parametro,res$order),])
}

sample_alloc_plot <- function(piloto,conf_level=0.95,max_rel_error=0.1,current_n){
  
  piloto<-pivot_longer(piloto[,c("Rep","Parc","N","G","V","h_media","dg","ho")],
               cols = c("N","G","V","h_media","dg","ho"),
               names_to = "parametro",values_to = "estimacion")
  
  means_sd <- group_by(piloto,parametro)|> 
    summarize(mean=mean(estimacion,na.rm=TRUE),
              sd=sd(estimacion,na.rm=TRUE),n=n())|>
    ungroup()
  
  max_y <- max(means_sd$n)
  
  means_sd$sd_n <- means_sd$sd/sqrt(means_sd$n)
  means_sd$q <- NA
  means_sd$q2 <- NA
  

  current_n <- ifelse(current_n<2,2,current_n)
  
  means_sd[["q"]]<- ifelse(means_sd$n<2,NA,qt(p = (1-conf_level)/2,df=means_sd$n-1))
  means_sd[["q2"]] <- ifelse(current_n<2,NA,qt(p=(1-conf_level)/2,df=current_n-1))
  
  means_sd$error <- (means_sd$q*means_sd$sd_n)
  means_sd$rel_error <- means_sd$error/means_sd$mean
  
  means_sd$error_curr <- means_sd$q*(means_sd$sd/sqrt(current_n))
  means_sd$rel_error_curr <- means_sd$error_curr/means_sd$mean
  
  means_sd$bound_min <- means_sd$mean*(1-max_rel_error)
  means_sd$bound_max <- means_sd$mean*(1+max_rel_error)
  
  means_sd$bound_min_curr <- means_sd$mean*(1-means_sd$rel_error_curr)
  means_sd$bound_max_curr <- means_sd$mean*(1+means_sd$rel_error_curr)

  variation_n <- prepare_error_pol(means_sd,conf_level)

  
  p<-ggplot(variation_n) +
    facet_wrap(.~parametro,scales="free_x") +
    geom_rect(data=means_sd,aes(xmin=bound_min,xmax=bound_max,ymin=0,ymax=50),
              col="black",fill="grey20",alpha=0.2)+
    geom_vline(data=means_sd,aes(xintercept=mean),col="black",linetype=2)+
    
    geom_polygon(aes(x=bound,y=n_samp),col="blue",fill="purple",alpha=0.1)+
    
    geom_point(data=piloto,aes(x=estimacion),y=max_y+0.2,col="red",shape=20,alpha=0.8)+
    
    geom_point(data=means_sd,aes(x=mean),y=max_y,col="blue",shape=20,alpha=0.5,size=3)+
    geom_linerange(data=means_sd,aes(xmin=(mean-error),xmax=(mean+error)),y=max_y,col="blue",alpha=0.5)
  
    if(current_n>=max_y){
      p <- p + geom_point(data=means_sd,aes(x=mean),y=current_n,col="darkgreen",shape=20,alpha=0.5,size=3)+
              geom_linerange(data=means_sd,aes(xmin=bound_min_curr,xmax=bound_max_curr),y=current_n,col="darkgreen",alpha=0.5)
    }
    
    p <- p +  ylim(0,50)
    p
 
}
