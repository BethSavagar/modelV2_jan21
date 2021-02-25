#Proportion of die-out

kvec <- c(0.1, 0.3, 0.5, 0.7, 1, Inf)
#kvec <- c(0.1, 0.3, 0.5, 0.7, 1, Inf)
MasterList <- list(k0.1List,
                   k0.3List,
                   k0.5List,
                   k0.7List,
                   k1List,
                   kInfList)
#MasterList <- list(k0.1List,k0.3List,k0.5List,k0.7List,k1List, kInfList)

propDieOut <- function(matrix){
  end <- matrix[,ncol(matrix)]
  DO <- sum(end == 0)
  pDO <- DO/length(end)
  return(pDO)
}

pDieOut <- vector()

for(k in 1:length(kvec)){
  
  kval <- kvec[k]
  setList <- MasterList[[k]]
  
  pDO <- sapply(setList, propDieOut)
  
  kDO<- cbind(kval,pars,pDO)
  
  pDieOut <- rbind(pDieOut, kDO)
}


ggplot(data = pDieOut, aes(x = factor(R0), y = pDO))+
  geom_point(aes(group = kval, colour = factor(kval)), alpha = 1
  )+
  scale_colour_discrete("k value")+
  scale_color_brewer(palette="Paired", name = "k value")+
  labs(x = "R0", y = "proportion die-out (1000G)")+
  facet_grid(kval~Rperiod, labeller = label_both)+
  theme_bw()

ggplot(data = pDieOut, aes(x = factor(R0), y = pDO))+
  geom_point(aes(group = Rperiod, colour = factor(Rperiod)), alpha = 1
  )+
  scale_colour_discrete("Rp value")+
  scale_color_brewer(palette="Paired", name = "Rp value")+
  labs(x = "R0", y = "proportion die-out (1000G)")+
  facet_grid(Rperiod~kval, labeller = label_both)+
  theme_bw()

ggplot(data = pDieOut, aes(x = factor(kval), y = pDO))+
  geom_point(aes(group = R0, colour = factor(R0)), alpha = 1
  )+
  scale_colour_discrete("R0 value")+
  scale_color_brewer(palette="Paired", name = "R0 value")+
  labs(x = "k val", y = "proportion die-out (1000G)")+
  facet_grid(R0~Rperiod, labeller = label_both)+
  theme_bw()

ggplot(data = pDieOut, aes(x = factor(kval), y = pDO))+
  geom_line(aes(group = R0, colour = factor(R0)), alpha = 1
  )+
  scale_colour_discrete("R0 value")+
  scale_color_brewer(palette="Paired", name = "R0 value")+
  labs(x = "k val", y = "proportion die-out (1000G)")+
  facet_wrap(~Rperiod, labeller = label_both)+
  theme_bw()


ggplot(data = pDieOut, aes(x = factor(R0), y = pDO))+
  geom_line(aes(group = kval, colour = factor(kval)), alpha = 1
  )+
  scale_colour_discrete("k value")+
  scale_color_brewer(palette="Paired", name = "k value")+
  labs(x = "R0", y = "proportion die-out (1000G)")+
  facet_wrap(~Rperiod, labeller = label_both)+
  theme_bw()



# CONTOUR PLOT
ggplot(pDieOut, aes(x = as.factor(Rperiod), y = R0, fill = pDO))+
  geom_tile()+
  facet_wrap(~kval)+
  scale_fill_distiller(palette = "Blues", name = "proportion Die Out")+
  xlab("Refractory Period")+
  scale_y_continuous(breaks = c(seq(0,5,0.5)),
                     labels =  c(as.character(seq(0,5,0.5))))+
  theme_bw()

ggplot(pDieOut, aes(x = as.factor(kval), y = R0, fill = pDO))+
  geom_tile()+
  facet_wrap(~Rperiod, labeller = label_both)+
  scale_fill_distiller(palette = "Blues", direction = -1, name = "proportion Die Out")+
  xlab("kval")+
  scale_y_continuous(breaks = c(seq(0,5,0.5)),
                     labels =  c(as.character(seq(0,5,0.5))))+
  theme_bw()
