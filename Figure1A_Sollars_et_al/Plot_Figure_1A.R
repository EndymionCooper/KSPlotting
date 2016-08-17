library(ggplot2)

M1 <- read.table("Coffee_KS.txt", header=FALSE)
M2 <- read.table("MonkeyFlower_KS.txt", header=FALSE)
M3 <- read.table("Olive_KS.txt", header=FALSE)
M4 <- read.table("Tomato_KS.txt", header=FALSE)
M5 <- read.table("Grape_KS.txt", header=FALSE)
M6 <- read.table("Bladderwort_KS.txt", header=FALSE)
M7 <- read.table("Ash_KS.txt", header=FALSE)

ggplot() +
  ggtitle("") +
  xlab("Substitutions per synonymous site (KS)") +
  ylab("Density") +
  coord_cartesian() +
  scale_x_continuous(breaks=seq(0,3,1/4), labels=c(0,"","","",1,"","","",2,"","","",3)) +
  layer(
    mapping=aes(x=M1$V1[M1$V1 <= 3],y=..density..,colour="Coffee",linetype="Coffee"),
    stat="bin", stat_params=list(binwidth=0.05,drop=TRUE),
    geom="line", geom_params=list()) +
  layer(
    mapping=aes(x=M2$V1[M2$V1 <= 3],y=..density..,colour="Monkey flower",linetype="Monkey flower"),
    stat="bin", stat_params=list(binwidth=0.05,drop=TRUE),
    geom="line", geom_params=list()) +
  layer(
    mapping=aes(x=M3$V1[M3$V1 <= 3],y=..density..,colour="Olive",linetype="Olive"),
    stat="bin", stat_params=list(binwidth=0.05,drop=TRUE),
    geom="line", geom_params=list()) +
  layer(
    mapping=aes(x=M4$V1[M4$V1 <= 3],y=..density..,colour="Tomato",linetype="Tomato"),
    stat="bin", stat_params=list(binwidth=0.05,drop=TRUE),
    geom="line", geom_params=list()) +
  layer(
    mapping=aes(x=M6$V1[M6$V1 <= 3],y=..density..,colour="Bladderwort",linetype="Bladderwort"),
    stat="bin", stat_params=list(binwidth=0.05,drop=TRUE),
    geom="line", geom_params=list()) +
  layer(
    mapping=aes(x=M7$V1[M7$V1 <= 3],y=..density..,colour="Ash",linetype="Ash"),
    stat="bin", stat_params=list(binwidth=0.05,drop=TRUE),
    geom="line", geom_params=list()) +
  layer(
    mapping=aes(x=M5$V1[M5$V1 <= 3],y=..density..,colour="Grape",linetype="Grape"),
    stat="bin", stat_params=list(binwidth=0.05,drop=TRUE),
    geom="line", geom_params=list()) +
  theme(
    text = element_text(family="Helvetica"),
    legend.text = element_text(size=12),
    axis.title = element_text(size=12),
    axis.text = element_text(size=10),
    panel.background= element_rect(fill = "white",colour = "black"),
    panel.grid.major= element_blank(),
    panel.grid.minor= element_blank(),
    panel.background = element_rect(colour = "white"),
    legend.key = element_rect(fill ="white"),
    legend.position = c(0.8,0.6),
    legend.background = element_rect(colour ="black")) +
  scale_colour_manual(guide=guide_legend(title = NULL),values=c("Coffee"= "brown",
"Monkey flower"= "yellow",
"Olive"= "black",
"Tomato"= "red",
"Bladderwort"= "black",
"Ash"="blue",
"Grape"="darkmagenta")) +
  
  scale_linetype_manual(guide=guide_legend(title = NULL),values=c("Coffee"= "solid",
"Monkey flower"= "solid",
"Olive"= "dashed",
"Tomato"= "solid",
"Bladderwort"= "solid",
"Ash"="solid",
"Grape"="solid"))
  
ggsave("kS.png", width = 15, height = 10, units = "cm")
