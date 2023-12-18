setwd("~/PostDoc Computer/Daphnia/Pasteuria")

library(dplyr)
library(tidyr)
library(ggplot2)
library(poppr)
library(pegas)
library("phangorn")
library("ape")
library(adegenet)
library(usedist)
library(cowplot)
library(ggtree)
library(hierfstat)
library(reshape2)
library(cluster)
library(vegan)
library(multivariance)

#Look at 2015 Data only

data<-read.csv("P.Microsats7March19.csv", header=TRUE, stringsAsFactors = FALSE)
data.coinf<-read.csv("P.Microsats7March19.coinf.csv", header=TRUE, stringsAsFactors = FALSE)
lakedistances<-read.table("LakeDistancesLong.txt", sep="", header=TRUE)

data2015<-filter(data, year==2015)

#dropping p17 - same across all samples
data2015.2<-select(data2015, sample, Species, Lake, year, collection.date, Julian, p1, p2, p3, p4, p7, p11, p12, p16, p18, p19)
data2015.3<-data2015.2[rowSums(is.na(data2015.2)) <= 2, ] # include samples that amplified at 8 or more out of 10 loci
data2015.3$species_lake_date<- paste(data2015.3$Species,data2015.3$Lake,data2015.3$collection.date,sep="_")
data2015.3$species_lake<- paste(data2015.3$Species,data2015.3$Lake,sep="_")

#Fix errors and drop ambiguous samples
data2015.3<-filter(data2015.3, sample!="B.R1.10.28.2015(6)")
data2015.3<-filter(data2015.3, sample!="B.R1.8.4.15(B.R1.10.28.15.1)")
data2015.3<-filter(data2015.3, sample!="B.R2.8.4.15(B.R2.10.28.15)")
data2015.3<-filter(data2015.3, sample!="B.R2.10.28.2015")
data2015.3$sample[data2015.3$sample=="B.D5.8.4.2015(B.D5.10.28.15)"]<-"B.D5.8.4.15"
data2015.3$sample[data2015.3$sample=="B.R5.8.4.15(B.R5.10.28.2015)"]<-"B.R5.8.4.15"

info<-unique(data2015.3[,c("species_lake_date","species_lake","Species","Lake","collection.date","Julian")])

#get into format for poppr
dat<-select(data2015.3, sample, species_lake_date, p1, p2, p3, p4, p7,p11,p12, p16,p18,p19 )
#dat<-select(data2015.3, sample, species_lake, p1, p2, p3, p4, p7,p11,p12, p16,p18,p19 )
write.table(dat, "Past2015.txt", col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE)

#Make datasets of denteifera and retrocurva hosts
dentifera<-filter(data2015.3, Species=="dentifera")
dent.data<-select(dentifera, sample, species_lake_date,p1, p2, p3, p4, p7,p11,p12, p16,p18,p19)
dentgroups<-dent.data%>%group_by(species_lake_date)%>%summarise(count=n())
dent.data<-filter(dent.data, species_lake_date!="dentifera_B_8.17.2015") #only one of these, so don't include in fst analysis

retrocurva<-filter(data2015.3, Species=="retrocurva")
retro.data<-select(retrocurva, sample, species_lake_date,p1, p2, p3, p4, p7,p11,p12, p16,p18,p19)
retrogroups<-retro.data%>%group_by(species_lake_date)%>%summarise(count=n())
retro.data<-filter(retro.data, species_lake_date!="retrocurva_B_8.4.2015") #only one of these, so don't include in fst analysis
retro.data<-filter(retro.data, species_lake_date!="retrocurva_CW_10.9.2015") #only one of these, so don't include in fst analysis
retro.data<-filter(retro.data, species_lake_date!="retrocurva_M_9.8.2015") #only one of these, so don't include in fst analysis

write.table(dent.data, "dent.data.txt", col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE)
write.table(retro.data, "retro.data.txt", col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE)

z <- read.loci("Past2015.txt", col.loci=3:12,row.names=1, col.pop=2, header=TRUE, allele.sep="-")
dent<-read.loci("dent.data.txt", col.loci=3:12,row.names=1, col.pop=2, header=TRUE, allele.sep="-")
retro<-read.loci("retro.data.txt", col.loci=3:12,row.names=1, col.pop=2, header=TRUE, allele.sep="-")

Z<-loci2genind(z, ploidy = 1)
DENT<-loci2genind(dent, ploidy = 1)
RETRO<-loci2genind(retro, ploidy = 1)

#######################################################################################
#Prepare data for tree and AMOVA
Y<-genind2genpop(Z)

genind2genalex(Z, filename = "Past2015genalex.csv", quiet = FALSE, pop = NULL,
               allstrata = TRUE, geo = FALSE, geodf = "xy", sep = ",",
               sequence = FALSE, overwrite=TRUE)

Past<-read.genalex("Past2015genalex.csv")

splitStrata(Past)<-~species/lake/date
setPop(Past)<-~species/lake/date

#splitStrata(Past)<-~species/lake
#setPop(Past)<-~species/lake

#all sites are polymorphic
iPast <- informloci(Past)
info_table(Past, plot = TRUE)
#This tells me that I shouldn't drop any more loci or even samples.

tab(Past) #look at the multilocus genotype table
nmll(Past) #count the number of multilocus genotypes 
mll(Past)#show the multilocus genotype definitions

PastDistance<-provesti.dist(Past) #This uses contracted MLGs.

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbP2<- c("hotpink3", "dodgerblue4", "darkorange2", "darkseagreen4")
HOST <- factor(data2015.3$Species)
LAKE <- factor(data2015.3$Lake)
library(RColorBrewer)
Color.Palette<-brewer.pal(4, "Set1")
Color.Pal<-brewer.pal(8, "Paired")
Color.Pal<-c("navyblue", "#1F78B4", "darkseagreen4", "#33A02C","deeppink3", "#E31A1C", "#FDBF6F", "#FF7F00")

tree<-hclust(PastDistance, method="average")
as.dendrogram(tree)
plot(tree)
tree<-aboot(Past, strata = NULL, tree = "upgma", distance = "prevosti.dist",
      sample = 100, cutoff = 30, showtree = TRUE, missing = "ignore",
      mcutoff = 0.5, quiet = FALSE, root = NULL)
      
tree<-as.phylo(tree)
byHost<-ggtree(tree, size=1, layout="rectangular")+geom_tiplab(color=Color.Palette[HOST], size=7, offset=0.01)+
  geom_treescale(width=0.05, x=0.1, y=-5, fontsize=5)+ggplot2::xlim(0, 0.7)+
  ggplot2::ylim(-5, 90)+geom_nodelab(size = 6, col= "black", nudge_y = 0.03, nudge_x = -0.01)

L<-as.data.frame(LAKE)
byLake<-ggtree(tree, size=1)%<+% L+geom_tiplab(color=Color.Pal[LAKE], size=7, offset=0.01)+
  geom_treescale(width=0.05, x=0.1, y=-5, fontsize=5)+ggplot2::xlim(0, 0.7)+  
  ggplot2::ylim(-5, 90)+geom_nodelab(size = 6, col= "black", nudge_y = 0.03, nudge_x = -0.01)

Figure3<-plot_grid(byLake, byHost, labels=c("A","B"), label_size=30)
save_plot("Pasteuria2015Figure3.jpg",Figure3,base_height = 20,base_width = 15)

setPop(Past)<-~species/lake/date
#setPop(Past)<-~species/lake

full2015.Amova<-poppr.amova(Past, ~species/lake/date, within=FALSE, filter=TRUE,threshold =0 , missing="ignore",
                            correction="lingoes", cutoff=0)
#full2015.Amova<-poppr.amova(Past, ~species/lake, within=FALSE, filter=TRUE,threshold =0 , missing="ignore",
       #                     correction="lingoes", cutoff=0)

sig.full2015<-randtest(full2015.Amova, nrepet=1000)

#Flip species and lake
setPop(Past)<-~lake/species/date
#setPop(Past)<-~lake/species

full2015.Amova.Flip<-poppr.amova(Past, ~lake/species/date, within=FALSE, filter=TRUE,threshold =0 , missing="ignore",
                                 correction="lingoes", cutoff=0.0)
#full2015.Amova.Flip<-poppr.amova(Past, ~lake/species, within=FALSE, filter=TRUE,threshold =0 , missing="ignore",
 #                                correction="lingoes", cutoff=0.0)
sig.full2015.Flip<-randtest(full2015.Amova.Flip, nrepet=1000)

###############################################################################################################
#Mantel tests
infolake<-select(data2015.3,Lake)
infolake$Lake<-as.factor(infolake$Lake)
distanceLake<-daisy(infolake, metric ="gower") #0 if same; 1 if different
infoSpecies<-select(data2015.3, Species)
infoSpecies$Species<-as.factor(infoSpecies$Species)
distanceSpecies<-daisy(infoSpecies, metric="gower")
infodate<-select(data2015.3, Julian)
distancedate<-daisy(infodate, metric ="manhattan") 

mantel(PastDistance, distanceLake, method="pearson", permutations=999)
mantel(PastDistance, distanceSpecies, method="pearson", permutations=999)
mantel(PastDistance, distancedate, method="pearson", permutations=999)

mantel.partial(PastDistance, distanceLake, distanceSpecies, method="pearson", permutations=999)
mantel.partial(PastDistance, distanceSpecies, distanceLake, method="pearson", permutations=999)
mantel.partial(PastDistance, distanceSpecies, distancedate, method="pearson", permutations=999)
mantel.partial(PastDistance, distanceLake, distancedate, method="pearson", permutations=999)

P.distances<-as.vector(PastDistance)
L.distances<-as.vector(distanceLake)
S.distance<-as.vector(distanceSpecies)
t.distance<-as.vector(distancedate)

Distances<-as.data.frame(cbind(P.distances, L.distances, S.distance, t.distance))
colnames(Distances)<-c("Shared","Lake","Species","time")

p1<-ggplot(Distances, aes(Lake, Shared, group=Lake))+geom_boxplot()+geom_jitter(width=0.25, height=0.025, aes(color=factor(Species)))+
  scale_x_continuous(breaks=c(0,1), labels=c("Same", "Different"))+
  xlab("Lake")+theme_bw()+theme(axis.text = element_text(color="black"))+ylab("Prevosti distance")+
  scale_color_manual(values=c("black","dodgerblue"), labels=c("Same host species", "Different host species"))+
  theme(legend.title = element_blank())+
  theme(legend.position = "bottom")

p2<-ggplot(Distances, aes(Species, Shared, group=Species))+geom_boxplot()+geom_jitter(width=0.25, height=0.025, aes(color=factor(Lake)))+
  scale_x_continuous(breaks=c(0,1), labels=c("Same", "Different"))+
  xlab("Host species")+theme_bw()+theme(axis.text = element_text(color="black"))+ylab("")+
  scale_color_manual(values=c("black","dodgerblue"), labels=c("Same lake", "Different lake"))+
  theme(legend.title = element_blank())+
  theme(legend.position = "bottom")

p3<-ggplot(Distances, aes(time, Shared))+geom_jitter(width=0, height=0.025)+
  xlab("Distance in time (days)")+theme_bw()+theme(axis.text = element_text(color="black"))+ylab("")

Fig4<-plot_grid(p1, p2, p3, labels=c("A","B","C"), nrow = 1, align="h", axis="b")
save_plot("PastPaperFig4.jpg",Fig4, base_width = 12, base_height = 5)
###############################################################################################################
#Compute pairwise Fst
#https://tomjenkins.netlify.app/tutorials/r-popgen-getting-started/

#dentifera Fst
Dent_Fst<-genet.dist(DENT, method="WC84")%>%round(digits=3)
Dent.df <- melt(as.matrix(Dent_Fst), varnames = c("Group1", "Group2"))

lab_order = c("dentifera_B_8.4.2015","dentifera_B_8.31.2015","dentifera_Ce_7.30.2015",
              "dentifera_Ce_11.6.2015","dentifera_CW_10.9.2015",
              "dentifera_G_9.1.2015","dentifera_L_7.21.2015","dentifera_L_9.30.2015", "dentifera_W_8.7.2015")
fst_mat=as.matrix(Dent_Fst)
fst_mat1=fst_mat[lab_order,]
fst_mat2=fst_mat1[,lab_order]
ind = which(upper.tri(fst_mat2), arr.ind = TRUE)
fst.df = data.frame(Group1 = dimnames(fst_mat2)[[2]][ind[,2]],
                    Group2 = dimnames(fst_mat2)[[1]][ind[,1]],
                    Fst = fst_mat2[ ind ])
fst.df$Group1 = factor(fst.df$Group1, levels = unique(fst.df$Group1))
fst.df$Group2 = factor(fst.df$Group2, levels = unique(fst.df$Group2))
fst.df$Fst[fst.df$Fst < 0] = 0
fst.label = expression(italic("F")[ST])

fst.info<-inner_join(fst.df, info, by=c("Group1"="species_lake_date"))
fst.info<-select(fst.info, Group1, Group2, Fst, Lake)
fst.info<-inner_join(fst.info, info, by=c("Group2"="species_lake_date"))
fst.info<-select(fst.info, Group1, Group2, Fst, Lake.x, Lake.y)
fst.info$samelakes<-fst.info$Lake.x==fst.info$Lake.y
fst.info.d<-fst.info
fst.info.d$species<-rep(c("dent"), length(fst.info.d$Fst))

# Extract middle Fst value for gradient argument
mid = max(fst.df$Fst) / 2
fst1<-ggplot(data = fst.df, aes(x = Group1, y = Group2, fill = Fst))+
  geom_tile(colour = "black")+
  geom_text(aes(label = Fst), color="black", size = 5)+
  scale_fill_gradient2(low = "blue", mid = "pink", high = "red", midpoint = mid, name = fst.label, limits = c(0, 0.35), breaks = c(0, 0.1, 0.2, 0.3))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle=90),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12)
  )

ggplot(fst.info, aes(samelakes, Fst))+geom_point()
same<-filter(fst.info, samelakes==TRUE)
dif<-filter(fst.info, samelakes==FALSE)
t.test(same$Fst, dif$Fst)

colnames(lakedistances)<-c("Lake.y", "Lake.x","distance")
fst.distances.d<-full_join(fst.info.d, lakedistances, by=c("Lake.x", "Lake.y"))
fst.distances.d$distance[is.na(fst.distances.d$distance)]<-0
fst.distances.d<-filter(fst.distances.d, !is.na(Group1))

dentfst<-ggplot(fst.distances.d, aes(distance,Fst))+geom_point()+
  xlab("Distance between lakes (km)")+theme_bw()+
  theme(axis.text = element_text(color="black"))
mod<-lm(data=fst.distances.d, Fst~distance)

d.fst.julian<-inner_join(fst.distances.d, info, by=c("Group1"="species_lake_date"))
d.fst.julian<-select(d.fst.julian, Group1, Group2, Fst, Julian)
d.fst.julian<-inner_join(d.fst.julian, info, by=c("Group2"="species_lake_date"))
d.fst.julian$timedif<-abs(d.fst.julian$Julian.x-d.fst.julian$Julian.y)

ggplot(d.fst.julian, aes(timedif, Fst))+geom_point()
mod<-lm(data=d.fst.julian, Fst~timedif)

#Retrocurva Fst
Retro_Fst<-genet.dist(RETRO, method="WC84")%>%round(digits=3)
Retro.df <- melt(as.matrix(Retro_Fst), varnames = c("Group1", "Group2"))

t.test(Dent.df$value, Retro.df$value)

lab_order = c("retrocurva_B_8.31.2015","retrocurva_M_8.24.2015",
              "retrocurva_M_9.21.2015",
              "retrocurva_N_9.9.2015","retrocurva_N_10.21.2015","retrocurva_W_8.7.2015")
fst_mat=as.matrix(Retro_Fst)
fst_mat1=fst_mat[lab_order,]
fst_mat2=fst_mat1[,lab_order]
ind = which(upper.tri(fst_mat2), arr.ind = TRUE)
fst.df = data.frame(Group1 = dimnames(fst_mat2)[[2]][ind[,2]],
                    Group2 = dimnames(fst_mat2)[[1]][ind[,1]],
                    Fst = fst_mat2[ ind ])
fst.df$Group1 = factor(fst.df$Group1, levels = unique(fst.df$Group1))
fst.df$Group2 = factor(fst.df$Group2, levels = unique(fst.df$Group2))
fst.df$Fst[fst.df$Fst < 0] = 0
fst.label = expression(italic("F")[ST])
# Extract middle Fst value for gradient argument
mid = max(fst.df$Fst) / 2
fst2<-ggplot(data = fst.df, aes(x = Group1, y = Group2, fill = Fst))+
  geom_tile(colour = "black")+
  geom_text(aes(label = Fst), color="black", size = 5)+
  scale_fill_gradient2(low = "blue", mid = "pink", high = "red", midpoint = mid, name = fst.label, limits = c(0, 0.35), breaks = c(0, 0.1, 0.2, 0.3))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle=90),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12)
  )
ggplot(fst.info, aes(samelakes, Fst))+geom_point()
same<-filter(fst.info, samelakes==TRUE)
dif<-filter(fst.info, samelakes==FALSE)
t.test(same$Fst, dif$Fst)

fst.info<-inner_join(fst.df, info, by=c("Group1"="species_lake_date"))
fst.info<-select(fst.info, Group1, Group2, Fst, Lake)
fst.info<-inner_join(fst.info, info, by=c("Group2"="species_lake_date"))
fst.info<-select(fst.info, Group1, Group2, Fst, Lake.x, Lake.y)
fst.info$samelakes<-fst.info$Lake.x==fst.info$Lake.y
fst.info.r<-fst.info
fst.info.r$species<-rep(c("retro"), length(fst.info.r$Fst))

t.test(fst.info.r$Fst, fst.info.d$Fst)

plot_grid(fst1,fst2)

fst.distances.r<-full_join(fst.info.r, lakedistances, by=c("Lake.x", "Lake.y"))
fst.distances.r$distance[is.na(fst.distances.r$distance)]<-0
fst.distances.r<-filter(fst.distances.r, !is.na(Group1))

retrofst<-ggplot(fst.distances.r, aes(distance,Fst))+geom_point()+
  xlab("Distance between lakes (km)")+theme_bw()+
  theme(axis.text = element_text(color="black"))
mod<-lm(data=fst.distances.r, Fst~distance)

Fig5<-plot_grid(dentfst, retrofst, labels=c("A","B"))
save_plot("PastPaperFig5.jpg",Fig5, base_width = 8, base_height = 4)

r.fst.julian<-inner_join(fst.distances.r, info, by=c("Group1"="species_lake_date"))
r.fst.julian<-select(r.fst.julian, Group1, Group2, Fst, Julian)
r.fst.julian<-inner_join(r.fst.julian, info, by=c("Group2"="species_lake_date"))
r.fst.julian$timedif<-abs(r.fst.julian$Julian.x-r.fst.julian$Julian.y)

ggplot(r.fst.julian, aes(timedif, Fst))+geom_point()
mod<-lm(data=r.fst.julian, Fst~timedif)

###############################################################################################################
#Make Field Graphs

#Prepare Density Data
fieldDens<-read.csv("2015DensityMichiganCLS25Oct18.csv", header=TRUE, stringsAsFactors = FALSE)
fieldDens$Lake[fieldDens$Lake=="LilAp"]<-"LittleAppleton"
LAKES<-c("North","Mill","LittleAppleton","CrookedW","Bishop","Cedar", "Walsh", "Gosling")
FD.2<-filter(fieldDens,Lake%in%LAKES)
FD.3<-select(FD.2, Lake, Rep, Julian, DentTot, RetroTot, CerioTot, ParvTot)
colnames(FD.3)<-c("Lake","Rep","Julian","dentifera","retrocurva","Cerio","parvula")
FD.4<-gather(FD.3, "Species","Density", 4:7)
FD.5<-FD.4%>%group_by(Lake, Julian, Species)%>%summarise(M.Density=mean(Density))
FD.5$M.Density.cor<-FD.5$M.Density/9 #Correct density
ggplot(FD.5, aes(Julian, M.Density.cor, color=Species))+geom_line()+facet_grid(~Lake)

#Prepare prevalence data
field<- read.table("fielddata2015.txt", header=TRUE, sep="\t", dec=".", strip.white=TRUE) # loads the data
field.2<-filter(field,Lake%in%LAKES)
field.3<-filter(field.2, Total>19)
SPECIES<-c("dentifera","parvula","retrocurva", "Cerio")
field.4<-filter(field.3, Species%in%SPECIES)
field.4$PastInf<-field.4$APasteuria+field.4$JPasteuria+field.4$MalePasteuria+field.4$EphipPastueria+field.4$JPastMetsch+field.4$ACaullPasteuria+field.4$JCaullPasteuria
field.4$PastPrev<-field.4$PastInf/field.4$Total
field.prev<-select(field.4, Species, Julian, Lake, PastPrev)

#combine density and prevalence
Field.dens<-full_join(field.prev, FD.5, by=c("Lake","Julian","Species"))
Field.dens$inf.dens<-Field.dens$M.Density.cor*Field.dens$PastPrev
Field.dens.1<-filter(Field.dens, !is.na(PastPrev))

densitygraph<-ggplot(Field.dens, aes(Julian,M.Density.cor, color=Species, shape=Species))+geom_point()+
  geom_line(size=1)+
  facet_grid(Lake~.)+
  theme_bw()+
  theme(strip.text = element_blank())+
  theme(strip.background = element_blank())+
  labs(x="Ordinal day", y=bquote('Host density'~(animals/m^2~'x'~10^5)))+
  scale_color_manual(values=c(Color.Palette), name="Host species", labels=c("Ceriodaphnia","D. dentifera","D. parvula","D. retrocurva"))+
  scale_shape_manual(values=c(16,17,0,1),name="", labels=c("Ceriodaphnia","D. dentifera","D. parvula","D. retrocurva"))+
  theme(axis.title = element_text(size=14))+
  theme(axis.text = element_text(size=12, color="black"))+
  scale_x_continuous(limits=c(205, 310), breaks = c(200,225,250,275,300), labels=c("Jul. 19","Aug. 13", "Sept. 7","Oct. 2","Oct. 27"))+
  scale_y_continuous(breaks=c(0, 100000, 200000,300000), labels=c(0,1,2,3))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle=90))+
  theme(legend.position = "none")

infdensitygraph<-ggplot(Field.dens.1, aes(Julian,inf.dens, color=Species, shape=Species))+geom_point()+
  geom_line(size=1, lty="dotdash")+
  facet_grid(Lake~.)+
  theme_bw()+
  labs(x="Ordinal day", y=bquote('Infected host density'~(animals/m^2~'x'~10^3)))+
  scale_color_manual(values=c(Color.Palette), name="Host species", labels=c("Ceriodaphnia","D. dentifera","D. parvula","D. retrocurva"))+
  scale_shape_manual(values=c(16,17,0,1),name="", labels=c("Ceriodaphnia","D. dentifera","D. parvula","D. retrocurva"))+
  theme(axis.title = element_text(size=14))+
  theme(axis.text = element_text(size=12, color="black"))+
  scale_x_continuous(limits=c(205, 310), breaks = c(200,225,250,275,300), labels=c("Jul. 19","Aug. 13", "Sept. 7","Oct. 2","Oct. 27"))+
  scale_y_continuous(breaks=c(0,5000,10000), labels=c(0,5,10))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle=90))+
  theme(legend.title = element_text(size=14))+
  theme(strip.text = element_text(size=14))+
  theme(strip.background = element_blank())+
  theme(legend.text = element_text(face="italic", size=12))

prevalencegraph<-ggplot(Field.dens.1, aes(Julian,PastPrev, color=Species, shape=Species))+geom_point()+
  geom_line(size=1, lty="longdash")+
  facet_grid(Lake~.)+
  theme_bw()+
  theme(strip.background = element_blank())+
  ylab("Infection prevalence")+
  scale_color_manual(values=c(Color.Palette), name="Host species", labels=c("Ceriodaphnia","D. dentifera","D. parvula","D. retrocurva"))+
  scale_shape_manual(values=c(16,17,0,1),name="", labels=c("Ceriodaphnia","D. dentifera","D. parvula","D. retrocurva"))+
  theme(axis.title = element_text(size=14))+
  theme(axis.text = element_text(size=12, color="black"))+
  scale_x_continuous(limits=c(205, 310), breaks = c(200,225,250,275,300), labels=c("Jul. 19","Aug. 13", "Sept. 7","Oct. 2","Oct. 27"))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")+
  theme(strip.text = element_blank())+
  theme(strip.background = element_blank())+
  theme(axis.text.x = element_text(angle=90))

Fig1<-plot_grid(densitygraph, prevalencegraph, infdensitygraph, ncol = 3, rel_widths = c(0.8,0.8,1.5), labels = c("A","B","C"))

save_plot("Pasteuria2015Figure1.jpg", Fig1, base_height = 11, base_width = 9)

for.wide<-select(Field.dens.1, Julian, Lake, Species, inf.dens)
wide.data<-spread(for.wide, Species, inf.dens)
wide.data$dentifera[is.na(wide.data$dentifera)]<-0 
wide.data$retrocurva[is.na(wide.data$retrocurva)]<-0 

Figure2<-ggplot(wide.data, aes(dentifera, retrocurva, color=Lake))+geom_point(size=2)+theme_bw()+
  labs(x=bquote('Infected D. dentifera density'~(animals/m^2)), y=bquote('Infected D. retrocurva density'~(animals/m^2)))+
  theme(axis.text = element_text(color="black"))+scale_y_continuous(limits=c(0,5000))+scale_x_continuous(limits=c(0,5000))
save_plot("PastPaperFig2.jpg",Figure2, base_width = 6, base_height = 5)

######################################################################################################################
##Repeat genetic analyses with coinfection data
coinf.1<-select(data.coinf, sample, Species, Lake, collection.date,Julian, p1, p2, p3, p4, p7, p11, p12, p16, p18, p19)
coinf.2<-coinf.1[rowSums(is.na(coinf.1)) <= 2, ] # include samples that amplified at 8 or more out of 10 loci
coinf.2$species_lake_date<- paste(coinf.2$Species,coinf.2$Lake,coinf.2$collection.date,sep="_")

coinf.2<-filter(coinf.2, sample!="B.R1.10.28.2015(2/6).co1")
coinf.2<-filter(coinf.2, sample!="B.R1.10.28.2015(2/6).co2")
coinf.2<-filter(coinf.2, sample!="B.R1.8.4.15(B.R1.10.28.15.1).co1")
coinf.2<-filter(coinf.2, sample!="B.R1.8.4.15(B.R1.10.28.15.1).co2")
coinf.2<-filter(coinf.2, sample!="B.R2.8.4.15(B.R2.10.28.15)")
coinf.2<-filter(coinf.2, sample!="B.R2.10.28.2015")
coinf.2$sample[coinf.2$sample=="B.D5.8.4.2015(B.D5.10.28.15)"]<-"B.D5.8.4.15"
coinf.2$sample[coinf.2$sample=="B.R5.8.4.15(B.R5.10.28.2015)"]<-"B.R5.8.4.15"

info<-unique(coinf.2[,c("species_lake_date","Species","Lake","collection.date","Julian")])

#get into format for poppr
dat<-select(coinf.2, sample, species_lake_date, p1, p2, p3, p4, p7, p11, p12, p16, p18, p19)

dentifera<-filter(coinf.2, Species=="dentifera")
dent.data<-select(dentifera, sample, species_lake_date,p1, p2, p3, p4, p7,p11,p12, p16,p18,p19)
dentgroups<-dent.data%>%group_by(species_lake_date)%>%summarise(count=n())
dent.data<-filter(dent.data, species_lake_date!="dentifera_B_8.17.2015") #only one of these, so don't include in fst analysis

retrocurva<-filter(coinf.2, Species=="retrocurva")
retro.data<-select(retrocurva, sample, species_lake_date,p1, p2, p3, p4, p7,p11,p12, p16,p18,p19)
retrogroups<-retro.data%>%group_by(species_lake_date)%>%summarise(count=n())
retro.data<-filter(retro.data, species_lake_date!="retrocurva_B_8.4.2015") #only one of these, so don't include in fst analysis
retro.data<-filter(retro.data, species_lake_date!="retrocurva_CW_10.9.2015") #only one of these, so don't include in fst analysis
retro.data<-filter(retro.data, species_lake_date!="retrocurva_M_9.8.2015") #only one of these, so don't include in fst analysis

write.table(dat, "Past2015coinf.txt", col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE)
write.table(dent.data, "dent.data.coinf.txt", col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE)
write.table(retro.data, "retro.data.coinf.txt", col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE)

z <- read.loci("Past2015coinf.txt", col.loci=3:12,row.names=1, col.pop=2, header=TRUE, allele.sep="-")
dent<-read.loci("dent.data.coinf.txt", col.loci=3:12,row.names=1, col.pop=2, header=TRUE, allele.sep="-")
retro<-read.loci("retro.data.coinf.txt", col.loci=3:12,row.names=1, col.pop=2, header=TRUE, allele.sep="-")

Z<-loci2genind(z, ploidy = 1)
DENT<-loci2genind(dent, ploidy = 1)
RETRO<-loci2genind(retro, ploidy = 1)

Y<-genind2genpop(Z)

genind2genalex(Z, filename = "Past2015coinf.genalex.csv", quiet = FALSE, pop = NULL,
               allstrata = TRUE, geo = FALSE, geodf = "xy", sep = ",",
               sequence = FALSE, overwrite=TRUE)

COINF<-read.genalex("Past2015coinf.genalex.csv")
splitStrata(COINF)<-~species/lake/date
setPop(COINF)<-~species/lake/date

tab(COINF) #look at the multilocus genotype table
nmll(COINF) #count the number of multilocus genotypes
mll(COINF)#show the multilocus genotype definitions

Coinf.Distance<-provesti.dist(COINF) #This seems to be using contracted MLGs.

HOST <- factor(coinf.2$Species)
LAKE <- factor(coinf.2$Lake)

full2015.Amova.co<-poppr.amova(COINF, ~species/lake/date, within=FALSE, filter=TRUE,threshold =0 , missing="ignore",
                            correction="lingoes", cutoff=0)

sig.full2015.co<-randtest(full2015.Amova.co, nrepet=1000)
#same as without coinfection

#Flip species and lake
setPop(COINF)<-~lake/species/date
full2015.Amova.Flip.co<-poppr.amova(COINF, ~lake/species/date, within=FALSE, filter=TRUE,threshold =0 , missing="ignore",
                                 correction="lingoes", cutoff=0)
sig.full2015.Flip.co<-randtest(full2015.Amova.Flip.co, nrepet=1000)
#same as without coinfection.

#Mantel tests
infolake<-select(coinf.2,Lake)
infolake$Lake<-as.factor(infolake$Lake)
distanceLake<-daisy(infolake, metric ="gower") #0 if same; 1 if different
infoSpecies<-select(coinf.2, Species)
infoSpecies$Species<-as.factor(infoSpecies$Species)
distanceSpecies<-daisy(infoSpecies, metric="gower")
infodate<-select(coinf.2, Julian)
distancedate<-daisy(infodate, metric ="manhattan") 

mantel(Coinf.Distance, distanceLake, method="pearson", permutations=999)
mantel(Coinf.Distance, distanceSpecies, method="pearson", permutations=999)
mantel(Coinf.Distance, distancedate, method="pearson", permutations=999)

mantel.partial(Coinf.Distance, distanceLake, distanceSpecies, method="pearson", permutations=999)
mantel.partial(Coinf.Distance, distanceSpecies, distanceLake, method="pearson", permutations=999)
mantel.partial(Coinf.Distance, distanceSpecies, distancedate, method="pearson", permutations=999)
mantel.partial(Coinf.Distance, distanceLake, distancedate, method="pearson", permutations=999)

###############################################################################################################
#Compute pairwise Fst

#dentifera Fst
Dent_Fst<-genet.dist(DENT, method="WC84")%>%round(digits=3)
Dent.df <- melt(as.matrix(Dent_Fst), varnames = c("Group1", "Group2"))

fst_mat=as.matrix(Dent_Fst)
ind = which(upper.tri(fst_mat), arr.ind = TRUE)
fst.df = data.frame(Group1 = dimnames(fst_mat)[[2]][ind[,2]],
                    Group2 = dimnames(fst_mat)[[1]][ind[,1]],
                    Fst = fst_mat[ ind ])
fst.df$Group1 = factor(fst.df$Group1, levels = unique(fst.df$Group1))
fst.df$Group2 = factor(fst.df$Group2, levels = unique(fst.df$Group2))
fst.df$Fst[fst.df$Fst < 0] = 0
fst.label = expression(italic("F")[ST])

fst.info<-inner_join(fst.df, info, by=c("Group1"="species_lake_date"))
fst.info<-select(fst.info, Group1, Group2, Fst, Lake)
fst.info<-inner_join(fst.info, info, by=c("Group2"="species_lake_date"))
fst.info<-select(fst.info, Group1, Group2, Fst, Lake.x, Lake.y)
fst.info$samelakes<-fst.info$Lake.x==fst.info$Lake.y
fst.info.d<-fst.info
fst.info.d$species<-rep(c("dent"), length(fst.info.d$Fst))
fst.info.d<-filter(fst.info.d, !is.na(Lake.y))

same<-filter(fst.info, samelakes==TRUE)
dif<-filter(fst.info, samelakes==FALSE)
t.test(same$Fst, dif$Fst)

colnames(lakedistances)<-c("Lake.y", "Lake.x","distance")
fst.distances.d<-full_join(fst.info.d, lakedistances, by=c("Lake.x", "Lake.y"))
fst.distances.d$distance[is.na(fst.distances.d$distance)]<-0
fst.distances.d<-filter(fst.distances.d, !is.na(Group1))

ggplot(fst.distances.d, aes(distance,Fst))+geom_point()+
  xlab("Distance between lakes (km)")+theme_bw()+
  theme(axis.text = element_text(color="black"))
mod<-lm(data=fst.distances.d, Fst~distance)

d.fst.julian<-inner_join(fst.distances.d, info, by=c("Group1"="species_lake_date"))
d.fst.julian<-select(d.fst.julian, Group1, Group2, Fst, Julian)
d.fst.julian<-inner_join(d.fst.julian, info, by=c("Group2"="species_lake_date"))
d.fst.julian$timedif<-abs(d.fst.julian$Julian.x-d.fst.julian$Julian.y)

ggplot(d.fst.julian, aes(timedif, Fst))+geom_point()
mod<-lm(data=d.fst.julian, Fst~timedif)

#Retrocurva Fst
Retro_Fst<-genet.dist(RETRO, method="WC84")%>%round(digits=3)
Retro.df <- melt(as.matrix(Retro_Fst), varnames = c("Group1", "Group2"))

t.test(Dent.df$value, Retro.df$value)

fst_mat=as.matrix(Retro_Fst)
ind = which(upper.tri(fst_mat), arr.ind = TRUE)
fst.df = data.frame(Group1 = dimnames(fst_mat)[[2]][ind[,2]],
                    Group2 = dimnames(fst_mat)[[1]][ind[,1]],
                    Fst = fst_mat[ ind ])
fst.df$Group1 = factor(fst.df$Group1, levels = unique(fst.df$Group1))
fst.df$Group2 = factor(fst.df$Group2, levels = unique(fst.df$Group2))
fst.df$Fst[fst.df$Fst < 0] = 0

fst.info<-inner_join(fst.df, info, by=c("Group1"="species_lake_date"))
fst.info<-select(fst.info, Group1, Group2, Fst, Lake)
fst.info<-inner_join(fst.info, info, by=c("Group2"="species_lake_date"))
fst.info<-select(fst.info, Group1, Group2, Fst, Lake.x, Lake.y)
fst.info$samelakes<-fst.info$Lake.x==fst.info$Lake.y
fst.info.r<-fst.info
fst.info.r$species<-rep(c("retro"), length(fst.info.r$Fst))

same<-filter(fst.info.r, samelakes==TRUE)
dif<-filter(fst.info.r, samelakes==FALSE)
t.test(same$Fst, dif$Fst)

t.test(fst.info.r$Fst, fst.info.d$Fst)

fst.distances.r<-full_join(fst.info.r, lakedistances, by=c("Lake.x", "Lake.y"))
fst.distances.r$distance[is.na(fst.distances.r$distance)]<-0
fst.distances.r<-filter(fst.distances.r, !is.na(Group1))
  

ggplot(fst.distances.r, aes(distance,Fst))+geom_point()+
  xlab("Distance between lakes (km)")+theme_bw()+
  theme(axis.text = element_text(color="black"))
mod<-lm(data=fst.distances.r, Fst~distance)

r.fst.julian<-inner_join(fst.distances.r, info, by=c("Group1"="species_lake_date"))
r.fst.julian<-select(r.fst.julian, Group1, Group2, Fst, Julian)
r.fst.julian<-inner_join(r.fst.julian, info, by=c("Group2"="species_lake_date"))
r.fst.julian$timedif<-abs(r.fst.julian$Julian.x-r.fst.julian$Julian.y)

ggplot(r.fst.julian, aes(timedif, Fst))+geom_point()
mod<-lm(data=r.fst.julian, Fst~timedif)
