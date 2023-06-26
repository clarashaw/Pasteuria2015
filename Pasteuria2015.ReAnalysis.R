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

#Take a look at 2015 Data only

data<-read.csv("P.Microsats7March19.csv", header=TRUE, stringsAsFactors = FALSE)
data.coinf<-read.csv("P.Microsats7March19.coinf.csv", header=TRUE, stringsAsFactors = FALSE)

data2015<-filter(data, year==2015)

#dropping p17 - the same all the way down
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

write.table(dat, "Past2015.txt", col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE)
z <- read.loci("Past2015.txt", col.loci=3:12,row.names=1, col.pop=2, header=TRUE, allele.sep="-")
Z<-loci2genind(z, ploidy = 1)
Y<-genind2genpop(Z)

genind2genalex(Z, filename = "Past2015genalex.csv", quiet = FALSE, pop = NULL,
               allstrata = TRUE, geo = FALSE, geodf = "xy", sep = ",",
               sequence = FALSE, overwrite=TRUE)

Past<-read.genalex("Past2015genalex.csv")

splitStrata(Past)<-~species/lake/date
setPop(Past)<-~species/lake/date

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
full2015.Amova<-poppr.amova(Past, ~species/lake/date, within=FALSE, filter=TRUE,threshold =0 , missing="ignore",
                            correction="lingoes", cutoff=0)

sig.full2015<-randtest(full2015.Amova, nrepet=1000)

#Flip species and lake
setPop(Past)<-~lake/species/date
full2015.Amova.Flip<-poppr.amova(Past, ~lake/species/date, within=FALSE, filter=TRUE,threshold =0 , missing="ignore",
                                 correction="lingoes", cutoff=0.0)
sig.full2015.Flip<-randtest(full2015.Amova.Flip, nrepet=1000)

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
##Repeat first analysis with coinfection data
coinf.1<-select(data.coinf, sample, Species, Lake, collection.date, p1, p2, p3, p4, p7, p11, p12, p16, p18, p19)
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

info<-unique(coinf.2[,c("species_lake_date","Species","Lake","collection.date")])

#get into format for poppr
dat<-select(coinf.2, sample, species_lake_date, p1, p2, p3, p4, p7,p11,p12, p16,p18,p19 )

write.table(dat, "Past2015coinf.txt", col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE)
z <- read.loci("Past2015coinf.txt", col.loci=3:12,row.names=1, col.pop=2, header=TRUE, allele.sep="-")
Z<-loci2genind(z, ploidy = 1)
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

#Save this with witdth of 1000 and height of 900 for things to line up. Legends still annoying
plot.phylo(upgma(Coinf.Distance),tip.color=cbPalette[HOST],cex=1.5,no.margin=TRUE)
legend(x=3, y=3, legend=c("Cerio", "D.dentifera","D. parvula", "D. retrocurva"),title="Host",
       col=c("#999999", "#E69F00", "#56B4E9", "#009E73"), lty=1, cex=2, text.font =3, bty="n")


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

