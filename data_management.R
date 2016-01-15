setwd("/media/data/final_bac/full")

##### How to sort based on multiple other variables########
om_otu_bor_omhorzonesort<-om_otu_bor[order(om_otu_bor[,"bac_om"],om_otu_bor[,"bac_hor"], om_otu_bor[ ,"bac_zone"]),]

#sort based on row.names
x.sort<-x[order(row.names(x)),]

#delete selected column dataframe
x$column.name<-NULL
#delete column in matrix
x<-x[,-column number]

####subset based on column in other matrix/dataframe
txb03f100_otu<-subset(b03_f100_1_nofire, bac_env_nofire$Zone=="TX")
org_env<-subset(fung_env, fung_env$horizon=="1")
min_env<-subset(fung_env, fung_env$horizon=="2")

###subset a distance matrix using logical variables
org<-fung_env$horizon=="1"
min<-fung_env$horizon=="2"
org_tyc03<-as.dist(f.thyc[org,org])
min_tyc03<-as.dist(f.thyc[min, min])

#add row and col names back to matrix after using full to convert lt dissimilarity matrix to full
f.tyc.full<-full(as.dist(f.tyc))
rownames(MD.2.tyc)<-c("BL026","BL028", "BL030", "BL032", "BL034", "BL038", "BL040", "BL044", "BL046", "BL048", "BR050", "BR052", "BR054", "BR056", "BR058", "BR060", "BR062", "BR064", "BR066", "BR068", "BR070", "BR072", "LH002", "LH004", "LH006", "LH008", "LH010", "LH012", "LH014", "LH016", "LH018", "LH020", "LH022", "LH024")



####remove rare species with specnumber in vegan
freq<-specnumber(bcSBSb03f100_otu, MARGIN=2)
bSBS_otu_rare2<-bcSBSb03f100_otu[,freq>6]

###another way to remove rares
maxab <- apply(otu, 2, max)
n1 <- names(which(maxab < 1000))
otu.ab <- otu[,-which(names(otu) %in% n1)]
#bonus of this approach, you can use that vector on another object like tax info of otus
taxa.ab <- taxa[-which(taxa$OTU %in% n1),]

##select otu's out of taxonomy
gawk 'FNR==NR {keys[$1]; next} $1 in keys' ../R/b03f100_om_inv_sel bac_final.an.0.03cons.taxonomy >b03f100_om_inv_sel.tax
#In R.  pull taxonomy info for selected otus
f.taxa <- read.table(file="../fung_real/fung.all.mar14.hil.0.945.cons.taxonomy", header=T)
f.glm.taxa <- f.taxa[f.taxa$OTU %in% f.sig.otu,]


#######make otu matrix with just selected otu
gawk '{print $1}' b03f100_hor_sel>b03f100_hor_bor.accnos
#transpose otu matrix in R
sed 's/"//g' b03f100_1_nofire_otu > b03f100_1_nofire_otu1
sed 's/,/ /g' b03f100_1_nofire_otu1 > b03f100_1_nofire_otu2
gawk 'FNR==NR {keys[$1]; next} $1 in keys' b03f100_zone_bor.accnos b03f100_1_nofire_otu2 > b03100_1_nofire_otu_zone_bor %paste sample names as top row in text editor, 
#or do it in R which is way easier
f.sig.otu <- names(f.sig)
f.glm.otu <- f.otu[, colnames(f.otu) %in% f.sig.otu]

#pull selected otu's from fasta
#create file with just otu
sort glm_selected.txt | uniq | sed "s/Otu0\{0,\}//" > glm_selected_filtered.txt 
cat glm_selected_filtered.txt| while read line; do grep -A 1 "\s${line}[|]" bac_final.0.03.rep.fasta; done >glm_sel.fasta


#find samples missing from one dataframe but present in second
setdiff(row.names(f.env), row.names(f.otu))


#######create color vector from data
taxa <- read.table(file="abundant.taxa.csv", sep=",", header=T)

###works with ggplot fill variable
bac.phyla.color<-c("Alphaproteobacteria" = "#A6CEE3", "Betaproteobacteria" = "#7DB4D5", 
                   "Deltaproteobacteria" = "#5C9FC9", "Gammaproteobacteria" = "#3A89BD", 
                   "Proteobacteria" = "#1F78B4", "Actinobacteria" = "#B2DF8A", "Acidobacteria" = "#79C360", 
                   "Bacteroidetes" = "#33A02C", "Firmicutes" = "#FF7F00", "Fusobacteria" = "#E31A1C", 
                   "Cyanobacteria" = "#FDBF6F", "Gemmatimonadetes" = "#FB9A99", "Planctomycetes" = "#CAB2D6",
                   "Verrucomicrobia" = "#6A3D9A", "Nitrospirae" = "#FFFF99", "other" = "#d3d3d3", "unclassified" = "#000000")

#### works with
phyla.f <- as.character(taxa$phyla)
phyla.color <- revalue(phyla.f, c("Alphaproteobacteria" = "#A6CEE3", "Betaproteobacteria" = "#7DB4D5",
                                  "Deltaproteobacteria" = "#5C9FC9", "Gammaproteobacteria" = "#3A89BD", 
                                  "Proteobacteria" = "#1F78B4", "Actinobacteria" = "#B2DF8A", 
                                  "Acidobacteria" = "#79C360", "Bacteroidetes" = "#33A02C", 
                                  "Firmicutes" = "#FF7F00", "Fusobacteria" = "#E31A1C", 
                                  "Cyanobacteria" = "#FDBF6F", "Gemmatimonadetes" = "#FB9A99", 
                                  "Planctomycetes" = "#CAB2D6", "Verrucomicrobia" = "#6A3D9A", 
                                  "Nitrospirae" = "#FFFF99", "other" = "#d3d3d3", "unclassified" = "#000000"))

### read in shared file
otu <- read.table(file="fancher.all.trim.contigs.good.unique.good.filter.precluster.pick.pick.an.unique_list.0.03.subsample.shared", header=T, row.names=2)
otu$label <- NULL
otu$nunumOtus <- NULL
