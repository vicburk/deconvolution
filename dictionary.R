
source("packages.R")

################################
# Map probes to gene symbols
################################

gene.expression <- read_xlsx("processed_data/gene_expression_clean.xlsx")

x = hgu219SYMBOL
mapped_probes = mappedkeys(x)

xx = as.list(x[mapped_probes])

gene.symbols = unlist(xx[names(xx) %in% gene.expression$probeset_id])
gene.symbols

##################################
# Map probes to entrez ID
##################################

y = hgu219ENTREZID
mapped_probesy = mappedkeys(y)

yy = as.list(y[mapped_probesy])

length(unlist(yy[names(yy) %in% gene.expression$probeset_id]))

gene.entrez = unlist(yy[names(yy) %in% gene.expression$probeset_id])

length(gene.entrez)
length(gene.symbols)

# Equal number of symbols and entrez IDs were mapped

##################################
# Which probes were unidentified?
##################################

uig = gene.expression$probeset_id[!gene.expression$probeset_id %in% names(xx)]

uig

nrow(gene.expression)
length(xx)

##################################
# Map unidentified probes to 
# PrimeView chip data Symbols
##################################

dictionary = read.delim("data/Human_AFFY_PrimeView_MSigDB.v7.2.chip",header = TRUE)
dictionary = dictionary[,1:2]
head(dictionary)

gene.symbols2 = dictionary[dictionary$Probe.Set.ID %in% uig,2]
names(gene.symbols2) = dictionary[dictionary$Probe.Set.ID %in% uig,1]

######################################
# Which probes are still unidentified?
######################################

uig2 = uig[!uig %in% dictionary$Probe.Set.ID]

gene.symbols2

#########################################
# Map unidentified symbols to entrez IDs
#########################################

hs = Homo.sapiens
keytypes(hs)
gene.entrez2 = HDO.db::select(hs,
                      columns= c("SYMBOL","ENTREZID"),
                      keys = gene.symbols2,
                      keytype = "SYMBOL")

dim(gene.entrez2)
length(gene.symbols2)

############################################
# Which symbols have multiple entrez IDs?
############################################

which(!table(gene.entrez2$SYMBOL) == table(gene.symbols2))

gene.entrez2[which(gene.entrez2$SYMBOL == "MEMO1"),]

gene.symbols2[which(gene.symbols2 == "MEMO1")]

# 51072 is the correct symbol, so remove 7795.

gene.entrez2 = gene.entrez2[gene.entrez2$ENTREZID != "7795",]

dim(gene.entrez2)
length(gene.symbols2)

gene.symbols2
gene.entrez2

any(duplicated(unique(gene.entrez2)$SYMBOL)) # No Symbols have multiple entrez IDs

############################################
# Which symbols were unidentified?
############################################

na.values = gene.symbols2[which(is.na(gene.entrez2$SYMBOL))]

#############################################
# Convert unidentified symbols to entrez ID
# Using alias
#############################################

alias = HDO.db::select(hs,
               columns = c("ALIAS","ENTREZID","SYMBOL","GENENAME","GENETYPE"),
               keys = na.values,
               keytype = "ALIAS")


head(alias,10)
head(na.values,10)

dim(alias)

#########################################
# Which unidentified genes still lack
# an entrez ID
#########################################

omitted = alias$ALIAS[is.na(alias$SYMBOL)]
omitted

na.values

# Remove unidentified entrez IDs

alias2 = na.omit(alias)
length(na.values)

dim(unique(alias2))

####################################
# Unique alias and ID combinations
####################################

alias3 = unique(alias2)
alias3

#########################################
# Which Aliases have multiple entrez IDs?
#########################################

duplicated.genes = unique(alias3[duplicated(alias3$ALIAS),1])
duplicated.genes

alias3[alias3$ALIAS %in% duplicated.genes,]

# Correct entrez IDs:
## GBA: 2629
## ELOA3: 162699
## H4-16: 121504

##########################################
# Remove duplicated aliases
##########################################

alias4 = alias3[alias3$ENTREZID %in% c("2629","162699", "121504"),]

alias5 = alias3[!alias3$ALIAS %in% duplicated.genes,]

alias6 = rbind(alias5,alias4)

alias6


#hgnc1 = gsub("\\s\\[.*","\\1",dictionary[dictionary$Gene.Symbol %in% u.alias, 3])
#hgnc2 = gsub(".*HGNC:(.+)].*","\\1",dictionary[dictionary$Gene.Symbol %in% u.alias, 3])

#cbind(dictionary[dictionary$Gene.Symbol %in% u.alias, 2],hgnc1)

#####################################################
# Replace NA values with correct Alias
#####################################################

remaining = na.values[!na.values %in% omitted]

gene.entrez2[which(names(gene.symbols2) %in% names(remaining)),1] = remaining

alias6

gene.entrez2[is.na(gene.entrez2$ENTREZID),]

remaining

####################################################
# Remove NA values
####################################################

gene.entrez.symbol2 = cbind(gene.entrez2, names(gene.symbols2))

gene.entrez.symbol3 = gene.entrez.symbol2[!is.na(gene.entrez.symbol2$SYMBOL),]

##########################################
# Match aliases with entrez IDs
##########################################

alias6$ALIAS
alias6$ENTREZID


gene.entrez.symbol3[gene.entrez.symbol3$SYMBOL == alias6$ALIAS[1],]


gene.entrez.symbol4 = gene.entrez.symbol3

# gene.entrez.symbol4[which(gene.entrez.symbol4$SYMBOL == alias6$ALIAS[1]),]
# alias6$ENTREZID[1]
for (i in 1:length(alias6$ALIAS)) {
  gene.entrez.symbol4[which(gene.entrez.symbol4$SYMBOL == alias6$ALIAS[i]),2] = alias6$ENTREZID[i]
}

gene.entrez.symbol4

all(names(gene.entrez)==names(gene.symbols))

gene.entrez.symbol = cbind(gene.symbols,gene.entrez,names(gene.symbols))
rownames(gene.entrez.symbol) = NULL
colnames(gene.entrez.symbol) = c("SYMBOL","ENTREZID","PROBE")
colnames(gene.entrez.symbol4) = c("SYMBOL","ENTREZID","PROBE")

gene.entrez.symbol.all = rbind(gene.entrez.symbol,gene.entrez.symbol4)

dim(gene.entrez.symbol.all)
dim(gene.expression)

# Only 617 probes were unidentified
nrow(gene.expression) - nrow(gene.entrez.symbol.all)

gene.entrez.symbol.all = gene.entrez.symbol.all[,c(3,1,2)]

write.xlsx(gene.entrez.symbol.all, "processed_data/dictionary.xlsx")
