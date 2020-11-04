library(ape)
library(spider)

setwd("Desktop/test/")

correcta <- read.dna('20180713_COI_correcta_edit.fas', format='fasta')
HAcorrecta <- haploAccum(correcta, method='random', permutations=1000)
plot(HAcorrecta, ylim=c(0,200), xlim=c(0,450), col='blue', ci.type=c('polygon'), ci.col='lightblue', ci.lty=0)

zonata <- read.dna('renamed_20180713_COI_zonata_edit.fas', format='fasta')
HAzonata <- haploAccum(zonata, method='random', permutations=1000)
plot(HAzonata, col='darkgreen', ci.type=c('polygon'), ci.col='lightgreen', ci.lty=0, add=TRUE)

cucurbitae <- read.dna('renamed_20180713_COI_cucurbitae_Edit.fas', format='fasta')
HAcucurbitae <- haploAccum(cucurbitae, method='random', permutations=1000)
plot(HAcucurbitae, col='red', ci.type=c('polygon'), ci.col='orange', ci.lty=0, add=TRUE)

chaoHaplo(zonata)
chaoHaplo(cucurbitae)
chaoHaplo(correcta)