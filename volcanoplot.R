setwd('D:/')
res = read.table('./data/89_F0405_res_0.3.txt',sep = '\t')
ccrg = read.table('./data2/Data_Sheet_1_Molecular Classification Based on Prognostic and Cell Cycle-Associated Genes in Patients With Colon Cancer.CSV', header=1)
#1857


ccrg_res =res[which(rownames(res) %in% ccrg$cellcyclemarker),]
dim(ccrg_res)#1584 6

range(ccrg_res$log2FoldChange)

cut_lfc <- 3
cut_pvalue <- 0.05

par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(ccrg_res)

# P values
with(topT, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", col='grey', cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value)))

with(subset(topT, pvalue<cut_pvalue & log2FoldChange>cut_lfc), points(log2FoldChange, -log10(pvalue), pch=20, col='brown1', cex=1.5))
with(subset(topT, pvalue<cut_pvalue & log2FoldChange<(-cut_lfc)), points(log2FoldChange, -log10(pvalue), pch=20, col='darkblue', cex=1.5))

## Add lines for FC and P-value cut-off
abline(v=0, col='black', lty=3, lwd=1.0)
abline(v=-cut_lfc, col='black', lty=4, lwd=2.0)
abline(v=cut_lfc, col='black', lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$pvalue<cut_pvalue], na.rm=TRUE)), col='black', lty=4, lwd=2.0)


#up gene, dn gene 수
dn_gene =rownames(ccrg_res[which(ccrg_res$pvalue<0.05 & ccrg_res$log2FoldChange< -3),])
up_gene =rownames(ccrg_res[which(ccrg_res$pvalue<0.05 & ccrg_res$log2FoldChange>3),])

write.table(dn_gene, './data/dn_gene_ccrg_res.txt')
write.table(up_gene,'./data/up_gene_ccrg_res.txt')


