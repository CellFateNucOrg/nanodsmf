# preparing dataset for package

MSssI_CpG<-read.delim("/Users/semple/Documents/MeisterLab/sequencingData/20171027_Minion_TMP_Meth/meth_calls/20171027_pass_barcode06_CpGcalls.tsv",
                      stringsAsFactors=F)
save(MSssI_CpG,file="data/MSssI_CpG.RData")


MSssI_GpC<-read.delim("/Users/semple/Documents/MeisterLab/sequencingData/20171027_Minion_TMP_Meth/meth_calls/20171027_pass_barcode06_GpCcalls.tsv",
                      stringsAsFactors=F)
save(MSssI_GpC,file="data/MSssI_GpC.RData")


MCviPI_CpG<-read.delim("/Users/semple/Documents/MeisterLab/sequencingData/20171027_Minion_TMP_Meth/meth_calls/20171027_pass_barcode07_CpGcalls.tsv",
                      stringsAsFactors=F)
save(MCviPI_CpG,file="data/MCviPI_CpG.RData")


MCviPI_GpC<-read.delim("/Users/semple/Documents/MeisterLab/sequencingData/20171027_Minion_TMP_Meth/meth_calls/20171027_pass_barcode07_GpCcalls.tsv",
                      stringsAsFactors=F)
save(MCviPI_GpC,file="data/MCviPI_GpC.RData")


Control_CpG<-read.delim("/Users/semple/Documents/MeisterLab/sequencingData/20171027_Minion_TMP_Meth/meth_calls/20171027_pass_barcode08_CpGcalls.tsv",
                       stringsAsFactors=F)
save(Control_CpG,file="data/Control_CpG.RData")


Control_GpC<-read.delim("/Users/semple/Documents/MeisterLab/sequencingData/20171027_Minion_TMP_Meth/meth_calls/20171027_pass_barcode08_GpCcalls.tsv",
                       stringsAsFactors=F)
save(Control_GpC,file="data/Control_GpC.RData")


genomeFile<-Biostrings::readDNAStringSet("/Users/semple/Documents/MeisterLab/GenomeVer/PCR_wPM28_32/PCR_wPM28_32.fasta")
save(genomeFile,file="data/genomeFile.RData")

genomeGR<-nanodsmf::DNAstringset2gr(genomeFile)
save(genomeGR,file="data/genomeGR.RData")

