# preparing dataset for package

# single molecule data
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




# frequency data
MSssI_freqCmG<-read.delim("/Users/semple/Documents/MeisterLab/sequencingData/20171027_Minion_TMP_Meth/meth_freq/20171027_pass_barcode06_freqCmG.tsv",
                      stringsAsFactors=F)
save(MSssI_freqCmG,file="data/MSssI_freqCmG.RData")

MSssI_freqGCm<-read.delim("/Users/semple/Documents/MeisterLab/sequencingData/20171027_Minion_TMP_Meth/meth_freq/20171027_pass_barcode06_freqGCm.tsv",
                      stringsAsFactors=F)
save(MSssI_freqGCm,file="data/MSssI_freqGCm.RData")

MCviPI_freqCmG<-read.delim("/Users/semple/Documents/MeisterLab/sequencingData/20171027_Minion_TMP_Meth/meth_freq/20171027_pass_barcode07_freqCmG.tsv",
                       stringsAsFactors=F)
save(MCviPI_freqCmG,file="data/MCviPI_freqCmG.RData")

MCviPI_freqGCm<-read.delim("/Users/semple/Documents/MeisterLab/sequencingData/20171027_Minion_TMP_Meth/meth_freq/20171027_pass_barcode07_freqGCm.tsv",
                       stringsAsFactors=F)
save(MCviPI_freqGCm,file="data/MCviPI_freqGCm.RData")

Control_freqCmG<-read.delim("/Users/semple/Documents/MeisterLab/sequencingData/20171027_Minion_TMP_Meth/meth_freq/20171027_pass_barcode08_freqCmG.tsv",
                        stringsAsFactors=F)
save(Control_freqCmG,file="data/Control_freqCmG.RData")

Control_freqGCm<-read.delim("/Users/semple/Documents/MeisterLab/sequencingData/20171027_Minion_TMP_Meth/meth_freq/20171027_pass_barcode08_freqGCm.tsv",
                        stringsAsFactors=F)
save(Control_freqGCm,file="data/Control_freqGCm.RData")




# PCR fragment ("genome") data
ttTi5605dna<-Biostrings::readDNAStringSet("/Users/semple/Documents/MeisterLab/GenomeVer/PCR_wPM28_32/PCR_wPM28_32.fasta")
save(ttTi5605dna,file="data/ttTi5605dna.RData")

ttTi5605gr<-nanodsmf::DNAstringset2gr(ttTi5605dna)
save(ttTi5605gr,file="data/ttTi5605gr.RData")

