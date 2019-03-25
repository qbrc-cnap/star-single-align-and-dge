#library(calibrate)
if(!require("calibrate", character.only=T)) stop("Please install the calibrate package first.")

Draw.dge.Volcano<-function(dge_result=NA,
                       pvalue.thresh=0.05,
                       log2fc.thresh=1.2,
                       title="",
                       outfile=NA,
                       res=150,
                       width=1200,
                       height=1200){

    #for(dataset in names(dge_result)){
    # Make a basic volcano plots
    # outfile=paste0(result_dir, dataset, ".volcano.png")
    if(!(is.na(outfile))){ png(outfile, res=res, width=width, height=height) }

    ##
    ##  Generate volcano plot for a deseq dge result table
    ##  The function assumes the table columns are:
    ##
    dge_res=dge_result

    with(dge_res, plot(log2FoldChange, -log10(pvalue), pch=20,
                       main=title,
                       xlim=c(min(dge_res$log2FoldChange)-1.5, max(dge_res$log2FoldChange)+1.5)),
                       ylim=c(0, max(-log10(pvalue))+1))
    # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
    with(subset(dge_res,  pvalue<pvalue.thresh & abs(log2FoldChange)>log2fc.thresh),
         points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
    with(subset(dge_res, pvalue<pvalue.thresh),
         points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
    # Label points with the textxy function from the calibrate plot
    with(subset(dge_res,  pvalue<pvalue.thresh & abs(log2FoldChange)>log2fc.thresh),
         textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=.8))
    if(!(is.na(outfile))){dev.off()}
}

