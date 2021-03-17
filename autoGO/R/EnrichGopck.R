  ##########################################################################################
# Data 2020-06-03
# Authot: Luis Soares (ldiass@live.com)
##########################################################################################

## listspecies
##' @return a data.frame containing the name of species and the code used on autoGO
##' @export
listspecies <- function() {
  species.df <-
    data.frame(
      code = c('hsa', 'mmu', 'rno', 'cfa'),
      species = c(
        'Homo sapiens',
        'Mus musculus',
        'Rattus norvegicus',
        'Canis lupus familiaris'
      )
    )
  return(species.df)
}


## geneIDAdvise
## Intern function to format the print output
geneIDAdvise <-
  function(x) {
    print(paste("Got gene ID on column:", x))
  }


## automatic_GO_enrich
##' The main function of the package that will iterate the list of file names to run the GO analysis
##' @param x the data.frame or vector that will contain the name/path of files of differentially expressed genes
##' @param spcode a char of the species code of the gene entries. Call listspecies() to look for the species code available
##' @param keytype a char for the key type of genes. The possible values are the following: ENSEMBL, ENTREZID, REFSEQ, SYMBOL, ENSEMBLTRANS, ENSEMBLPROT. The default value is ENSEMBL.
##' @param orderby char containin name of column for which the gene keys must be ordered. If left empty, the order of the input list will be used.
##' @param dotplotgenes integer containing the number of genes to the dotplot, default = 30
##' @param fontSize integer containing the point of font of the plot, default = 8
##' @param genekeyPos integer containing the number of column of the gene identifiers. If no number set, we will try to find it accordingly to the keytype.
##' @param maxnumberofgenes integer containing the maximum number of genes to be evaluated by the GO algorithms. The default value is 5000.
##' @param GO_pvalue number containg the threshold for the GO enrichment algorithms. The default value is 0.05.
##' @param KEGG, GOALL, GOBP, GOMF, GOCC, writeTable and writeplot are logical values. They determine wheter enrichment will be runned and if the dotplot and table will be writed. The default of all values is TRUE.
##' @export
##' @author Soares, Luis

automatic_GO_enrich <-
  function(x = dir(),
           spcode="",
           keytype = "ENSEMBL",
           orderby="",
           dotplotgenes = 30,
           fontSize = 8,
           genekeyPos = 0,
           maxnumberofgenes = 5000,
           GO_pvalue = 0.05,
           KEGG = TRUE,
           GOALL = TRUE,
           GOBP = TRUE,
           GOMF = TRUE,
           GOCC = TRUE,
           writeTable = TRUE,
           writeplot = TRUE,
           tableseps=",") {
    #Confirm filelist type
    filelist<-x
    filelist <- as.character(as.vector(filelist))
    if (!(typeof(filelist) == "character")) {
      stop(
        "It was not possible to read the file list, if let empty, all the files on the wd will be use. Type getwd() to check it."
      )
    }

    #Check key type
    keytype <- toupper(keytype)
    if (match(keytype,
              c('ENSEMBL', 'ENTREZID', 'REFSEQ', 'SYMBOL', 'ENSEMBLTRANS', 'ENSEMBLPROT')) == 0) {
      stop("The key type must be Ensembl, EntrezID, RefSeq, Symbol, EnsemblTrans or EnsemblProt")
    }

    #Check integers
    dotplotgenes = as.integer(dotplotgenes)
    fontSize = as.integer(fontSize)
    genekeyPos = as.integer(genekeyPos)
    maxnumberofgenes = as.integer(maxnumberofgenes)

    if (!(is.integer(
      c(
        dotplotgenes,
        fontSize,
        genekeyPos,
        maxnumberofgenes,
        maxnumberofgenes
      )
    ))) {
      stop(
        "All the following parameters must be integer:
           dotplotgenes, fontSize, genekeyPos and maxnumberofgenes"
      )
    }

    #Checkspecies
    if (c('mmu') == c(spcode)) {
      if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
        BiocManager::install("org.Mm.eg.db")
      }
      library(org.Mm.eg.db)
      myOrgDb <- org.Mm.eg.db
    }

    else if (c('hsa') == c(spcode)) {
      if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        BiocManager::install("org.Hs.eg.db")
      }
      library(org.Hs.eg.db)
      myOrgDb <- org.Hs.eg.db
    }

    else if (c(spcode) == c('rno')) {
      if (!requireNamespace("org.Rn.eg.db", quietly = TRUE)) {
        BiocManager::install("org.Rn.eg.db")
      }
      library(org.Rn.eg.db)
      myOrgDb <- org.Rn.eg.db
    }

    else if (c(spcode) == c('cfa')) {
      if (!requireNamespace("org.Cf.eg.db", quietly = TRUE)) {
        BiocManager::install("org.Cf.eg.db")
      }
      library(org.Cf.eg.db)
      myOrgDb <- org.Cf.eg.db
    }
    else{
      stop("Species not found, type listspecies() to show all species codes available.")
    }

    #Create the gene annotation for the specific species
    keys <- AnnotationDbi::keys(myOrgDb)
    geneAnotation <-
      AnnotationDbi::select(
        myOrgDb,
        keys = keys,
        keytype = 'ENTREZID',
        columns = c('ENTREZID', keytype)
      )

    for (filename in filelist) {
      tryCatch({
      print(paste("File loaded:", filename))
      #Input datafile
      DEGTable <- read.delim(filename, stringsAsFactors = FALSE, sep=tableseps)

      #Check if the table must be reordered
      if(!orderby==""){
      #Find number of orderby column
      orderbycol <-
        grep(orderby, colnames(DEGTable), ignore.case = TRUE)

      #Order the table
      DEGTable <- DEGTable[order(DEGTable[, orderbycol]),]
    }

      #Cut on the maximum number of genes
      if (maxnumberofgenes < nrow(DEGTable)) {
        DEGTable <- DEGTable[1:maxnumberofgenes,]
      }

      #Find the column of Gene IDs on the dataset
      if (!(genekeyPos == 0)) {
        geneIDAdvise(genekeyPos)
      } else if (length(grep("ens", keytype, ignore.case = TRUE)) > 0) {
        genekeyPos = as.numeric(grep("ens", colnames(DEGTable), ignore.case = TRUE))
        geneIDAdvise(genekeyPos)
      } else if (keytype == 'ENTREZID' &&
                 length(grep("entrez", colnames(DEGTable), ignore.case = TRUE)) > 0) {
        genekeyPos = as.numeric(grep("entrez", colnames(DEGTable)))
        geneIDAdvise(genekeyPos)
      } else if (keytype == 'SYMBOL' &&
                 length(grep("symbol", colnames(DEGTable), ignore.case = TRUE)) > 0) {
        genekeyPos = as.numeric(grep("symbol", colnames(DEGTable)))
        geneIDAdvise(genekeyPos)
      } else if (keytype == 'REFSEQ' &&
                 length(grep("refseq", colnames(DEGTable), ignore.case = TRUE)) > 0) {
        genekeyPos = as.numeric(grep("refseq", colnames(DEGTable)))
        geneIDAdvise(genekeyPos)
      }
      else{
        stop("Not possible to find the column of geneID automatically, please set it manually with the parameter genekeyPos")
      }

      #Get the Entrez ID
      DEGTable$Entrez <-
        geneAnotation$ENTREZID[match(DEGTable[, genekeyPos], geneAnotation[, 2])]

      #Create a vector of Entrez ID values, removing NA and duplicates
      genes <-
        DEGTable$Entrez[!(is.na(DEGTable$Entrez))]
      print(paste(nrow(genes)," genes mapped"))

      genes <- dplyr::distinct(as.data.frame(genes), .keep_all = TRUE)

      if (KEGG) {
        print("Running KEGG")
        #Run Keggprofile
        keggResult <-
          KEGGprofile::find_enriched_pathway(
            as.matrix(genes),
            species = spcode,
            download_latest = TRUE,
            returned_pvalue = GO_pvalue
          )
        print(paste(nrow(keggResult$stastic), " KEGG pathways found"))
        if(nrow(keggResult$stastic)==0){
          print("No enriched KEGG pathways found") 
          }else{
          print(paste(nrow(keggResult$stastic), " KEGG pathways found"))
          keggResult$stastic <-
          keggResult$stastic[order(keggResult$stastic$pvalueAdj),]

        #Adding species to the ID
        keggResult$stastic$KeggID <-
          unlist(lapply(rownames(keggResult$stastic), function(x) {
            paste(spcode, x, sep = "")
          }))

        #Create the column where the geneList will be stored
        keggResult$stastic$GeneList <- 'empty'

        #Add the gene list to the static table
        GOsets <- keggResult$detail
        for (name in names(GOsets)) {
          #Getting back the original key types
          keyColumn <- geneAnotation[, 2]
          geneIDoriginal <-
            keyColumn[match(GOsets[[name]], geneAnotation$ENTREZID)]

          keggResult$stastic$GeneList[grep(name, keggResult$stastic$KeggID)] <-
            (paste(geneIDoriginal, collapse = "/ "))
        }
        if (writeTable) {
          #Write the table
          write.table(
            keggResult$stastic,
            file = paste(filename, "_keggGO.txt", sep = ""),
            sep = "\t"
          )
        }

        if(writeplot){
          #Ordering the result by Padj
          mydata <-
            keggResult$stastic[order(keggResult$stastic$pvalueAdj),]

          #Converting the Padj to -log(Padj)
          mydata$TransP <-
            as.numeric(lapply(mydata$pvalueAdj, function(x) {
              -log10(x)
            }))

          #Setting the infinite numbers as the other minimals
          InfiniteValues <- is.infinite(mydata$TransP)
          maxtransp <- max(mydata$TransP[!InfiniteValues])
          mydata$TransP[InfiniteValues] <- maxtransp

          #Write the dotplot
          pdfname <- paste(filename, "_kegg.pdf", sep = "")
          pdf(pdfname)

          mydata$Pathway_Name <-
            factor(mydata$Pathway_Name, levels = mydata$Pathway_Name[order(mydata$TransP)])
          ggplot2::theme_set(
            ggplot2::theme_linedraw() +
              ggplot2::theme(legend.position = "right")
          )
          myplot <-
            ggplot2::ggplot(
              mydata[1:dotplotgenes,],
              ggplot2::aes_string(
                x = 'Percentage',
                y = "Pathway_Name",
                size = 'Gene_Found',
                color = 'TransP'
              )
            ) +
            ggplot2::geom_point() +
            ggplot2::scale_color_continuous(
              low = "red",
              high = "blue",
              name = '-log(FDR)',
              guide = ggplot2::guide_colorbar(reverse = TRUE)
            ) +
            ggplot2::ylab(NULL) + ggplot2::theme(axis.text.y = ggplot2::element_text(
              size = 12,
              angle = 0,
              hjust = 1,
              vjust = 0,
              face = "plain"
            ))   + ggplot2::scale_size(range = c(3, 8))
          print(myplot)

          dev.off()
        }
        }
        print("KEGG done!")
      }

      #Convert gene list to vector
      genes <- unlist(genes)

      if (GOALL) {
        print("Runnnig enrichGO all")

        ego <- clusterProfiler::enrichGO(
          gene = genes,
          keyType = "ENTREZID",
          OrgDb = myOrgDb,
          ont = 'ALL',
          pAdjustMethod = 'BH',
          pvalueCutoff = GO_pvalue,
          qvalueCutoff = 0.1
        )
        if (writeTable) {
          write.table(ego,
                      file = paste(filename, "_GO.csv", sep = ""),
                      sep = "\t")
        }
        if (writeplot) {
          pdf(paste(filename, "_GO_ALL.pdf", sep = ""))
          myplot <-
            clusterProfiler::dotplot(ego, showCategory = dotplotgenes, font.size = 8)
          print(myplot)
          dev.off()
        }
        print("enrichGO all done")
      }

      if (GOBP && writeplot) {
        print("Running enrichGO BP dotplot")
        egobp <- clusterProfiler::enrichGO(
          gene = genes,
          OrgDb = myOrgDb,
          ont = 'BP',
          pAdjustMethod = 'BH',
          pvalueCutoff = GO_pvalue,
          qvalueCutoff = 0.1
        )

        print("enrichGO BP done")
        egobp <-
          clusterProfiler::simplify(egobp,
                   cutoff = 0.7,
                   by = "p.adjust",
                   select_fun = min)
        print("enrichGO BP simplify done")
        pdf(paste(filename, "_bp.pdf", sep = ""))
        myplot <-
          clusterProfiler::dotplot(
            egobp,
            showCategory = dotplotgenes,
            font.size = fontSize
          )
        print(myplot)
        dev.off()
      }

      if (GOMF && writeplot) {
        print("Running enrichGO MF dotplot")
        egomf <- clusterProfiler::enrichGO(
          gene = genes,
          OrgDb = myOrgDb,
          ont = 'MF',
          pAdjustMethod = 'BH',
          pvalueCutoff = GO_pvalue,
          qvalueCutoff = 0.1
        )
        print("enrichGO MF done")
        egomf <-
          clusterProfiler::simplify(egomf,
                   cutoff = 0.7,
                   by = "p.adjust",
                   select_fun = min)
        print("enrichGO MF simplify done")
        pdf(paste(filename, "_mf.pdf", sep = ""))
        myplot <-
          clusterProfiler::dotplot(
            egomf,
            showCategory = dotplotgenes,
            font.size = fontSize
          )
        print(myplot)
        dev.off()
      }

      if (GOCC && writeplot) {
        egocc <- clusterProfiler::enrichGO(
          gene = genes,
          OrgDb = myOrgDb,
          ont = 'CC',
          pAdjustMethod = 'BH',
          pvalueCutoff = GO_pvalue,
          qvalueCutoff = 0.1
        )
        print("enrichGO CC done")
        egocc <-
          clusterProfiler::simplify(egocc,
                   cutoff = 0.7,
                   by = "p.adjust",
                   select_fun = min)
        print("enrichGO CC simplify done")
        pdf(paste(filename, "_cc.pdf", sep = ""))
        myplot <-
          clusterProfiler::dotplot(
            egocc,
            showCategory = dotplotgenes,
            font.size = fontSize
          )
        print(myplot)
        dev.off()
      }
      #End of try catch
      },
      error=function(x){cat("ERROR with file ",filename," : ", conditionMessage(x), "\n")})
      }
    }



