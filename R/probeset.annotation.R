library(biomaRt)

#' Maps probe names to gene names and descriptions
#'
#' @export
probeset.annotation <- function
(
 probesets,
 dataset="hsapiens_gene_ensembl",
 mapping="ensg",
 symbol.idx="hgnc_symbol",
 na.rm=F
)
{
  mart <- useMart(biomart="ensembl", dataset=dataset)
  probesetsMap <- gsub("_at","",probesets)
  GS <- {
    if ( grepl("entrez",mapping) ) 
      getGene( id=probesetsMap, type="entrezgene", mart=mart)
    else if ( grepl( "ens",mapping) )
      getGene( id=probesetsMap, type="ensembl_gene_id", mart=mart)
    else {
      getGene( id=probesetsMap, type=mapping, mart=mart)
    }
  }
  ann <- data.frame(symbol=probesets,description="No information for gene",stringsAsFactors=F)
  rownames(ann) <- probesets

  ## use probeIDs for missing gene symbols and 'No info' for missing descriptions
  ##
  sNA.idx <- is.na(GS[,symbol.idx]) | GS[,symbol.idx]=="" # missing gene symbols
  dNA.idx <- is.na(GS[,"description"]) | GS[,"description"]=="" # missing descriptions
  GS[sNA.idx,symbol.idx] <- GS[sNA.idx,1]
  GS[dNA.idx,"description"] <- "No information for gene"

  verbose("\t",sum(sNA.idx)," probesets w/o gene symbol\n",sep="")
  verbose("\t",sum(dNA.idx)," probesets w/o description\n",sep="")

  ## abort if less than 30% probesets match to gene symbols
  ##
  if ( sum(!is.na(match.idx <- match(probesetsMap,GS[,1])))<round(length(probesets)*.3) )
    stop( "less than 30% probesets annotated, something wrong!" )

  ## create map with same order as input (unless na.rm)
  ##
  ann[!is.na(match.idx),] <- GS[match.idx[!is.na(match.idx)],c(symbol.idx,"description")]
  if ( na.rm ) ann <- ann[!is.na(match.idx),,drop=F]

  if ( any(colnames(ann)!=c("symbol","description")) )
    stop( "colnames(ann)!=c(symbol,description)" )

  ann
}
annotate.data <- function
(
 dat,
 release=c("current","2012"),
 probe.postfix="",
 map=FALSE,
 merge=TRUE
 )
{
  ## find mapping to gene symbols and rename dataset accordingly (to
  ## make things more fun we do the mapping from rat ensembl to hgnc
  ## gene symbols)

  release <- match.arg(release)
  
  require(biomaRt)
  
  if ( release=="current" ) {
    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    rat <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  }
  else if ( release=="2012" ) {
    ## WARNING: biomart release 69 has a known bug and will be fixed in Jan 2013, using the release of May 2012:
    human <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",host="may2012.archive.ensembl.org")
    rat <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "rnorvegicus_gene_ensembl",host="may2012.archive.ensembl.org")
  }
  else {
    stop( "unrecognized release:", release )
  }
  ## establish mapping between Ensemble ID in rats and Gene Symbol in human
  ratList <- getLDS(attributes=c("ensembl_gene_id"),
                    filters="ensembl_gene_id",
                    values=gsub(probe.postfix,"",rownames(dat)),
                    mart = rat,
                    attributesL=c("hgnc_symbol","ensembl_gene_id","description"),
                    martL=human)
  IDX <- grep('Ensembl.Gene.ID',colnames(ratList)); if (length(IDX)!=2) stop( "two Ensembl.Gene.ID's expected: ", colnames(ratList) )
  colnames(ratList)[IDX] <- c('Ensembl.Gene.ID.Rat','Ensembl.Gene.ID.Human')

  ## eliminate unmapped entries, and keep only first of each many-to-one entries
  ratList <- ratList[ratList[,2]!='',]
  ratList <- ratList[match(unique(ratList[,1]),ratList[,1]),]

  if ( any(is.na(match.idx <- match(ratList[,1],sub(paste(probe.postfix,"$",sep=""),'',rownames(dat))))) )
    stop( "missing probes from dat" )

  if (map) {
    ratList[,1] <- paste(ratList[,1],probe.postfix,sep="")
    return(ratList)
  }
  ## ELSE ..
  dat <- dat[match.idx,,drop=F]
  if ( merge ) {
    verbose("merging multiple probes ..")
    dat <- map.many2one(dat,map=cbind(from=paste(ratList[,1],probe.postfix,sep=""),to=ratList[,2]))$dat
    ratList <- ratList[match(rownames(dat),ratList[,2]),]
    verbose("done.\n")
  }
  else
    rownames(dat) <- ratList[,2]

  list(dat=dat,map=ratList)
}
