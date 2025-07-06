#' handle column names of either CIRCexplorer2 or CIRI2
#'
#' @param circ.method specify circRNA calling method
#'
#' @return a vector of column names
#' @export
circNames<-function(circ.method=c('CIRCexplorer2','CIRI2')){
  circ.method=match.arg(circ.method)
  circexplorer2.names=c('chrom','start','end','name','score','strand','thickStart','thinkEnd','itemRgb','exonCount',
                        'exonSizes','exonOffsets','readNumber','circType','genename','isoformName','exonIndex/intronIndex','flankIntron')

  ciri2.names=c('circRNA_ID','chr','circRNA_start','circRNA_end','junction_reads','SM_MS_SMS','non_junction_reads','junction_reads_ratio','circRNA_type','gene_id','strand','junction_reads_ID')

  if(circ.method=='CIRCexplorer2'){
    return(circexplorer2.names)
  }else if(circ.method=='CIRI2'){
    return(ciri2.names)
  }
}

#' create S4 circObj object directly using output files from either CIRCexplorer2 or CIRI2
#'
#' @param samplefiles files path
#' @param conditions specify number of samples for each condition
#' @param circ.method specify method used for calling circRNA
#' @param cutoff  minimal reads for circRNA, default is 2
#'
#' @return circObj S4 object
#' @export
#' @importFrom GenomicRanges GRanges
#' @import methods
#' @importFrom utils read.table
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors Rle mcols
makecircObj<-function(samplefiles, conditions = c(2,2), circ.method=c('CIRCexplorer2','CIRI2'),cutoff=2){
  require(GenomicRanges)

  circ.method=match.arg(circ.method)

  circObj=new(
    Class = "circObj",
    circ.method = circ.method
  )

  nsample=length(samplefiles)
  conds = c(rep(1, conditions[1]),rep(2, conditions[2]))
  circRNA=list()

  if(circ.method=='CIRCexplorer2'){
    for(isample in 1:nsample){
      message("Process ",samplefiles[isample])
      m=read.table(samplefiles[isample],stringsAsFactors  = F)
      colnames(m)=circNames(circ.method)
      m=m[m$circType=='circRNA',]
      m=m[m$chrom!='chrM',]
      m=m[m$readNumber>=cutoff,]
      circ=GRanges(seqnames = Rle(m$chrom), ranges = IRanges(m$start,m$end),strand=m$strand)
      mcols(circ)=m[,c('exonCount','exonSizes','exonOffsets','readNumber','genename','isoformName','exonIndex/intronIndex')]

      leftintrons.chr=rightintrons.chr=character(nrow(m))
      leftintrons.start=leftintrons.end=rightintrons.start=rightintrons.end=numeric(nrow(m))

      introns=unlist(sapply(1:length(m$flankIntron),function(x) strsplit(m$flankIntron[x], split = '\\|')))

      leftintrons=introns[1:length(introns)%%2==1]
      left.rm=which(leftintrons=='None')
      leftintrons.chr[left.rm]='None'
      leftintrons.start[left.rm]= -1
      leftintrons.end[left.rm]= -1
      tmp=unlist(sapply(1:length(leftintrons[-left.rm]),function(x) strsplit(leftintrons[-left.rm][x], split = ':')))
      leftintrons.chr[-left.rm]=tmp[1:length(tmp)%%2==1]
      tmp=tmp[1:length(tmp)%%2==0]
      tmp=unlist(sapply(1:length(tmp),function(x) strsplit(tmp[x], split = '-')))
      leftintrons.start[-left.rm]=as.numeric(tmp[1:length(tmp)%%2==1])
      leftintrons.end[-left.rm]=as.numeric(tmp[1:length(tmp)%%2==0])

      rightintrons=introns[1:length(introns)%%2==0]
      right.rm=which(rightintrons=='None')
      rightintrons.chr[right.rm]='None'
      rightintrons.start[right.rm]= -1
      rightintrons.end[right.rm]= -1
      tmp=unlist(sapply(1:length(rightintrons[-right.rm]),function(x) strsplit(rightintrons[-right.rm][x], split = ':')))
      rightintrons.chr[-right.rm]=tmp[1:length(tmp)%%2==1]
      tmp=tmp[1:length(tmp)%%2==0]
      tmp=unlist(sapply(1:length(tmp),function(x) strsplit(tmp[x], split = '-')))
      rightintrons.start[-right.rm]=as.numeric(tmp[1:length(tmp)%%2==1])
      rightintrons.end[-right.rm]=as.numeric(tmp[1:length(tmp)%%2==0])

      mcols(circ)$leftintrons.chr=leftintrons.chr
      mcols(circ)$leftintrons.start=leftintrons.start
      mcols(circ)$leftintrons.end=leftintrons.end
      mcols(circ)$rightintrons.chr=rightintrons.chr
      mcols(circ)$rightintrons.start=rightintrons.start
      mcols(circ)$rightintrons.end=rightintrons.end

      mcols(circ)$sampleid = isample
      mcols(circ)$condid = conds[isample]

      circRNA[[isample]]=circ
    }

  }else if(circ.method=='CIRI2'){
    for(isample in 1:nsample){
      message("Process ",samplefiles[isample])
      m=read.table(samplefiles[isample],stringsAsFactors  = F,sep='\t',skip=1)
      colnames(m)=circNames(circ.method)
      m=m[,-12]
      m=m[m$chr!='chrM',]
      m=m[m$junction_reads>=cutoff,]
      circ=GRanges(seqnames = Rle(m$chr), ranges = IRanges(m$circRNA_start,m$circRNA_end),strand=m$strand)
      mcols(circ)=m[,c('junction_reads','non_junction_reads','junction_reads_ratio','circRNA_type')]
      mcols(circ)$sampleid = isample
      mcols(circ)$condid = conds[isample]
      circRNA[[isample]]=circ
    }
  }

  circRNA(circObj)=circRNA
  samples(circObj)=samplefiles
  circObj@nsample=nsample
  circObj@nsample.g1 = conditions[1]
  circObj@nsample.g2 = conditions[2]

  circs.all = unlist(GRangesList(circRNA))
  circs.all.uniq= unique(circs.all)

  circid.uniq = seq(length(circs.all.uniq))
  circid=circid.uniq[GenomicRanges::match(circs.all,circs.all.uniq)]
  circs.all$circid = circid

  circs.uniq.df = data.frame(seqnames=circs.all.uniq@seqnames,
                             start=circs.all.uniq@ranges@start,
                             end=end(circs.all.uniq),
                             width=width(circs.all.uniq),
                             strand=circs.all.uniq@strand,
                             gene=circs.all.uniq$genename,
                             circid=circid.uniq)

  circRNA.all(circObj)=circs.all
  circRNA.all.uniq(circObj) = circs.uniq.df
  return(circObj)
}

#' create S4 circObj object directly from circRNA GRanges
#'
#' @param GRanges circRNA in GRanges format
#' @param sparse_filter sparsity filtering of circRNA, based on proportion sample has the circRNA and the mean expression of circRNA.
#' @param metadata metadata information regarding the samples, covariates such as age, gender...
#'
#' @return circObj S4 object
#' @export
#' @import methods
#' @importFrom stats median
makecircObjfromGRanges<-function(GRanges,sparse_filter=0.05,metadata=NULL){

  circObj=new(
    Class = "circObj"
  )

  circObj@samples = unique(GRanges$sampleid)
  circObj@nsample = length(unique(GRanges$sampleid))
  circObj@nsample.g1 = length(unique(GRanges$sampleid[GRanges$condid == 1]))
  circObj@nsample.g2 = length(unique(GRanges$sampleid[GRanges$condid == 2]))

  circs.all = GRanges
  circs.all.uniq= unique(circs.all)

  circid.uniq = seq(length(circs.all.uniq))
  circid=circid.uniq[GenomicRanges::match(circs.all,circs.all.uniq)]
  circs.all$circid = circid
  circs.uniq.df = data.frame(seqnames=circs.all.uniq@seqnames,
                             start=circs.all.uniq@ranges@start,
                             end=end(circs.all.uniq),
                             width=width(circs.all.uniq),
                             strand=circs.all.uniq@strand,
                             gene=circs.all.uniq$genename,
                             circid=circid.uniq)

  # counts matrix for circRNA
  counts = matrix(0,nrow=length(unique(circs.all$circid)),ncol=length(unique(circs.all$sampleid)))
  juncread.all = circs.all$readNumber
  colidx = circs.all$sampleid
  rowidx = circs.all$circid
  colnames(counts) = unique(circs.all$sampleid)
  rownames(counts) = paste0("circ-",unique(circs.all$circid))
  for (i in 1:length(circs.all)){
    #print(i)
    counts[rowidx[i],colidx[i]] = juncread.all[i]
  }

  # drop sparse circRNA based on preset criteria
  # filter circRNA. proportion of sample has this circRNA and mean express level 10%
  xx.m1=rowMeans(counts!=0)
  xx.m2=rowMeans(counts)
  threshold=sparse_filter
  id=which(xx.m2>threshold & xx.m1>threshold) #filter circRNA

  circs.all = circs.all[circs.all$circid %in% id]
  circs.uniq.df = circs.uniq.df[circs.uniq.df$circid %in% id, ]
  counts = counts[id, ]

  # normalize read counts
  sfs=colSums(counts);sfs=sfs/median(sfs)   #normalize
  counts=sweep(counts,2,sfs,FUN='/') #normalize
  counts = ceiling(counts)

  adj.juncread.all = numeric(length(circs.all))
  colidx = circs.all$sampleid
  rowidx = paste0("circ-",circs.all$circid)
  for (i in 1:length(circs.all)){
    # print(i)
    adj.juncread.all[i] = counts[rowidx[i],colidx[i]]
  }

  circs.all$readNumber = adj.juncread.all
  ######################################
  circRNA.all(circObj)=circs.all
  circRNA.all.uniq(circObj) = circs.uniq.df
  circObj@metadata = metadata
  return(circObj)
}


#' clustering individual circRNA based on A5BS and A3BS events
#'
#' @param circObj circObj object
#'
#' @return a updated circObj that contains clustering information under slots A5BS.cluster and A3BS.cluster
#' @export
#' @importFrom GenomicRanges seqnames start end strand
#' @importFrom dplyr %>% group_by
getCircCluster<-function(circObj){
  require(dplyr)
  circs.all = circObj@circRNA.all
  ####
  bsj1=paste(seqnames(circs.all),start(circs.all),sep=':')
  bsj2=paste(seqnames(circs.all),end(circs.all),sep=':')
  plus=which(as.character(strand(circs.all))=='+')
  minus=which(as.character(strand(circs.all))=='-')

  # 5' to 3' start  # + direction and start pos
  bsj51=c(bsj1[plus]); bsj51.pos=c(start(circs.all)[plus])
  # 5' to 3' end  # + direction and end pos
  bsj52=c(bsj2[plus]); bsj52.pos=c(end(circs.all)[plus])
  # 3' to 5' start
  bsj31=c(bsj1[minus]); bsj31.pos=c(start(circs.all)[minus])
  # 3' to 5' end
  bsj32=c(bsj2[minus]); bsj32.pos=c(end(circs.all)[minus])

  #5' to 3' data frame
  A5BS=data.frame(from=bsj51,to=bsj52,from.pos=bsj51.pos,end.pos=bsj52.pos,
                  sampleid=circs.all$sampleid[c(plus)],condid=circs.all$condid[c(plus)],
                  circid=circs.all$circid[c(plus)],nread=circs.all$readNumber[c(plus)],
                  genename=circs.all$genename[plus])

  # dim(A5BS)
  # length(unique(A5BS$circid))
  # sum(A5BS$end.pos>A5BS$from.pos)

  #3' to 5' data frame
  A3BS=data.frame(from=bsj31,to=bsj32,from.pos=bsj31.pos,end.pos=bsj32.pos,
                  sampleid=circs.all$sampleid[c(minus)],condid=circs.all$condid[c(minus)],
                  circid=circs.all$circid[c(minus)],nread=circs.all$readNumber[c(minus)],
                  genename=circs.all$genename[minus])

  # dim(A3BS)
  # sum(A3BS$end.pos>A3BS$from.pos)
  # length(unique(A3BS$circid))

  # clustering based on A5BS event and A3BS event
  A5BS.order=A5BS[order(A5BS$from,A5BS$to,A5BS$sampleid,A5BS$circid),]

  A5BS.clus=A5BS.order %>% group_by(from)
  A5BS.tb=table(A5BS.clus$from)
  juncid=unlist(sapply(1:length(A5BS.tb),function(x) rep(x,A5BS.tb[x])))
  A5BS.clus$juncid=juncid

  A3BS.order=A3BS[order(A3BS$to,A3BS$from,A3BS$sampleid,A3BS$circid),]
  A3BS.clus=A3BS.order %>% group_by(to)
  A3BS.tb=table(A3BS.clus$to)
  juncid=unlist(sapply(1:length(A3BS.tb),function(x) rep(x,A3BS.tb[x])))
  A3BS.clus$juncid=juncid

  circObj@A5BS.cluster = A5BS.clus
  circObj@A3BS.cluster = A3BS.clus
  return(circObj)
}

#' individual circRNA DE method
#'
#' @param dat circRNA by sample counts matrix for circRNA
#'
#' @return DE results for individual circRNA
#' @export
#' @importFrom stats pnorm
runPois.ztest<-function(dat){
  sfs=colSums(dat$counts);sfs=sfs/min(sfs)   #normalize
  dat$counts=sweep(dat$counts,2,sfs,FUN='/') #normalize
  n0=sum(dat$designs==0)
  n1=sum(dat$designs==1)
  m0=rowMeans(dat$counts[,dat$designs==0]) #cerebellum
  m1=rowMeans(dat$counts[,dat$designs==1]) #diencephalon
  n=nrow(dat$counts)
  pval=rep(1,n)
  for(i in 1:n){
    z=(m1[i]-m0[i])/sqrt(m1[i]/n1+m0[i]/n0)
    pval[i]=2*pnorm(-abs(z))
  }
  list(pval=pval,m0=m0+0.01,m1=m1+0.01,fc=(m1+0.01)/(m0+0.01))
}

#' individual circRNA DE function
#'
#' @param circObj circObj contains circRNA.all.uniq and circRNA.all
#' @param DEmethod Pois method for small sample no covariates; GLM for large sample with covariates.
#' @param formula_str specify formular used for poisson GLM regression. e.g. "readNumber ~ condid + age + sex"
#'
#' @return updated circObj that contains inidividual circRNA DE result under circRNA.DE slot
#' @export
#' @importFrom progress progress_bar
#' @importFrom stats p.adjust glm formula poisson coef
circRNADE<-function(circObj,DEmethod=c('Pois', 'GLM', 'edgeR', 'DESeq2'),formula_str=NULL){

  if(DEmethod=='Pois'){

    tmp = as.data.frame(circObj@circRNA.all)
    circid.uni=unique(tmp$circid)
    ncirc=length(circid.uni) # number of unique circRNA
    nsample = circObj@nsample
    nsample.g1 = circObj@nsample.g1
    nsample.g2 = circObj@nsample.g2

    # Vectorized matrix assignment
    m=matrix(0, ncirc, nsample)

    # Create mappings once
    circid_to_row <- setNames(1:ncirc, circid.uni)

    # Get row and column indices for all data at once
    row_indices <- circid_to_row[tmp$circid]
    col_indices <- tmp$sampleid

    # Single vectorized assignment
    m[cbind(row_indices, col_indices)] <- tmp$readNumber

    dat=NULL
    dat$counts=m
    dat$designs= c(rep(0,nsample.g1), rep(1,nsample.g2))
    res=runPois.ztest(dat)  # read counts matrix is normalized in runPois.ztest function
    fdr=p.adjust(res$pval,method='fdr')

    # Convert matrix back to data frame for results
    m=as.data.frame(m)
    m$circid=circid.uni
    m$genename=tmp$genename[match(circid.uni,tmp$circid)]
    m$m0 = res$m0
    m$m1 = res$m1
    m$log2fc=log2(res$fc)
    m$pvalue=res$pval
    m$fdr=fdr
    m$direction = "NULL"
    m[m$log2fc>0,]$direction = "up"
    m[m$log2fc<0,]$direction = "dw"

    circObj@circRNA.DE = m
    return(circObj)
  }

  if(DEmethod=='GLM'){
    tmp <- data.frame(circObj@circRNA.all)
    circid.uni <- unique(tmp$circid)
    ncirc <- length(circid.uni) # number of unique circRNAs
    formula_vars = all.vars(formula(formula_str))
    meta = circObj@metadata
    unique.circRNAs = circObj@circRNA.all.uniq
    res=matrix(0,ncirc,4)

    # Check if all formula variables exist in data
    required_vars <- setdiff(formula_vars, "readNumber")
    if(!all(required_vars %in% colnames(meta))) {
      missing_vars <- setdiff(required_vars, colnames(meta))
      stop("Missing variables in metadata: ", paste(missing_vars, collapse=", "))
    }

    # Pre-split data once
    tmp_by_circid <- split(tmp[,c("sampleid", "readNumber")], tmp$circid)

    # Pre-compute sample mapping
    sample_to_idx <- setNames(seq_len(nrow(meta)), meta$sampleid)

    # Reuse model_data object
    model_data <- meta
    if(!"readNumber" %in% colnames(model_data)) {
      model_data$readNumber <- 0
    }

    # Pre-compile formula object
    formula_obj <- formula(formula_str)

    pb <- progress_bar$new(
      format = "  [:bar] :percent :elapsedfull",
      total = ncirc, clear = FALSE, width = 60
    )

    for(icirc in 1:ncirc){
      circ_id <- as.character(circid.uni[icirc])

      # Get pre-split data instead of filtering entire dataframe
      circ_data <- tmp_by_circid[[circ_id]]

      # Reset readNumber for each circRNA (fresh start)
      model_data$readNumber[] <- 0  # Use for in-place assignment

      # Robust alignment
      if(!is.null(circ_data) && nrow(circ_data) > 0) {
        # Use pre-computed mapping instead of match()
        sample_indices <- sample_to_idx[circ_data$sampleid]
        valid_indices <- !is.na(sample_indices)
        if(any(valid_indices)) {
          model_data$readNumber[sample_indices[valid_indices]] <- circ_data$readNumber[valid_indices]
        }
      }

      # GLM fitting
      tryCatch({
        fit = glm(formula_obj, data=model_data, family = poisson(link = "log"))

        coefficients <- coef(fit)
        m1 <- exp(coefficients["(Intercept)"])

        # Dynamic coefficient detection
        cond_coef_name <- grep("condid", names(coefficients), value=TRUE)[1]
        if(length(cond_coef_name) > 0 && !is.na(coefficients[cond_coef_name])) {
          m2 <- exp(coefficients["(Intercept)"] + coefficients[cond_coef_name])
          coef_summary <- summary(fit)$coefficients
          cond_row <- which(rownames(coef_summary) == cond_coef_name)
          pval <- coef_summary[cond_row, 4]
          log2fc <- coef_summary[cond_row, 1] / log(2)
        } else {
          m2 <- m1
          pval <- 1
          log2fc <- 0
        }

        res[icirc,] = c(pval, m1, m2, log2fc)

      }, error = function(e) {
        warning(paste("GLM failed for circRNA", circ_id, ":", e$message))
        res[icirc,] = c(NA, NA, NA, NA)
      })

      pb$tick()
    }

    # Result processing
    res = data.frame(res)
    colnames(res) = c("pvalue","m0","m1","log2fc")
    res$fdr=p.adjust(res$pvalue,method='fdr')
    res$circid = circid.uni
    res$direction = "NULL"
    res[!is.na(res$log2fc) & res$log2fc>0,]$direction = "up"
    res[!is.na(res$log2fc) & res$log2fc<0,]$direction = "dw"
    res2 = merge(unique.circRNAs, res, by="circid")
    circObj@circRNA.DE = res2
    return(circObj)
  }
}


#' circ-cluster DE function
#'
#' @param circObj updated circObj that contains individual circRNA DE results
#' @param circ.cutoff minimal circRNA number within the circ-cluster
#' @param DEmethod cluster DE methods: Meta-novel statistical method, edgeR, DESeq2
#'
#' @return circ-cluster DE result stored in list, A5BS cluster DE and A3BS cluster DE
#' @export
#' @importFrom progress progress_bar
#' @importFrom stats qnorm pnorm p.adjust
circClusterDE<-function(circObj, circ.cutoff=2, DEmethod=c('Meta', "edgeR", 'DESeq2')){

  cluster_types <- c('A5BS.cluster', 'A3BS.cluster')
  results <- list()

  for (cluster_type in cluster_types) {
    if (cluster_type %in% slotNames(circObj)) {

      print(cluster_type)

      clus.dat <- slot(circObj, cluster_type)
      circRNADE.out <- circObj@circRNA.DE
      circid=unique(clus.dat$circid)
      ncirc=length(circid)
      juncid=unique(clus.dat$juncid)
      njunc=length(juncid)
      ncirc_in_junc=numeric(njunc)
      pvals_junc=numeric(njunc)
      fc = numeric(njunc)
      m0 = numeric(njunc)
      m1 = numeric(njunc)

      if(DEmethod=="Meta"){
        pvals.alls = circRNADE.out$pvalue
        names(pvals.alls) = circRNADE.out$circid
        m0.all = circRNADE.out$m0
        names(m0.all) = circRNADE.out$circid
        m1.all = circRNADE.out$m1
        names(m1.all) = circRNADE.out$circid

        # method combined pvalue from single circRNA
        pb <- progress_bar$new(
          format = "  [:bar] :percent :elapsedfull",
          total = njunc, clear = FALSE, width = 60
        )
        for(ijunc in 1:njunc){
          tmpcircid=unique(clus.dat$circid[clus.dat$juncid==juncid[ijunc]])
          ncirc_in_junc[ijunc]=length(tmpcircid)
          fc[ijunc] = sum(m1.all[match(tmpcircid,names(m1.all))])/
            sum(m0.all[match(tmpcircid,names(m0.all))])
          m0[ijunc] = sum(m0.all[match(tmpcircid,names(m0.all))])
          m1[ijunc] = sum(m1.all[match(tmpcircid,names(m1.all))])

          # determine cluster regulate direction
          if(fc[ijunc] > 1){
            clus.dir = "up"
          } else {
            clus.dir ="dw"
          }

          # cluster contains only 1 circRNA case
          if(length(tmpcircid)==1){
            pvals_junc[ijunc]= pvals.alls[match(tmpcircid,names(pvals.alls))]
          } else {
            # cluster contains more than 1 circRNA case
            clus.counts = numeric(length(tmpcircid))
            for(i in 1:length(tmpcircid)){
              count = sum(clus.dat[clus.dat$circid==tmpcircid[i],"nread"])
              clus.counts[i] = count
            }

            # read counts weight
            log.weight = log(clus.counts)
            log.weight.norm = log.weight/sum(log.weight)

            # directional weight
            dir.weight = NULL
            orig.circ.pval = NULL
            for(i in 1:length(tmpcircid)){
              circ.dir =  circRNADE.out[circRNADE.out$circid==tmpcircid[i],"direction"]
              circ.pval = circRNADE.out[circRNADE.out$circid==tmpcircid[i],"pvalue"]
              if(circ.dir == clus.dir){
                dir.weight = c(dir.weight, 2)
                orig.circ.pval = c(orig.circ.pval,circ.pval)
              } else {
                dir.weight = c(dir.weight, 0.5)
                orig.circ.pval = c(orig.circ.pval,circ.pval)
              }
            }
            #convert pvalue to Z-score; two tails test (testing wether coeff is zero not not)
            z.scores = -qnorm(orig.circ.pval/2)
            comp.weight = log.weight.norm*dir.weight
            weighted.z.scores = z.scores * comp.weight
            # stouffe's method
            weighted.sum.z = sum(weighted.z.scores)/sqrt(sum(comp.weight^2))
            adjusted.pval = 2*pnorm(-abs(weighted.sum.z))
            pvals_junc[ijunc]= adjusted.pval
          }

          pb$tick()
        }
        juncDE = data.frame(juncid = juncid, numcircs =ncirc_in_junc, m0=m0, m1=m1, pvalue = pvals_junc, log2fc=log2(fc))
        juncDE$fdr = p.adjust(juncDE$pvalue,method='fdr')
        juncDE = juncDE[juncDE$numcircs>=circ.cutoff,]


        results[[cluster_type]] <- juncDE
        circObj@circCluster.DE = results

      }
    }
  }
  return(circObj)
}

