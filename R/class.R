
setClass('circObj',
         slots = c(
           samples='vector',
           nsample='numeric',
           nsample.g1 = 'numeric',
           nsample.g2 = 'numeric',
           circRNA = 'list',
           circ.method='character',
           species='character',
           circRNA.all = "GRanges",
           circRNA.all.uniq = "data.frame",
           circRNA.DE = 'data.frame',
           A5BS.cluster = "data.frame",
           A3BS.cluster = "data.frame",
           circCluster.DE = 'list',
           metadata = "data.frame"
         )
)

setMethod(
  f = 'show',
  signature = 'circObj',
  definition = function(object) {
    cat(
      'An object of class',
      class(object),
      'with circRNAs predicted by',
      object@circ.method,
      '\n',
      'The object consists of ', length(object@circRNA),'\n'
    )
    for(i in 1:length(object@circRNA)){
      cat(
        object@samples[i],
        ': There are ',length(x = object@circRNA[[i]]),' circRNAs'
        #     length(x = object@ABS5[[i]]),
        #     ' ABS5; ',
        #     length(x = object@ABS3[[i]]),
        #     ' ABS3.',
        #     '\n'
      )
    }
    invisible(x = NULL)
  }
)

setGeneric("samples", function(x) standardGeneric("samples"))
setGeneric("samples<-", function(x, value) standardGeneric("samples<-"))
setMethod('samples', 'circObj', function(x) x@samples)
setMethod('samples<-', 'circObj', function(x, value) {
  x@samples <- value
  x
})

setGeneric("samples", function(x) standardGeneric("samples"))
setGeneric("samples<-", function(x, value) standardGeneric("samples<-"))
setMethod('samples', 'circObj', function(x) x@samples)
setMethod('samples<-', 'circObj', function(x, value) {
  x@samples <- value
  x
})

setGeneric("circRNA", function(x) standardGeneric("circRNA"))
setGeneric("circRNA<-", function(x, value) standardGeneric("circRNA<-"))
setMethod('circRNA', 'circObj', function(x) x@circRNA)
setMethod('circRNA<-', 'circObj', function(x, value) {
  x@circRNA <- value
  x
})

setGeneric("circRNAFlankingIntron", function(x) standardGeneric("circRNAFlankingIntron"))
setGeneric("circRNAFlankingIntron<-", function(x, value) standardGeneric("circRNAFlankingIntron<-"))
setMethod('circRNAFlankingIntron', 'circObj', function(x) x@circRNAFlankingIntron)
setMethod('circRNAFlankingIntron<-', 'circObj', function(x, value) {
  x@circRNAFlankingIntron <- value
  x
})

setGeneric("circClus", function(x) standardGeneric("circClus"))
setGeneric("circClus<-", function(x, value) standardGeneric("circClus<-"))
setMethod('circClus', 'circObj', function(x) x@circClus)
setMethod('circClus<-', 'circObj', function(x, value) {
  x@circClus <- value
  x
})

setGeneric("leftFlankingIntron", function(x) standardGeneric("leftFlankingIntron"))
setGeneric("leftFlankingIntron<-", function(x, value) standardGeneric("leftFlankingIntron<-"))
setMethod('leftFlankingIntron', 'circObj', function(x) x@leftFlankingIntron)
setMethod('leftFlankingIntron<-', 'circObj', function(x, value) {
  x@leftFlankingIntron <- value
  x
})

setGeneric("rightFlankingIntron", function(x) standardGeneric("rightFlankingIntron"))
setGeneric("rightFlankingIntron<-", function(x, value) standardGeneric("rightFlankingIntron<-"))
setMethod('rightFlankingIntron', 'circObj', function(x) x@rightFlankingIntron)
setMethod('rightFlankingIntron<-', 'circObj', function(x, value) {
  x@rightFlankingIntron <- value
  x
})

setGeneric("circRNA.all", function(x) standardGeneric("circRNA.all"))
setMethod("circRNA.all", "circObj", function(x) {x@circRNA.all})
setGeneric("circRNA.all<-", function(x, value) standardGeneric("circRNA.all<-"))
setMethod("circRNA.all<-", "circObj", function(x, value) {
  x@circRNA.all <- value
  x
})

setGeneric("circRNA.all.uniq", function(x) standardGeneric("circRNA.all.uniq"))
setMethod("circRNA.all.uniq", "circObj", function(x) {x@circRNA.all.uniq})
setGeneric("circRNA.all.uniq<-", function(x, value) standardGeneric("circRNA.all.uniq<-"))
setMethod("circRNA.all.uniq<-", "circObj", function(x, value) {
  x@circRNA.all.uniq <- value
  x
})

