library(devtools)


devtools::load_all()  # Load the updated package
devtools::document()  # Update documentation
devtools::check()     # Run package checks

circexplorers=c(
  '../cerebellum/CIRCexplorer/SRR3192427.txt',
  '../cerebellum/CIRCexplorer/SRR3192428.txt',
  '../frontalcortex/CIRCexplorer/SRR3192424.txt',
  '../frontalcortex/CIRCexplorer/SRR3192425.txt'
)

circ.obj <- makecircObj(
  samplefiles = circexplorers,
  conditions = c(2, 2),
  circ.method = 'CIRCexplorer2',
  cutoff = 2
)

circ.obj <- getCircCluster(circObj = circ.obj)
circ.obj <- circRNADE(circObj = circ.obj, DEmethod = 'Pois')
circ.obj <- circClusterDE(circObj = circ.obj, circ.cutoff = 2, DEmethod = 'Meta')

A5BS_DE =  getClusterDE(circ.obj, cluster_type = "A5BS")
head(A5BS_DE)
A3BS_DE = getClusterDE(circ.obj, cluster_type = "A3BS")
head(A3BS_DE)

# Step 2: Drill down into specific interesting clusters
getClusterAnnotation(circ.obj, cluster_type = "A5BS", juncid = c(6, 14, 17))



devtools::load_all()  # Load the updated package
devtools::document()  # Update documentation
devtools::check()     # Run package checks


data("demo.circs", package = "circMeta2")
data("metainfo", package = "circMeta2")

circ.obj = makecircObjfromGRanges(GRanges = demo.circs, metadata=metainfo)

circ.obj = getCircCluster(circObj = circ.obj)

circ.obj = circRNADE(circObj = circ.obj, DEmethod = 'GLM', formula_str = "readNumber ~ condid + age + sex")

circ.obj = circClusterDE(circObj=circ.obj, circ.cutoff=2, DEmethod='Meta')

A5BS_DE =  getClusterDE(circ.obj, cluster_type = "A5BS")
head(A5BS_DE)
A3BS_DE = getClusterDE(circ.obj, cluster_type = "A3BS")
head(A3BS_DE)

# Step 2: Drill down into specific interesting clusters
getClusterAnnotation(circ.obj, cluster_type = "A5BS", juncid = c(6, 14, 17))




