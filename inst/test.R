library(devtools)


devtools::load_all()  # Load the updated package
devtools::document()  # Update documentation
devtools::check()     # Run package checks

circexplorers=c(
  '/Users/fengdi/Library/CloudStorage/OneDrive-UniversityofFlorida/Collaborations/circRNA_project/Rpackage/circRNA/cerebellum/CIRCexplorer/SRR3192427.txt',
  '/Users/fengdi/Library/CloudStorage/OneDrive-UniversityofFlorida/Collaborations/circRNA_project/Rpackage/circRNA/cerebellum/CIRCexplorer/SRR3192428.txt',
  '/Users/fengdi/Library/CloudStorage/OneDrive-UniversityofFlorida/Collaborations/circRNA_project/Rpackage/circRNA/frontalcortex/CIRCexplorer/SRR3192424.txt',
  '/Users/fengdi/Library/CloudStorage/OneDrive-UniversityofFlorida/Collaborations/circRNA_project/Rpackage/circRNA/frontalcortex/CIRCexplorer/SRR3192425.txt'
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

# # change lisence
# usethis::use_mit_license("Fengdi Zhao")


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

#
# # make data samller
# load("data/BM10.circs.rda")
# # modify as needed
# sampleid = sample(c(1:length(BM10.circs)),100000)
# demo.circs <- BM10.circs[sampleid]  # Example: keep first 100 rows
#
# usethis::use_data(demo.circs)
#
