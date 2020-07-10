source("Step0_setup.R")
seu = readRDS(seurat_obj_rds)
Idents(seu) = get_cluster_rename()[Idents(seu)]
library(patchwork)
vln_plots = lapply(seq_len(nrow(cluster_identities)), function(i){
    name = paste(cluster_identities$seurat_clusters[i], ":", cluster_identities$bio_name[i])
    p = Seurat::VlnPlot(seu, signature_genes[[i]], cols = get_clusters_colors(), pt.size = .3)
    p + plot_annotation(title = name)
})
names(vln_plots) = cluster_identities$bio_name
pdf("violins.pdf", width = 10, height = 7)
lapply(vln_plots, function(p){plot(p); NULL})
dev.off()
