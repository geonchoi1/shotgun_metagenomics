#!/usr/bin/env Rscript
# MAG circular tree (final decided layout).
# inside->out: ML tree | MAG-ID filled with phylum colour (label bg) | MIMAG shape (HQ/NC/MQ)
#              | TPM heatmap sqrt single-hue (genome-mode CoverM) | ARG / VFG rings. genus/order->supp.
suppressMessages({
  library(ggtree); library(ggtreeExtra); library(ggplot2); library(treeio)
  library(dplyr); library(ggnewscale); library(scales)
})
base <- "/home/gchoi/wwtp_plasmidome"; fp <- file.path(base,"analysis/mag")
tree <- read.tree(file.path(fp,"inputs/tree.nwk"))
meta <- read.delim(file.path(fp,"inputs/tree_meta.tsv"))
tpm  <- read.delim(file.path(fp,"inputs/tpm_long.tsv"))
av   <- read.delim(file.path(fp,"inputs/arg_vfg_per_mag.tsv"))
meta <- meta[meta$label %in% tree$tip.label,]
meta$phylum <- sub("_[A-Z]$","",meta$phylum)   # GTDB suffix 병합 (Myxococcota_A -> Myxococcota)
tpm  <- tpm[tpm$label %in% tree$tip.label,]
av   <- av[av$label %in% tree$tip.label,]
tpm$sample <- factor(tpm$sample, levels=c("IN","Anaerobic","Anoxic","Oxic","RAS","EF"))
levels(tpm$sample) <- c("Inf","Ana","Anx","Oxi","RAS","Eff")
tpm$tpm_log <- sqrt(tpm$tpm)
meta$mimag3 <- factor(meta$mimag3, levels=c("HQ","NC","MQ"))
qdf <- data.frame(label=meta$label, metric="MIMAG",
                  qual=factor(meta$mimag3, levels=c("HQ","NC","MQ")))

# ARG / VFG rings (per-MAG gene counts)
arg_d <- data.frame(label=av$label, metric="ARG", val=av$ARG)
vfg_d <- data.frame(label=av$label, metric="VFG", val=av$VFG)

ph_order <- names(sort(table(meta$phylum), decreasing=TRUE))
pal <- c("#e6ab02","#e7298a","#1b9e77","#7570b3","#d95f02","#66a61e","#a6761d",
         "#1f78b4","#fb9a99","#b2df8a","#cab2d6","#fdbf6f","#8dd3c7","#bc80bd",
         "#999999","#80b1d3","#bebada","#ffed6f")[seq_along(ph_order)]
names(pal) <- ph_order
meta$phylum <- factor(meta$phylum, levels=ph_order)

p <- ggtree(tree, layout="fan", open.angle=16, size=0.22) %<+% meta

# (1) MAG-ID label filled with phylum colour (text black, background = phylum)
p <- p + geom_tiplab(aes(label=mag_short, fill=phylum), geom="label",
                     size=0.9, label.padding=unit(0.4,"pt"),
                     offset=0.01, color="grey10") +
  scale_fill_manual(values=pal, name="Phylum", na.value="grey85",
                    guide=guide_legend(order=1, ncol=1, keywidth=0.45, keyheight=0.45))

# (2) MIMAG quality — colorstrip tile (HQ green / NC blue / MQ grey)
p <- p + new_scale_fill() +
  geom_fruit(data=qdf, geom=geom_tile, mapping=aes(y=label, x=metric, fill=qual),
             width=0.04, offset=0.14, color="white", size=0.04) +
  scale_fill_manual(values=c(HQ="#1b7837", NC="#2166ac", MQ="#dddddd"), name="MIMAG",
                    guide=guide_legend(order=2, keywidth=0.45, keyheight=0.45))

# (3) TPM heatmap — log10(x+1), single-hue
p <- p + new_scale_fill() +
  geom_fruit(data=tpm, geom=geom_tile, mapping=aes(y=label, x=sample, fill=tpm_log),
             width=0.05, offset=0.04, color="white", size=0.04,
             axis.params=list(axis="x", text.size=1.7, text.angle=-90, line.size=0)) +
  scale_fill_gradient(low="#f7f7f7", high="#000000", na.value="white",
                      limits=c(0,158), oob=scales::squish, name="TPM\n(sqrt)",
                      guide=guide_colorbar(order=3, barwidth=0.6, barheight=4))

# (4) ARG ring (red gradient)
p <- p + new_scale_fill() +
  geom_fruit(data=arg_d, geom=geom_tile, mapping=aes(y=label, x=metric, fill=val),
             width=0.045, offset=0.05, color="white", size=0.04,
             axis.params=list(axis="x", text.size=1.6, text.angle=-90, line.size=0)) +
  scale_fill_gradient(low="#fff5f0", high="#cb181d", na.value="grey92", name="ARG (n)",
                      guide=guide_colorbar(order=4, barwidth=0.6, barheight=3.5))
# (5) VFG ring (purple gradient, sqrt)
p <- p + new_scale_fill() +
  geom_fruit(data=vfg_d, geom=geom_tile, mapping=aes(y=label, x=metric, fill=sqrt(val)),
             width=0.045, offset=0.04, color="white", size=0.04,
             axis.params=list(axis="x", text.size=1.6, text.angle=-90, line.size=0)) +
  scale_fill_gradient(low="#fcfbfd", high="#6a51a3", na.value="grey92", name="VFG (sqrt n)",
                      guide=guide_colorbar(order=5, barwidth=0.6, barheight=3.5))

p <- p + theme(legend.position="right", legend.text=element_text(size=6),
               legend.title=element_text(size=7.5), legend.key.size=unit(0.32,"cm"),
               plot.title=element_text(size=12),
               plot.subtitle=element_text(size=7, color="grey40")) +
  labs(title="MAG phylogenomic tree (WWTP community context)",
       subtitle=paste0(length(tree$tip.label),
         " MAGs (bac120 IQ-TREE ML (Q.YEAST+F+R5, UFBoot1000)) | label bg = phylum | shape = MIMAG (HQ/NC/MQ) | ",
         "rings = relative abundance (%) IN->EF | outer = ARG / VFG gene counts. genus/order in Suppl."))

ggsave(file.path(fp,"mag_tree.png"), p, width=19, height=17, dpi=300, bg="white", limitsize=FALSE)
cat("saved analysis/mag/mag_tree.png | tips",length(tree$tip.label),"\n")
