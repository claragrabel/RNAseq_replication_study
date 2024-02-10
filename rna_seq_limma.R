####################################################################
## Análisis de los resultados de HISAT2 y STRINGTIE               ##
## usando los paquetes de R ballgown y limma.                     ##  
##                                                                ##
## Autor: Francisco J. Romero-Campero fran@us.es                  ##
####################################################################

## El paquete de bioconductor ballgown proporciona las funciones necesarias para 
## realizar un análisis de expresión génica diferencial y visualizar los resultados
## a partir de procesamiento de los datos brutos de secuenciación realizados con 
## hisat2 y stringtie. 

## Para ejecutar con éxito este script es necesario descargar la carpeta samples
## completa a tu ordenador, mover este script a la carpeta samples y fijar el 
## Working Directory To Source File Location. 

## Instalación y carga de los paquetes necesarios. Sólo es necesario instalar los
## paquetes la primera vez que se ejecuta este script en un ordenador el resto de las
## veces bastará cargar los paquetes simplemente. 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ballgown")

library(ballgown)


## Para cargar los datos es necesario crear previamente un fichero tabular
## que contenga como primera columna los nombres de las carpetas donde se guarda
## cada muesra típicamente sample1, sample2, etc ... El resto de columnas
## haran referencia al genotipo, tratamiento y demás caracteriśticas de cada muestra. 
pheno.data <- read.csv("pheno_data.csv")
pheno.data

## La función ballgown se usa para cargar o leer los datos. Es necesario especificar
## el directorio donde se encuentra las muestras. En nuestro caso especificamos .
## para indicar que se encuentran en el actual directorio. 
bg.data <- ballgown(dataDir = ".", samplePattern = "sample", pData=pheno.data)
bg.data
sampleNames(bg.data)

## La función gexpr extrae los niveles de expresión génicos en cada muestra
## medidos como FPKM (Fragments Per Kilobase of exon and Million of mapped reads)
gene.expression <- gexpr(bg.data)
head(gene.expression)
dim(gene.expression)
gene.names <- rownames(gene.expression)

## Nombramos las columnas con los nombres de nuestras muestras. 
colnames(gene.expression) <- c("col0_1","col0_2","abc_1","abc_2")

## Por motivos técnicos sumamos 1 a todos los niveles de expresión. 
## El problema viene provocado por x < 1 --> log2(x) < 0
gene.expression.1 <- gene.expression + 1

## Guardamos los datos de expresión génica sin procesar
write.table(x = gene.expression.1,file = "pcg_gene_expression.tsv",
            quote = F,row.names = F,
            sep = "\t")

## Representación de la distribución global de la expresión génica
boxplot(gene.expression, outline=F,col=rainbow(4),ylab="Gene Expression (FPKM)",
        cex.lab=1.5)

boxplot(log2(gene.expression.1), outline=F,col=rainbow(4),
        ylab="log2(Gene Expression)",
        cex.lab=1.5)

## En este punto ballgown no realiza ninguna normalización de los datos más alla
## del cálculo de los niveles de expresión por FPKM.
## Utilizamos el paquete de R NormalyzerDE para esta tarea. Para ello es neceario generar
## un fichero con un formato específico.


BiocManager::install("NormalyzerDE")

library(NormalyzerDE)

design <- data.frame(sample=colnames(gene.expression),
                     group=c(rep("col0",2),rep("abc",2)))

write.table(x = design,file = "normalyzer_design.tsv",quote = F,row.names = F,
            sep = "\t")

normalyzer(jobName = "PCG",designPath = "normalyzer_design.tsv",
           dataPath = "pcg_gene_expression.tsv",outputDir = ".")


normalized.gene.expression <- read.table(file="PCG/Quantile-normalized.txt", header=T)
head(normalized.gene.expression)
rownames(normalized.gene.expression) <- gene.names
 
## Representación de la distribución global de la expresión génica tras 
## normalización
boxplot(normalized.gene.expression, outline=F,col=rainbow(4),
        ylab="log2(FPKM + 1)",cex.lab=1.5)

## Previsualizamos la similitud entre las réplicas
plot(x = normalized.gene.expression[,"col0_1"],
     y = normalized.gene.expression[,"col0_2"],
     pch=19,col="grey",xlab="Col0_1",ylab="Col0_2",cex=0.5)
text(x=3,y=14,
     labels = paste(c(
      "cor = ",
      round(100*cor(normalized.gene.expression[,"col0_1"],
                    normalized.gene.expression[,"col0_2"]),
            digits = 2),
      "%"), collapse=""))


plot(x = normalized.gene.expression[,"abc_1"],
     y = normalized.gene.expression[,"abc_2"],
     pch=19,col="grey",xlab="abc_1",ylab="abc_2",cex=0.5)
text(x=3,y=14,
     labels = paste(c(
      "cor = ",
      round(100*cor(normalized.gene.expression[,"abc_1"],
                    normalized.gene.expression[,"abc_2"]),
            digits = 2),
      "%"), collapse=""))

## Realizamos un análisis de componentes principales y clustering 
## jerárquico para continuar con la exploración de los datos.
library(FactoMineR)
library(factoextra)

## Los dos paquetes anteriores son geniales para análisis estadísticos en
## datos de secuencaición de nueva generación. Tienen una documentación muy 
## completa: http://www.sthda.com/english/articles/tag/factominer/
## Documentación genérica: http://www.sthda.com/english/

## Por ejemplo para PCA y clustering jerárquico:
## http://www.sthda.com/english/articles/22-principal-component-methods-videos/65-pca-in-r-using-factominer-quick-scripts-and-videos/
## http://www.sthda.com/english/articles/22-principal-component-methods-videos/74-hcpc-using-factominer-video/

pca.gene.expression <- data.frame(colnames(normalized.gene.expression),
                                  t(normalized.gene.expression))
colnames(pca.gene.expression)[1] <- "Sample"
head(pca.gene.expression)

res.pca <- PCA(pca.gene.expression, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
res.hcpc <- HCPC(res.pca, graph=FALSE,nb.clust = 2)   
fviz_dend(res.hcpc,k=2,
          cex = 0.75,                       # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle",
          labels_track_height = 1400      # Augment the room for labels
)

fviz_pca_ind(res.pca, col.ind = c("Col0","Col0","abc","abc"), 
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE, 
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)


## Calculamos la matrix de expresión media. 
col0 <- (normalized.gene.expression[,"col0_1"] + normalized.gene.expression[,"col0_2"])/2
abc <- (normalized.gene.expression[,"abc_1"] + normalized.gene.expression[,"abc_2"])/2

mean.expression <- matrix(c(col0,abc),ncol=2)
colnames(mean.expression) <- c("col0","abc")
rownames(mean.expression) <- rownames(normalized.gene.expression)
head(mean.expression)

## Previsualizamos el efecto de la mutación en un scatterplot.
plot(col0,abc,pch=19,cex=0.7,xlab="Col0",
     ylab=substitute(italic("abc")),
     cex.lab=1.25,
     col="grey")

##El paquete **limma** (Linear Models for Microarray Analysis) proporciona las 
##funciones necesarias para determinar los genes expresados de forma 
##diferencial (DEGs). 

library(limma)

## Especificamos el diseño experimental

experimental.design <- model.matrix(~ -1+factor(c(1,1,2,2)))
colnames(experimental.design) <- c("col0","abc")

##A continuación, ajustamos la estimación de los niveles de expresión de cada
##gen a un modelo lineal teniendo en cuenta el diseño experimental. Este paso
##fundamentalmente se calcula la media de las réplicas en cada condición.

linear.fit <- lmFit(normalized.gene.expression, experimental.design)

##Para especificar los constrastes a realizar utilizamos la función
##*makeContrasts* que recibe como entrada los contrastes a realizar separados 
##por comas y especificados con los nombres de las dos condiciones 
##correspondientes separadas por un guión -. También recibe el argumento 
##levels, un vector con el nombre de las condiciones:

contrast.matrix <- makeContrasts(abc-col0,levels=c("col0","abc"))

##Calculamos el fold-change y los p-valores correspondientes para cada gen en
##cada uno de los constrastes especificados utilizando las funciones *constrasts.fit* 
##y *eBayes*.

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

nrow(normalized.gene.expression)

col0.abc <- topTable(contrast.results, number=7507,coef=1,sort.by="logFC")
head(col0.abc)

log.fold.change <- col0.abc$logFC
q.value <- col0.abc$adj.P.Val
genes.ids <- rownames(col0.abc)
names(log.fold.change) <- genes.ids
names(q.value) <- genes.ids

activated.genes <- genes.ids[log.fold.change > 2 & q.value < 0.05]
repressed.genes <- genes.ids[log.fold.change < - 2 & q.value < 0.05]

length(activated.genes)
length(repressed.genes)

## Volcano plot
log.q.val <- -log10(q.value)
plot(log.fold.change,log.q.val,pch=19,col="grey",cex=0.8,
     xlim=c(-12,12),ylim = c(0,4), 
     xlab="log2(Fold-chage)",ylab="-log10(q-value)",cex.lab=1.5)

points(x = log.fold.change[activated.genes],
       y = log.q.val[activated.genes],col="red",cex=0.8,pch=19)
points(x = log.fold.change[repressed.genes],
       y = log.q.val[repressed.genes],col="blue",cex=0.8,pch=19)


## Enriquecimiento funcional. 
library(clusterProfiler)
library(org.At.tair.db)

activated.atha.enrich.go <- enrichGO(gene          = activated.genes,
                                     OrgDb         = org.At.tair.db,
                                     ont           = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable      = FALSE,
                                     keyType = "TAIR")


BiocManager::install("GOSemSim", force=T)
library(GOSemSim)


term_sim <- pairwise_termsim(activated.atha.enrich.go)

barplot(activated.atha.enrich.go,showCategory = 20)
dotplot(activated.atha.enrich.go,showCategory = 20)
emapplot(activated.atha.enrich.go,showCategory = 20)
cnetplot(activated.atha.enrich.go,showCategory = 20)

BiocManager::install("org.At.tair.db", force=T)
library(org.At.tair.db)


repressed.atha.enrich.go <- enrichGO(gene          = repressed.genes,
                                     OrgDb         = org.At.tair.db,
                                     ont           = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable      = FALSE,
                                     keyType = "TAIR")

barplot(repressed.atha.enrich.go,showCategory = 20)
dotplot(repressed.atha.enrich.go,showCategory = 20)
emapplot(repressed.atha.enrich.go,showCategory = 20)
cnetplot(repressed.atha.enrich.go,showCategory = 20)


activated.atha.enrich.kegg <- enrichKEGG(gene  = activated.genes,
                                         organism = "ath",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05)
df.activated.atha.enrich.kegg <- as.data.frame(activated.atha.enrich.kegg)
head(df.activated.atha.enrich.kegg)


repressed.atha.enrich.kegg <- enrichKEGG(gene  = repressed.genes,
                                         organism = "ath",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05)

df.repressed.atha.enrich.kegg <- as.data.frame(repressed.atha.enrich.kegg)
head(df.repressed.atha.enrich.kegg)

## Podemos visualizar las rutas KEGG usando el paquete pathveiw
library(pathview)
pathview(gene.data = sort(log.fold.change,decreasing = TRUE),
         pathway.id = "ath00941",
         species = "ath",
         limit = list(gene=max(abs(log.fold.change)), cpd=1),gene.idtype = "TAIR")




## Código para desarrollar una función gráfico de barras

gene <- activated.genes[1]

expression.matrix <- 2^normalized.gene.expression - 1
cond.names <- c("col0","abc")

gene.expression.barplot <- function(gene, expression.matrix,cond.names)
{
  expr.1 <- unlist(c(expression.matrix[gene, 1:2]))
  expr.2 <- unlist(c(expression.matrix[gene, 3:4]))
  
  mean.1 <- mean(expr.1)
  mean.2 <- mean(expr.2)
  
  sd.1 <- sd(expr.1)
  sd.2 <- sd(expr.2)
  
  means <- c(mean.1, mean.2)
  sds <- c(sd.1, sd.2)
  
  arrow.top <- means + sds
  arrow.bottom <- means - sds
  
  
  xpos <- barplot(means,ylim=c(0,1.5*max(arrow.top)),col=rainbow(2),
                  main=gene,names.arg = cond.names,
                  ylab="FPKM")
  arrows(xpos, arrow.top, xpos, arrow.bottom,code = 3,angle=90,length=0.05)
}

## Análisis de expresión génica diferencial con DESeq2
library(DESeq2)

pheno.data

gene.count.matrix <- read.table(file = "gene_count_matrix.csv",header = T,sep = ",")
head(gene.count.matrix)

sapply(X = strsplit(x = gene.count.matrix$gene_id,split = "\\|"),FUN = function(x){return(x[1])})

gene.ids <- sapply(X = strsplit(x = gene.count.matrix$gene_id,split = "\\|"),FUN = function(x){return(x[1])})

gene.count.matrix <- gene.count.matrix[,-1]
rownames(gene.count.matrix) <- gene.ids
head(gene.count.matrix)

dds <- DESeqDataSetFromMatrix(countData=gene.count.matrix, colData=pheno.data, design = ~ genotype)

dds$genotype <- relevel(dds$genotype, ref = "col0")

dds <- DESeq(dds) 
res <- results(dds) 
res

log.fold.change <- res$log2FoldChange
q.value <- res$padj
names(log.fold.change) <- genes.ids
names(q.value) <- genes.ids

activated.genes.deseq2 <- genes.ids[log.fold.change > 2 & q.value < 0.05]
activated.genes.deseq2 <- activated.genes.deseq2[!is.na(activated.genes.deseq2)]

repressed.genes.deseq2 <- genes.ids[log.fold.change < - 2 & q.value < 0.05]
repressed.genes.deseq2 <- repressed.genes.deseq2[!is.na(repressed.genes.deseq2)]

length(activated.genes.deseq2)
length(repressed.genes.deseq2)

activated.atha.enrich.go <- enrichGO(gene          = activated.genes.deseq2,
                                     OrgDb         = org.At.tair.db,
                                     ont           = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable      = FALSE,
                                     keyType = "TAIR")

barplot(activated.atha.enrich.go,showCategory = 20)
dotplot(activated.atha.enrich.go,showCategory = 20)
emapplot(activated.atha.enrich.go,showCategory = 20)
cnetplot(activated.atha.enrich.go,showCategory = 20)


