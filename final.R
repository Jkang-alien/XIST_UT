setwd('/home/jun/XIST_UT/')
load('methy.RData')

data[1:5,1:5]
index <- duplicated(rownames(data)) == FALSE
data_unique = data[index,]
data_var <- data_unique[,colSums(is.na(data_unique)) == 0]
methy = data.frame(ID = rownames(data_var), data_var)

rm(list=setdiff(ls(), "methy"))


library(cgdsr)

mycgds = CGDS("http://www.cbioportal.org/public-portal/")
cancerstudy = getCancerStudies(mycgds)
mycancerstudy = getCancerStudies(mycgds)[146,1]
caselist = getCaseLists(mycgds,mycancerstudy)
mycaselist = getCaseLists(mycgds,mycancerstudy)[7,1]
geneticprofile = getGeneticProfiles(mycgds,mycancerstudy)
mut_profile = getGeneticProfiles(mycgds,mycancerstudy)[8,1]
methyl_profile = getGeneticProfiles(mycgds,mycancerstudy)[7,1]
mut_BRCA1 = getProfileData(mycgds,'BRCA1',mut_profile,mycaselist)
methyl_BRCA1 = getProfileData(mycgds,'BRCA1',methyl_profile,mycaselist)
sum(methyl_BRCA1$BRCA1>0.4, na.rm = TRUE)
sum(methyl_BRCA1$BRCA1>0.4, na.rm = TRUE)

hist(methyl_BRCA1$BRCA1)
is.na(mut_BRCA1) -> na_mut_BRCA1
factor(na_mut_BRCA1, levels = c(FALSE,TRUE), labels = c('mutant', 'wild type')) -> BRCA1_m
BRCA1 = data.frame(ID = gsub('-01', '', gsub('\\.', '-', rownames(mut_BRCA1))),
                   BRCA1 = BRCA1_m,
                   BRCA1_meth = methyl_BRCA1$BRCA1)

rna_profile = getGeneticProfiles(mycgds,mycancerstudy)[3,1]
rna_XIST = getProfileData(mycgds,'XIST',rna_profile,mycaselist)
XIST = data.frame(ID = gsub('-01', '', gsub('\\.', '-', rownames(rna_XIST))),
                   XIST = rna_XIST)

myclinicaldata <- getClinicalData(mycgds,mycaselist)

set_4_subtype = getCancerStudies(mycgds)[89,1]
caselist4subtype = getCaseLists(mycgds,set_4_subtype)[2,1]
clinicaldata4subtype <- getClinicalData(mycgds,caselist4subtype)
subtype = data.frame(ID =  gsub('-01', '', gsub('\\.', '-', rownames(clinicaldata4subtype))),
                     subtype = clinicaldata4subtype$SUBTYPE)

myclinicaldata$ID <- gsub('-01', '', gsub('\\.', '-', rownames(myclinicaldata)))


myclinicaldata <- as.data.frame(unclass(myclinicaldata))
myclinicaldata$CLINICAL_STAGE <- as.character(myclinicaldata$CLINICAL_STAGE)
myclinicaldata$CLINICAL_STAGE[myclinicaldata$CLINICAL_STAGE == ''] <- NA
myclinicaldata$CLINICAL_STAGE[myclinicaldata$CLINICAL_STAGE %in% 
                                c('Stage I', 'Stage IA',
                                  'Stage IB', 'Stage IC',
                                  'Stage II', 'Stage IIA',
                                  'Stage IIB')] <- 'Stage I-II'

myclinicaldata$CLINICAL_STAGE[myclinicaldata$CLINICAL_STAGE %in% 
                                c('Stage III', 'Stage IIIA',
                                  'Stage IIIB', 'Stage IIIC', 
                                  'Stage IIIC1',
                                  'Stage IIIC2')] <- 'Stage III'

myclinicaldata$CLINICAL_STAGE[myclinicaldata$CLINICAL_STAGE %in% 
                                c('Stage IV', 'Stage IVA',
                                  'Stage IVB', 'Stage IVC')] <- 'Stage IV'
                              
myclinicaldata$CLINICAL_STAGE <- factor(myclinicaldata$CLINICAL_STAGE)

myclinicaldata$ETHNICITY[myclinicaldata$ETHNICITY == ''] <- NA
myclinicaldata$ETHNICITY[myclinicaldata$ETHNICITY == '[Not Evaluated]'] <- NA
myclinicaldata$RACE[myclinicaldata$RACE == ''] <- NA
myclinicaldata$RACE[myclinicaldata$RACE == 'AMERICAN INDIAN OR ALASKA NATIVE'] <- 'ASIAN'
myclinicaldata$RACE[myclinicaldata$RACE == 'NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER'] <- 'ASIAN'
myclinicaldata$RACE <- factor(myclinicaldata$RACE)
myclinicaldata$MENOPAUSE_STATUS[myclinicaldata$MENOPAUSE_STATUS == 'Indeterminate (neither Pre or Postmenopausal)'] <- NA
myclinicaldata$MENOPAUSE_STATUS[myclinicaldata$MENOPAUSE_STATUS == ''] <- NA
myclinicaldata$RESIDUAL_TUMOR[myclinicaldata$RESIDUAL_TUMOR == ''] <- NA
myclinicaldata$RESIDUAL_TUMOR <- factor(myclinicaldata$RESIDUAL_TUMOR)

summary(myclinicaldata$RACE)


########################################################################


myclinicaldata <- merge(myclinicaldata, subtype, 
                        by = 'ID',
                        all.x = TRUE)


########################################################################

data_clinical <- merge(myclinicaldata, methy, by = 'ID', all.y = TRUE)
data_clinical <- merge(BRCA1, data_clinical, by = 'ID', all.y = TRUE)
data_clinical <- merge(XIST, data_clinical, by = 'ID', all.y = TRUE)

data_clinical <- subset(data_clinical, 
                        ICD_O_3_HISTOLOGY %in% 
                          c('8380/3',
                            '8441/3'))

data_clinical$histological_type <- 
  factor(data_clinical$ICD_O_3_HISTOLOGY, 
         labels = c('Endometrioid adenocarcinoma',
                    'Serous cystadenocarcinoma, NOS'))

names(data_clinical) <- tolower(names(data_clinical)) 

dim(data_clinical)

##############################subset###############################
data_serous <- subset(data_clinical, 
                      histological_type == 'Serous cystadenocarcinoma, NOS')
data_endo <- subset(data_clinical, 
                      histological_type == 'Endometrioid adenocarcinoma')

grep('cg',colnames(data_serous))
data_hc_s = as.matrix(data_serous[,grep('cg',colnames(data_serous))])
rownames(data_hc_s) <- data_serous$ID

library(IlluminaHumanMethylation450k.db)

x <- IlluminaHumanMethylation450kCHRLOCEND

colv <- order(abs(unlist(lapply(as.list(x[colnames(data_hc_s)]), max))))


x <- IlluminaHumanMethylation450kMAP
ann_col = c()
ann_col[grep('p', as.list(x[colnames(data_hc_s)]))] <- 'p arm'
ann_col[grep('q', as.list(x[colnames(data_hc_s)]))] <- 'q arm'
ann_col = data.frame(Probe_location = ann_col)


library(NMF)
hc = hclust(dist(data_hc_s), 'complete')


memb_3 <- factor(cutree(hc, k = 3), levels = 1:3,
                 labels = c('A', 'B', 'C'))

data_serous$group <- memb_3

ann_s <- data.frame(Subtype_h = data_serous$histological_type,
                    cluster = data_serous$group,                  
                    Subtype_m = data_serous$SUBTYPE,
                    XIST = log(data_serous$XIST),
                    BRCA1 = data_serous$BRCA1,
                    BRCA1_M = factor(data_serous$BRCA1_meth > 0.4))

ann_colors_s <- list(Subtype_h = c('yellow', 'blue'),
                     cluster = c('darkolivegreen','deepskyblue', 'darkslategray1'),
                     Subtype_m = c('lightblue', 'orange', 'purple', 'green'),
                     BRCA1 = c('greenyellow', 'lightskyblue'),
                     BRCA1_M = c('seagreen', 'mediumpurple')
)

library(NMF)
library(Cairo)
CairoPDF(file = '/home/jun/XIST_UT/Figures/Figure1.pdf',
    width =7.5, height = 7.5, pointsize = 16)
layout(matrix(c(1,1,2,2), ncol = 2, byrow = TRUE),
       widths = c(1,1),
       heights = c(276,95)) 
#########################################################################

ann_e <- data.frame(Subtype_h = data_endo$histological_type,
                    cluster = data_endo$group,                  
                    Subtype_m = data_endo$subtype,
                    XIST = log(data_endo$xist))

ann_colors_e <- list(Subtype_h = c('yellow', 'blue'),
                     cluster = c('chartreuse', 'gold', 'cyan'),
                     Subtype_m = c('lightblue', 'orange', 'purple', 'green'))
ah<- aheatmap(data_hc_e, 
              hclustfun=function(d) hclust(d, method="complete"),
              Colv = colv,
              annRow = ann_e,
              #annCol = ann_col,
              annColors = ann_colors_e,
              #cex = 2,
              labRow = rep('',dim(data_hc_e)[1]),
              labCol = rep('',dim(data_hc_e)[2])
              #fontsize = 12,
              #cexCol = ,
              #fontsize = 16,
              #labCol = rep('',dim(data_hc_e)[2]),
              #Colv = colv
              #reorderfun = function(d, w) reorder(d, 10)
)

#######################################################################
ann_s <- data.frame(Subtype_h = data_serous$histological_type,
                    cluster = data_serous$group,                  
                    Subtype_m = data_serous$subtype,
                    XIST = log(data_serous$xist))

ann_colors_s <- list(Subtype_h = c('yellow', 'blue'),
                     cluster = c('darkolivegreen','deepskyblue', 'darkorchid'),
                     Subtype_m = c('lightblue', 'orange', 'purple', 'green'))

ah<- aheatmap(data_hc_s, 
              hclustfun=function(d) hclust(d, method="complete"),
              Colv = colv,
              annRow = ann_s,
              annCol = ann_col,
              annColors = ann_colors_s,
              #cex = 2,
              labRow = rep('',dim(data_hc_s)[1]),
              labCol = rep('',dim(data_hc_s)[2])
              #fontsize = 12,
              #cexCol = ,
              #fontsize = 16,
              #labCol = rep('',dim(data_hc_s)[2]),
              #Colv = colv
              #reorderfun = function(d, w) reorder(d, 10)
)

dev.off()
###########################################################################

library(sm)
CairoPDF(file = '/home/jun/XIST_UT/Figures/density.pdf',
    width =7.5, height =5, pointsize = 16)
par(mfrow = c(1,2))

sm.density.compare(log(data_endo[is.na(data_endo$xist) == FALSE,]$xist), 
                   data_endo[is.na(data_endo$xist) == FALSE,]$group, 
                   xlab="XIST level (log RSEM)",
                   col = ann_colors_e$cluster,
                   h = 0.5,
                   lty = rep(1,3),
                   ylab = '',
                   xlim = c(0,12),
                   ylim = c(0, 0.6),
                   ylab = 'Density'
)

legend('topright', levels(data_endo$group), 
       lty = rep(1,3),
       bty = 'n',
       col = ann_colors_e$cluster,
)

##############################################################################

sm.density.compare(log(data_serous[is.na(data_serous$xist) == FALSE,]$xist), 
                   data_serous[is.na(data_serous$xist) == FALSE,]$group,
                   #xlab="XIST level (log RSEM)",
                   h = 0.5,
                   lty = rep(1,3),
                   col = ann_colors_s$cluster,
                   xlim = c(0,12),
                   ylim = c(0, 0.6)
                  )

legend('topright', levels(data_serous$group), 
       lty = rep(1,3),
       bty = 'n',
       col = ann_colors_s$cluster
)


dev.off()

#################serous reduced##########################

pdf(file = '/home/jun/XIST_UT/Figures/hm_s_2_reduced.pdf',
    width =7.5, height = 7.5, pointsize = 16)
ah<- aheatmap(data_hc_s[,1:100], 
              hclustfun=function(d) hclust(d, method="complete"),
              Colv = 1:100,
              annRow = ann_s,
              annCol = ann_col[1:100,],
              annColors = ann_colors_s,
              #cex = 2,
              labRow = rep('',dim(data_hc_s)[1]),
              labCol = rep('',100)
              #fontsize = 12,
              #cexCol = ,
              #fontsize = 16,
              #labCol = rep('',dim(data_hc_e)[2]),
              #Colv = colv
              #reorderfun = function(d, w) reorder(d, 10)
)

dev.off()
ah$rowInd

IGV_S <- data.frame(ID = data_serous$ID[ah$rowInd],
                    order = 1:length(data_serous$ID))
write.table(IGV_S,
            file = 'IGV_S.txt',
            row.names=F,
            col.names=T,sep="\t", quote=FALSE)

library(sm)

pdf(file = '/home/jun/XIST_UT/Figures/XIST_s.pdf',
    width =3, height = 3, pointsize = 10)
sm.density.compare(log(data_serous[is.na(data_serous$XIST) == FALSE,]$XIST), 
                   data_serous[is.na(data_serous$XIST) == FALSE,]$group,
                   xlab="XIST level (log RSEM)",
                   h = 0.5,
                   lty = rep(1,3),
                   col = ann_colors_s$cluster,
                   ylab = 'Density')

legend('topright', levels(data_serous$group), 
       lty = rep(1,3),
       bty = 'n',
       col = ann_colors_s$cluster
)
dev.off()

###############  ENDOMETRIOID  ###############################
grep('cg',colnames(data_endo))
data_hc_e = as.matrix(data_endo[,grep('cg',colnames(data_endo))])
rownames(data_hc_e) <- data_endo$ID



library(IlluminaHumanMethylation450k.db)

x <- IlluminaHumanMethylation450kCHRLOCEND

colv <- order(abs(unlist(lapply(as.list(x[colnames(data_hc_e)]), max))))

ann_col = c()
ann_col[grep('p', as.list(x[colnames(data_hc_e)]))] <- 'p arm'
ann_col[grep('q', as.list(x[colnames(data_hc_e)]))] <- 'q arm'
ann_col = data.frame(Probe_location = ann_col)


library(NMF)
hc = hclust(dist(data_hc_e), 'complete')


memb_3 <- factor(cutree(hc, k = 3), levels = 1:3,
                 labels = c('D', 'E', 'F'))

data_endo$group <- memb_3

library(sm)
sm.density.compare(log(data_endo[is.na(data_endo$xist) == FALSE,]$xist), 
                   data_endo[is.na(data_endo$xist) == FALSE,]$group, 
                   xlab="XIST level (log RSEM)",
                   h = 0.5,
                   lty = rep(1,3),
                   ylab = ''
)

legend('topright', levels(data_endo$group), 
       lty = rep(1,3),
       bty = 'n',
       col = 2:5
)

ann_e <- data.frame(Subtype_h = data_endo$histological_type,
                    cluster = data_endo$group,                  
                    Subtype_m = data_endo$subtype,
                    XIST = log(data_endo$xist))

ann_colors_e <- list(Subtype_h = c('yellow', 'blue'),
                     cluster = c('chartreuse', 'gold', 'cyan'),
                     Subtype_m = c('lightblue', 'orange', 'purple', 'green'))

library(NMF)
par(mfrow = c(1,1))
#par(new = TRUE)
pdf(file = '/home/jun/XIST_UT/Figures/hm_e_2.pdf',
    width =7.5, height = 7.5, pointsize = 16)
ah<- aheatmap(data_hc_e, 
              hclustfun=function(d) hclust(d, method="complete"),
              Colv = colv,
              annRow = ann_e,
              #annCol = ann_col,
              annColors = ann_colors_e,
              #cex = 2,
              labRow = rep('',dim(data_hc_e)[1]),
              labCol = rep('',dim(data_hc_e)[2])
              #fontsize = 12,
              #cexCol = ,
              #fontsize = 16,
              #labCol = rep('',dim(data_hc_e)[2]),
              #Colv = colv
              #reorderfun = function(d, w) reorder(d, 10)
)

dev.off()


pdf(file = '/home/jun/XIST_UT/Figures/hm_e_2_reduced.pdf',
    width =7.5, height = 7.5, pointsize = 16)
ah<- aheatmap(data_hc_e[,1:100], 
              hclustfun=function(d) hclust(d, method="complete"),
              #Colv = colv[1:100],
              annRow = ann_e,
              #annCol = ann_col,
              annColors = ann_colors_e,
              #cex = 2,
              labRow = rep('',dim(data_hc_e)[1]),
              labCol = rep('',100)
              #fontsize = 12,
              #cexCol = ,
              #fontsize = 16,
              #labCol = rep('',dim(data_hc_e)[2]),
              #Colv = colv
              #reorderfun = function(d, w) reorder(d, 10)
)

dev.off()


##################################################################
IGV_E <- data.frame(ID = data_endo$ID[ah$rowInd],
                    order = 1:length(data_endo$id))
write.table(IGV_E,
            file = 'IGV_E.txt',
            row.names=F,
            col.names=T,sep="\t", quote=FALSE)


###################################################################
data_table_group_endo <- data_endo[,-grep('cg', colnames(data_endo))]
data_table_group_serous <- data_serous[,-grep('cg', colnames(data_serous))]
data_table_group <- rbind(data_table_group_endo, data_table_group_serous)
summary(data_table_group$group)
data_table_group$group_x <- data_table_group$group
data_table_group$group_x [data_table_group$id %in% c('TCGA-B5-A5OD',
                                               'TCGA-B5-A1N2',
                                               'TCGA-FI-A2EW',
                                               'TCGA-K6-A3WQ',
                                               'TCGA-EY-A1GS',
                                               'TCGA-AX-A3GI')] <- 'A'
summary(data_table_group$group_x)



library(reshape)
factor(combine_factor(data_table_group$group_x, c(2,1,3,1,2,3)), 
       labels = c('Preserved Xi',
                  'Partial reactivation of Xi',
                  'Two copies of Xa'))-> data_table_group$group_2 

summary(data_table_group$brca1_meth)

data_table_group$logxist <- log(data_table_group$xist)
data_table_group$methy_brca1 <- factor(data_table_group$brca1_meth > 0.4,
                                      levels = c(TRUE, FALSE),
                                      labels = c('methylated', 'unmethylated'))

###################Survival#########################
library(survival)
library(rms)



diff = survdiff(Surv(os_months, os_status == 'DECEASED')~ group_2, 
                data = data_table_group)
diff

cox <- coxph(Surv(os_months, os_status == 'DECEASED')~ 
              clinical_stage + group_2, 
             data = data_table_group)
s_cox <- summary(cox)
s_cox
paste(paste('(',
            paste(round(s_cox$conf.int, 2)[,3], round(s_cox$conf.int, 2)[,4], sep = '-'),
            sep = ''), ')', sep = '') -> CI
HR <- paste(round(s_cox$coefficients[,2],2), CI, sep = ' ')
p_value <- round(s_cox$coefficients[,5],3)

library(xtable)

# Bind columns together, and select desired rows
res <- cbind(HR, p_value)

# Print results in a LaTeX-ready form
write(print(xtable(res), type = 'html'), file = 'cox_os.html')

setwd('/home/jun/XIST_UT/Figures/')
svg(file = "Figure3.svg", pointsize = 10,
    width = 7.5 , height = 4,)
layout(matrix(c(1,2), ncol = 2, byrow = TRUE))
par(mar=c(5,3,1,4), mgp = c(2, 1, 0))

fit = npsurv(Surv(os_months, os_status == 'DECEASED')~ group_2, 
             data = data_table_group)

strata = levels(data_table_group$group_2)

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:3),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 1.0, strata, lty = c(1:3), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*10, 0.8, 'P-value: 0.052', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

fit = npsurv(Surv(dfs_months, dfs_status == 'Recurred/Progressed')~ group_2, 
             data = data_table_group)
fit

strata = levels(data_table_group$group_2)
par(mar=c(5,3,1,4), mgp = c(2, 1, 0))
survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         ylab = '',
         lty = c(1:3),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)

legend(40, 1.0, strata, lty = c(1:3), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*10, 0.8, 'P-value: 0.305', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

dev.off()

par(mfrow = c(1,1))
xtable(cox)
################# DFS #####################################
fit = npsurv(Surv(dfs_months, dfs_status == 'Recurred/Progressed')~ group_2, 
             data = data_table_group)
fit

strata = levels(data_table_group$group_2)

setwd('/home/jun/XIST_UT/Figures/')
svg(file = "DFS.svg", pointsize = 14)
par(mar=c(6,4,2,8), mgp = c(2, 1, 0))
survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:3),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 1.0, strata, lty = c(1:3), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*10, 0.8, 'P-value: 0.031', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

dev.off()




diff = survdiff(Surv(dfs_months, dfs_status == 'Recurred/Progressed')~ group_2, 
                data = data_table_group)
diff

fit

cox <- coxph(Surv(dfs_months, dfs_status == 'Recurred/Progressed')~ 
               clinical_stage +group_2, 
             data = data_table_group)
s_cox <- summary(cox)
s_cox
paste(paste('(',
            paste(round(s_cox$conf.int, 2)[,3], round(s_cox$conf.int, 2)[,4], sep = '-'),
            sep = ''), ')', sep = '') -> CI
HR <- paste(round(s_cox$coefficients[,2],2), CI, sep = ' ')
p_value <- round(s_cox$coefficients[,5],3)

library(xtable)

# Bind columns together, and select desired rows
res <- cbind(HR, p_value)

# Print results in a LaTeX-ready form
write(print(xtable(res), type = 'html'), file = 'cox_dfs.html')


######################## Survival subtype #####################


####################### Endometrioid ##########################


diff = survdiff(Surv(os_months, os_status == 'DECEASED')~ group_2, 
                data = subset(data_table_group, histological_type == 'Endometrioid adenocarcinoma'))
diff

cox <- coxph(Surv(os_months, os_status == 'DECEASED')~ 
               clinical_stage + group_2, 
             data = subset(data_table_group, histological_type == 'Endometrioid adenocarcinoma'))
s_cox <- summary(cox)
s_cox
paste(paste('(',
            paste(round(s_cox$conf.int, 2)[,3], round(s_cox$conf.int, 2)[,4], sep = '-'),
            sep = ''), ')', sep = '') -> CI
HR <- paste(round(s_cox$coefficients[,2],2), CI, sep = ' ')
p_value <- round(s_cox$coefficients[,5],3)

library(xtable)

# Bind columns together, and select desired rows
res <- cbind(HR, p_value)

# Print results in a LaTeX-ready form
write(print(xtable(res), type = 'html'), file = './Figures/cox_os_e.html')

library(Cairo)
setwd('/home/jun/XIST_UT/Figures/')
CairoSVG(file = "survival_subtype.svg", pointsize = 10,
    width = 7.5 , height = 8,)
layout(matrix(c(1:4), ncol = 2, byrow = TRUE))
par(mar=c(6,3,2,4), mgp = c(2, 1, 0))

fit = npsurv(Surv(os_months, os_status == 'DECEASED')~ group_2, 
             data = subset(data_table_group, histological_type == 'Endometrioid adenocarcinoma'))

strata = levels(subset(data_table_group, histological_type == 'Endometrioid adenocarcinoma')$group_2)

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:3),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 1.0, strata, lty = c(1:3), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*10, 0.8, 'P-value: 0.052', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

fit = npsurv(Surv(dfs_months, dfs_status == 'Recurred/Progressed')~ group_2, 
             data = subset(data_table_group, histological_type == 'Endometrioid adenocarcinoma'))
fit

strata = levels(subset(data_table_group, histological_type == 'Endometrioid adenocarcinoma')$group_2)
par(mar=c(6,3,2,4), mgp = c(2, 1, 0))
survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         ylab = '',
         lty = c(1:3),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)

legend(40, 1.0, strata, lty = c(1:3), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*10, 0.8, 'P-value: 0.305', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

fit = npsurv(Surv(os_months, os_status == 'DECEASED')~ group_2, 
             data = subset(data_table_group, histological_type == 'Serous cystadenocarcinoma, NOS'))

strata = levels(subset(data_table_group, histological_type == 'Serous cystadenocarcinoma, NOS')$group_2)

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:3),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 1.0, strata, lty = c(1:3), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*10, 0.8, 'P-value: 0.755', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

fit = npsurv(Surv(dfs_months, dfs_status == 'Recurred/Progressed')~ group_2, 
             data = subset(data_table_group, histological_type == 'Serous cystadenocarcinoma, NOS'))
fit

strata = levels(subset(data_table_group, histological_type == 'Serous cystadenocarcinoma, NOS')$group_2)
par(mar=c(6,3,2,4), mgp = c(2, 1, 0))
survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         ylab = '',
         lty = c(1:3),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)

legend(40, 1.0, strata, lty = c(1:3), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*10, 0.8, 'P-value: 0.077', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')


dev.off()

diff = survdiff(Surv(dfs_months, dfs_status == 'Recurred/Progressed')~ group_2, 
                data = subset(data_table_group, histological_type == 'Endometrioid adenocarcinoma'))
diff

fit

cox <- coxph(Surv(dfs_months, dfs_status == 'Recurred/Progressed')~ 
               clinical_stage +group_2, 
             data = subset(data_table_group, histological_type == 'Endometrioid adenocarcinoma'))
s_cox <- summary(cox)
s_cox
paste(paste('(',
            paste(round(s_cox$conf.int, 2)[,3], round(s_cox$conf.int, 2)[,4], sep = '-'),
            sep = ''), ')', sep = '') -> CI
HR <- paste(round(s_cox$coefficients[,2],2), CI, sep = ' ')
p_value <- round(s_cox$coefficients[,5],3)

library(xtable)

# Bind columns together, and select desired rows
res <- cbind(HR, p_value)

# Print results in a LaTeX-ready form
write(print(xtable(res), type = 'html'), file = 'cox_dfs_e.html')

########################### Serous ###########################

diff = survdiff(Surv(os_months, os_status == 'DECEASED')~ group_2, 
                data = subset(data_table_group, histological_type == 'Serous cystadenocarcinoma, NOS'))
diff

cox <- coxph(Surv(os_months, os_status == 'DECEASED')~ 
               clinical_stage + group_2, 
             data = subset(data_table_group, histological_type == 'Serous cystadenocarcinoma, NOS'))
s_cox <- summary(cox)
s_cox
paste(paste('(',
            paste(round(s_cox$conf.int, 2)[,3], round(s_cox$conf.int, 2)[,4], sep = '-'),
            sep = ''), ')', sep = '') -> CI
HR <- paste(round(s_cox$coefficients[,2],2), CI, sep = ' ')
p_value <- round(s_cox$coefficients[,5],3)

library(xtable)

# Bind columns together, and select desired rows
res <- cbind(HR, p_value)

# Print results in a LaTeX-ready form
write(print(xtable(res), type = 'html'), file = 'cox_os_s.html')

library(Cairo)
setwd('/home/jun/XIST_UT/Figures/')
CairoSVG(file = "oa_s.svg", pointsize = 10,
         width = 7.5 , height = 4,)
layout(matrix(c(1,2), ncol = 2, byrow = TRUE))
par(mar=c(5,3,1,4), mgp = c(2, 1, 0))

fit = npsurv(Surv(os_months, os_status == 'DECEASED')~ group_2, 
             data = subset(data_table_group, histological_type == 'Serous cystadenocarcinoma, NOS'))

strata = levels(subset(data_table_group, histological_type == 'Serous cystadenocarcinoma, NOS')$group_2)

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:3),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 1.0, strata, lty = c(1:3), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*10, 0.8, 'P-value: 0.755', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

fit = npsurv(Surv(dfs_months, dfs_status == 'Recurred/Progressed')~ group_2, 
             data = subset(data_table_group, histological_type == 'Serous cystadenocarcinoma, NOS'))
fit

strata = levels(subset(data_table_group, histological_type == 'Serous cystadenocarcinoma, NOS')$group_2)
par(mar=c(5,3,1,4), mgp = c(2, 1, 0))
survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         ylab = '',
         lty = c(1:3),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)

legend(40, 1.0, strata, lty = c(1:3), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*10, 0.8, 'P-value: 0.077', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

dev.off()

diff = survdiff(Surv(dfs_months, dfs_status == 'Recurred/Progressed')~ group_2, 
                data = subset(data_table_group, histological_type == 'Serous cystadenocarcinoma, NOS'))
diff

fit

cox <- coxph(Surv(dfs_months, dfs_status == 'Recurred/Progressed')~ 
               clinical_stage +group_2, 
             data = subset(data_table_group, histological_type == 'Serous cystadenocarcinoma, NOS'))
s_cox <- summary(cox)
s_cox
paste(paste('(',
            paste(round(s_cox$conf.int, 2)[,3], round(s_cox$conf.int, 2)[,4], sep = '-'),
            sep = ''), ')', sep = '') -> CI
HR <- paste(round(s_cox$coefficients[,2],2), CI, sep = ' ')
p_value <- round(s_cox$coefficients[,5],3)

library(xtable)

# Bind columns together, and select desired rows
res <- cbind(HR, p_value)

# Print results in a LaTeX-ready form
write(print(xtable(res), type = 'html'), file = 'cox_dfs_s.html')


        
######################## Table ################################

data_table_group$menopause_status[data_table_group$menopause_status == 'Indeterminate (neither Pre or Postmenopausal)'] <- NA

library(compareGroups)
compareGroups(group_2~ 
                age+
                ethnicity +
                race+
                menopause_status +
                residual_tumor +
                clinical_stage+
                histological_type+
                subtype+
                logxist +
                brca1 +
                methy_brca1,,
              data = data_table_group) -> table_a

createTable(table_a,
            hide = NA, show.all = TRUE,
            digits.ratio = 1,digits = 1) -> table_a
table_a

export2csv(table_a, '/home/jun/XIST_UT/Figures/table_a.csv')     
    

################ RNA #####################################

setwd('/home/jun/XIST_UT/')

rna <- read.delim('./gdac.broadinstitute.org_UCEC.Merge_rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2014120600.0.0/UCEC.rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt')
#row_number_XIST <- grep('XIST', rna$Hybridization.REF)

rna = t(rna)
data_rna <- (rna[-1, -1])
dim(data_rna)
matrix_rna <- matrix(as.numeric(data_rna), nrow = 381, byrow = FALSE )
dim(matrix_rna)
length(unlist(strsplit(rna[1,-1], '\\|'))[seq(1,2*(dim(rna)[2]-1), 2)])
colnames(matrix_rna) <- unlist(strsplit(rna[1,-1], '\\|'))[seq(1,2*(dim(rna)[2]-1), 2)]
matrix_rna_X <- matrix_rna[, colnames(matrix_rna) %in% my.regions$hgnc_symbol]
ID <- gsub('.[0-9A-Z]{3}.[0-9A-Z]{3}.[0-9A-Z]{4}.[0-9A-Z]{2}', '', rownames(data_rna))
index_dup <- duplicated(ID)
length(index_dup)

df_rna = data.frame(id = gsub('\\.', '-',ID), matrix_rna_X)
df_rna[1:5,100:105]
rm(matrix_rna, matrix_rna_X, data_rna)


library(dplyr)

data_table_group %>%
  select(id, group_2, histological_type) %>%
  filter((group_2 != 'Partial reactivation of Xi') &
           histological_type == 'Serous cystadenocarcinoma, NOS') %>%
  select(id, group_2) -> group_serous

data_rna_group_serous <- merge(group_serous, df_rna, by.x = 'id', by.y = 'ID')


p_value <- c()

for (i in 3:dim(data_rna_group_serous)[2]){
  a <- wilcox.test(data_rna_group_serous[,i]~data_rna_group_serous$group_2)
  b <- a$p.value
  #c <- a$estimate
  p_value <- rbind(p_value, b)
  #diff <- rbind(diff, c)
}


tapply(data_rna_group_serous[,3], data_rna_group_serous$group_2, quantile)$A[4]

diff <- c()
for (i in 3:dim(data_rna_group_serous)[2]){
  q <- tapply(data_rna_group_serous[,i], data_rna_group_serous$group_2, mean)
  d <- q[3] - q[1]
  diff <- rbind(diff, d)
}

loc <- matrix( c(diff[abs(diff) > 1000 & p_value < 0.1], 
                 -log(p_value, base =10)[abs(diff) > 1000 & p_value <0.1]), 
               nrow = 2,
               byrow = TRUE)

colnames(loc) <- colnames(matrix_rna_X)[abs(diff) > 1000 & 
                                          p_value < 0.1]

svg(file = '/home/jun/XIST_UT/Figures/volcano_mean.svg',
    width =7.5, height = 7.5, pointsize = 16)                
vp <- plot(diff, -log(p_value, base = 10),
           xlab = 'Difference of mRNA expression',
           ylab = '-Log10(p-value)')
text(loc[1,]+500, loc[2,], colnames(loc), adj = c(0,0.5),
     cex =0.5)
#arrows(0.3, 4, 0.4, 3.25, angle = 30, length = 0.05)
dev.off()

p_value[is.na(p_value)]<-1
colnames(matrix_rna_X)[p_value < 0.01]

colnames(matrix_rna_X)[abs(diff) > 1000 & 
  p_value < 0.05]




f_MAGE <- genes_X$hgnc_symbol [grep('MAGE', genes_X$hgnc_symbol) ]
f_GAGE <- genes_X$hgnc_symbol [grep('GAGE', genes_X$hgnc_symbol) ]
f_SAGE <- genes_X$hgnc_symbol [grep('SAGE', genes_X$hgnc_symbol) ]
f_SSX <- genes_X$hgnc_symbol [grep('SSX', genes_X$hgnc_symbol) ]
f_CTAG <- genes_X$hgnc_symbol [grep('CTAG', genes_X$hgnc_symbol) ]

XCTs <- c(f_CTAG, f_GAGE, f_MAGE, f_SAGE, f_SSX)

matrix_rna_XCTs <- matrix_rna[, colnames(matrix_rna) %in% XCTs]
log(apply(matrix_rna_XCTs, 1, mean))
XCT_family <- c(rep('CTAG', 2), rep('GAGAE', 8), rep('MAGEA', 11),
                rep('MAGEB', 8), rep('MAGEC', 3), rep('MAGED', 4),
                rep('MAGEE', 2), rep('MAGEH', 1),
                rep('SAGE', 1), rep('SSX', 8))


df_XCTs = data.frame(id = gsub('\\.', '-',ID), 
                     matrix_rna_XCTs)
data_table_group$subtype_2 <-factor(data_table_group$subtype, 
                                     levels(data_table_group$subtype)[c(3,4,2,1)])


data_XCTs <- merge(data_table_group, df_XCTs, by = 'id')



levels(data_XCTs$group_2)

Group <- data_XCTs$group_2


library(NMF)

CairoSVG(file = '/home/jun/XIST_UT/Figures/bp_R1.svg',
    width =9, height = 7.5, pointsize = 10)
layout(matrix(c(1,1,1,1,2:5), ncol = 4, byrow = TRUE),
       widths = c(1,1,1), heights = rep(1.5,1), respect = FALSE)
par (mar=c(4.1, 5.1, 4.1, 1.1))

aheatmap(log(1+data_XCTs[,49:96]),
         hclustfun=function(d) hclust(d, method="complete"),
         annRow = Group,
         #labRow = NA, #rep('',dim(data_XCTs)[1]),
         Rowv = order(data_XCTs$group_2),
         Colv = order(colnames(data_XCTs[,48:95])))

par (mar=c(5.1, 5, 1, 1.1))
bp <- boxplot(logXCTs~subtype_2, data_XCT_clinical,
              ylab = 'mean RNA expression of\nCTAs (log RSEM)',
              lab = '',
              #border = "white",
              frame = FALSE,
              names = rep('',4),
              xaxt='n',
              outpch = NA,
              ylim = c(3,8.5)
              
)

axis(mgp=c(3.5, 2, 0),
     side = 1,
     at = 1:4,
     pos = 2.5,
     tck = -0.01,
     labels=c('MSI',
              "POLE",
              'CN low',
              "CN high"
              ))

stripchart(logXCTs~subtype_2, data_XCT_clinical,
           vertical=T,pch=1, method="jitter",cex=1,add=T)

x <- c(1, 1, 3, 3)
y <- c(7.9, 8.0, 8.0, 7.9)-0.9
lines (x,y)
x <- c(2, 2, 4, 4)
y <- c(7.9, 8.0, 8.0, 7.9)-0.1
lines (x,y)

text (3, 8.2, '***')
text (2, 7.4, 'NS')

par (mar=c(5.1, 0, 1, 3.1))
bp <- boxplot(logXCTs~histological_type, data_XCT_clinical,
              xaxt = 'n',
              yaxt = 'n',
              #ylab = 'mean RNA expression of CTAs (log RSEM)',
              outpch = NA,
              #border = "white",
              frame = FALSE,
              ylim = c(3,8.5)
)

axis(mgp=c(3.5, 2, 0),
     side = 1,
     at = 1:2,
     pos = 2.5,
     tck = -0.01,
     labels=c('Enodmetrioid\nadenocarcinoma',
              "Serous\nadenocarcinoma"))

stripchart(logXCTs~histological_type, data_XCT_clinical,
           vertical=T,pch=1, method="jitter",cex=1,add=T)

x <- c(1, 1, 2, 2)
y <- c(7.9, 8.0, 8.0, 7.9)+0.1
lines (x,y)
text (1.5, 8.3, '**')

x <- c(1.5, 1.5, 2.5, 2.5, 1.5)
y <- c(3,8,8,3,3)-0.1
lines (x,y)


par (mar=c(5.1, 1, 1, 1.1))
bp <- boxplot(logXCTs ~ group_2, 
              subset(data_XCT_clinical,
                     histological_type == 'Endometrioid adenocarcinoma'),
              #ylab = 'mean RNA expression of CTAs (log RSEM)',
              frame = FALSE,
              xaxt = 'n',
              yaxt = 'n',
              outpch = NA,
              ylim = c(3,8.5)
)

axis(mgp = c(3.5,3,0),
     side = 1,
     at = 1:3,
     pos = 2.5,
     labels = c('Preserved Xi',
                'Partial\nreactivation\nof Xi',
                'Two Xa')
     
)

stripchart(logXCTs ~ group_2, 
           subset(data_XCT_clinical,
                  histological_type == 'Endometrioid adenocarcinoma'),
           vertical=T,pch=1, method="jitter",cex=1,add=T)


x <- c(2.02, 2.02, 3, 3)
y <- c(7.2,7.3,7.3,7.2)+0.8
lines (x,y)

x <- c(1, 1, 1.98, 1.98)
y <- c(7.2, 7.3, 7.3, 7.2)+0.8
lines (x,y)

x <- c(1, 1, 3, 3)
y <- c(7.8, 7.9, 7.9, 7.8)+0.8
lines (x,y)

text (1.5, 7.45+0.8, '***')
text (2.5, 7.45+0.8, 'NS')
text (2, 8+0.8, '***')

par (mar=c(5.1, 1, 1, 1.1))
bp <- boxplot(logXCTs ~ group_2, 
              subset(data_XCT_clinical,
                     histological_type == 'Serous cystadenocarcinoma, NOS'),
              #ylab = 'mean RNA expression of CTAs (log RSEM)',
              frame = FALSE,
              xaxt = 'n',
              yaxt = 'n',
              outpch = NA,
              ylim = c(3,8.5)
)

axis(mgp = c(3.5,3,0),
      side = 1,
      at = 1:3,
      pos = 2.5,
      labels = c('Preserved Xi',
                 'Xa+',
                 'Two Xa')
      
      )

stripchart(logXCTs ~ group_2, 
           subset(data_XCT_clinical,
                  histological_type == 'Serous cystadenocarcinoma, NOS'),
           vertical=T,pch=1, method="jitter",cex=1,add=T)


x <- c(2.02, 2.02, 3, 3)
y <- c(7.2,7.3,7.3,7.2)
lines (x,y)

x <- c(1, 1, 1.98, 1.98)
y <- c(7.2, 7.3, 7.3, 7.2)
lines (x,y)

x <- c(1, 1, 3, 3)
y <- c(7.8, 7.9, 7.9, 7.8)
lines (x,y)

text (1.5, 7.45, 'NS')
text (2.5, 7.45, 'NS')
text (2, 8, '*')


dev.off()


############################################################
layout(matrix(1))

boxplot(SSX7~group_x, data_rna_group_serous_2)
a1 <- aov(logXCTs ~ subtype, data_XCT_clinical)
posthoc <- TukeyHSD(x=a1, 'subtype', conf.level=0.95)
summary(a1)
posthoc

t.test(logXCTs ~ data_XCT_clinical$subtype == 'Copy-number high (Serous-like)', data_XCT_clinical)
t.test(logXCTs~histological_type, data_XCT_clinical)
a2 <- aov(logXCTs ~ group_2, 
    subset(data_XCT_clinical,
           histological_type == 'Serous cystadenocarcinoma, NOS'))
posthoc2 <- TukeyHSD(x=a2, 'group_2', conf.level=0.95)
summary(a2)
posthoc2

a3 <- aov(logXCTs ~ group_2, 
          subset(data_XCT_clinical,
                 histological_type == 'Endometrioid adenocarcinoma'))
posthoc3 <- TukeyHSD(x=a3, 'group_2', conf.level=0.95)
summary(a3)
posthoc3


######################################################################

f_BAGE <- c('BAGE','BAGE2','BAGE3','BAGE4','BAGE5')
f_CAGE <- my.regions$hgnc_symbol [grep('CAGE', my.regions$hgnc_symbol) ]
f_CTCFL <- my.regions$hgnc_symbol [grep('CTCFL', my.regions$hgnc_symbol) ]
f_BRDT <- my.regions$hgnc_symbol [grep('BRDT', my.regions$hgnc_symbol) ]
f_DDX43 <- my.regions$hgnc_symbol [grep('DDX43', my.regions$hgnc_symbol) ]
f_ACRBP <- my.regions$hgnc_symbol [grep('ACRBP', my.regions$hgnc_symbol) ]
f_SPO11 <- my.regions$hgnc_symbol [grep('SPO11', my.regions$hgnc_symbol) ]
f_SYCP1 <- my.regions$hgnc_symbol [grep('SYCP1', my.regions$hgnc_symbol) ]

n_XCTs <- c(f_BAGE, f_CAGE, f_CTCFL, f_BRDT, f_DDX43, f_ACRBP, f_SPO11, f_SYCP1)

matrix_rna_n_XCTs <- matrix_rna[, colnames(matrix_rna) %in% n_XCTs]
log(apply(matrix_rna_n_XCTs, 1, mean))
n_XCT_family <- c(rep('CTAG', 2), rep('GAGAE', 8), rep('MAGEA', 11),
                rep('MAGEB', 8), rep('MAGEC', 3), rep('MAGED', 4),
                rep('MAGEE', 2), rep('MAGEH', 1),
                rep('SAGE', 1), rep('SSX', 8))


df_n_XCTs = data.frame(id = gsub('\\.', '-',ID), 
                     matrix_rna_n_XCTs,
                     lognxct = log(apply(matrix_rna_n_XCTs, 1, mean)))

data_n_XCTs <- merge(data_table_group, df_n_XCTs, by = 'id')

boxplot(lognxct~group_2, data_n_XCTs)

levels(data_n_XCTs$group_2)

Group <- data_n_XCTs$group_2

library(NMF)


###################################################################
#####################################################################



CairoSVG(file = '/home/jun/XIST_UT/Figures/somatic_R1.svg',
    width =7.5, height = 8.7, pointsize = 10)
layout(matrix(c(1,1,1,1,2:9), ncol = 4, byrow = TRUE),
       widths = c(1,1,1), heights = rep(1.5,1,1), respect = FALSE)
par (mar=c(4.1, 5.1, 4.1, 1.1))

aheatmap(log(1+(data_n_XCTs[,49:56])),
         hclustfun=function(d) hclust(d, method="complete"),
         annRow = Group,
         #labRow = NA, #rep('',dim(data_XCTs)[1]),
         Rowv = order(data_n_XCTs$group_2),
         Colv = order(colnames(data_n_XCTs[,49:56])))

par (mar=c(5.1, 5, 1, 1.1))
bp <- boxplot(lognxct~subtype_2, data_n_XCTs,
              ylab = 'mean RNA expression of\nCTAs (log RSEM)',
              lab = '',
              #border = "white",
              frame = FALSE,
              names = rep('',4),
              xaxt='n',
              outpch = NA,
              ylim = c(2,8)
              
)

axis(mgp=c(3.5, 2, 0),
     side = 1,
     at = 1:4,
     pos = NA,
     tck = -0.01,
     labels=c('MSI',
              "POLE",
              "CN low",
              'CN high'))

stripchart(lognxct~subtype_2, data_n_XCTs,
           vertical=T,pch=1, method="jitter",cex=1,add=T)

x <- c(1, 1, 3, 3)
y <- c(7.9, 8.0, 8.0, 7.9)-0.9
lines (x,y)

x <- c(2, 2, 4, 4)
y <- c(7.9, 8.0, 8.0, 7.9)-0.1
lines (x,y)

text (3, 8.1, '***')
text (2, 7.4, 'NS')

par (mar=c(5.1, 0, 1, 3.1))
bp <- boxplot(lognxct~histological_type, data_n_XCTs,
              xaxt = 'n',
              yaxt = 'n',
              #ylab = 'mean RNA expression of CTAs (log RSEM)',
              outpch = NA,
              #border = "white",
              frame = FALSE,
              ylim = c(2,8)
)

axis(mgp=c(3.5, 2, 0),
     side = 1,
     at = 1:2,
     pos = NA,
     tck = -0.01,
     labels=c('Enodmetrioid\nadenocarcinoma',
              "Serous\nadenocarcinoma"))

stripchart(lognxct~histological_type, data_n_XCTs,
           vertical=T,pch=1, method="jitter",cex=1,add=T)

x <- c(1, 1, 2, 2)
y <- c(7.9, 8.0, 8.0, 7.9)-0.5
lines (x,y)
text (1.5, 7.7, '***')

x <- c(1.5, 1.5, 2.5, 2.5, 1.5)
y <- c(3,8,8,3,3)-1
lines (x,y)


par (mar=c(5.1, 1, 1, 1.1))
bp <- boxplot(lognxct ~ group_2, 
              subset(data_n_XCTs,
                     histological_type == 'Endometrioid adenocarcinoma'),
              #ylab = 'mean RNA expression of CTAs (log RSEM)',
              frame = FALSE,
              xaxt = 'n',
              yaxt = 'n',
              outpch = NA,
              ylim = c(2,8)
)

axis(mgp = c(3.5,3,0),
     side = 1,
     at = 1:3,
     pos = NA,
     labels = c('Preserved Xi',
                'Partial\nreactivation\nof Xi',
                'Two Xa')
     
)

stripchart(lognxct ~ group_2, 
           subset(data_n_XCTs,
                  histological_type == 'Endometrioid adenocarcinoma'),
           vertical=T,pch=1, method="jitter",cex=1,add=T)


x <- c(2.02, 2.02, 3, 3)
y <- c(7.2,7.3,7.3,7.2)
lines (x,y)

x <- c(1, 1, 1.98, 1.98)
y <- c(7.2, 7.3, 7.3, 7.2)
lines (x,y)

x <- c(1, 1, 3, 3)
y <- c(7.8, 7.9, 7.9, 7.8)
lines (x,y)

text (1.5, 7.45, 'NS')
text (2.5, 7.45, 'NS')
text (2, 8, '**')


par (mar=c(5.1, 1, 1, 1.1))
bp <- boxplot(lognxct ~ group_2, 
              subset(data_n_XCTs,
                     histological_type == 'Serous cystadenocarcinoma, NOS'),
              #ylab = 'mean RNA expression of CTAs (log RSEM)',
              frame = FALSE,
              xaxt = 'n',
              yaxt = 'n',
              outpch = NA,
              ylim = c(2,8)
)

axis(mgp = c(3.5,3,0),
     side = 1,
     at = 1:3,
     pos = NA,
     labels = c('Preserved Xi',
                'Partial\nreactivation\nof Xi',
                'Two Xa')
     
)

stripchart(lognxct ~ group_2, 
           subset(data_n_XCTs,
                  histological_type == 'Serous cystadenocarcinoma, NOS'),
           vertical=T,pch=1, method="jitter",cex=1,add=T)


x <- c(2, 2, 3, 3)
y <- c(7.9, 8.0, 8.0, 7.9)-0.9
lines (x,y)

x <- c(1, 1, 2.5, 2.5)
y <- c(7.9, 8.0, 8.0, 7.9)-0.1
lines (x,y)

text (1.75, 8.1, '*')
text (2.5, 7.4, 'NS')

par (mar=c(5.1, 5, 1, 1.1))
bp <- boxplot(methy_somatic~subtype_2, somatic,
              ylab = 'methylation of\nsomatic CpG island (beta value)',
              lab = '',
              #border = "white",
              frame = FALSE,
              names = rep('',4),
              xaxt='n',
              outpch = NA,
              ylim = c(0.25,0.5)
              
)

axis(mgp=c(3.5, 2, 0),
     side = 1,
     at = 1:4,
     pos = NA,
     tck = -0.01,
     labels=c('MSI',
              "POLE",
              "CN low",
              'CN high'))

stripchart(methy_somatic~subtype_2, somatic,
           vertical=T,pch=1, method="jitter",cex=1,add=T)

x <- c(1, 1, 3, 3)
y <- (c(7.9, 8.0, 8.0, 7.9)-0.9)/10-0.26
lines (x,y)
x <- c(2, 2, 4, 4)
y <- (c(7.9, 8.0, 8.0, 7.9)-0.9)/10-0.23
lines (x,y)

text (3, 0.49, '***')
text (2, 0.46, 'NS')

par (mar=c(5.1, 0, 1, 3.1))
bp <- boxplot(methy_somatic~histological_type, somatic,
              xaxt = 'n',
              yaxt = 'n',
              #ylab = 'mean RNA expression of CTAs (log RSEM)',
              outpch = NA,
              #border = "white",
              frame = FALSE,
              outpch = NA,
              ylim = c(0.25,0.5)
              
)

axis(mgp=c(3.5, 2, 0),
     side = 1,
     at = 1:2,
     pos = NA,
     tck = -0.01,
     labels=c('Enodmetrioid\nadenocarcinoma',
              "Serous\nadenocarcinoma"))

stripchart(methy_somatic~histological_type, somatic,
           vertical=T,pch=1, method="jitter",cex=1,add=T)

x <- c(1, 1, 2, 2)
y <- (c(7.9, 8.0, 8.0, 7.9)-0.9)/10-0.26
lines (x,y)
text (1.5, 0.46, '***')

x <- c(1.5, 1.5, 2.5, 2.5, 1.5)
y <- c(0.27, 0.42, 0.42, 0.27, 0.27)
lines (x,y)

par (mar=c(5.1, 1, 1, 1.1))
bp <- boxplot(methy_somatic ~ group_2, 
              subset(somatic,
                     histological_type == 'Endometrioid adenocarcinoma'),
              #ylab = 'mean RNA expression of CTAs (log RSEM)',
              frame = FALSE,
              xaxt = 'n',
              yaxt = 'n',
              outpch = NA,
              ylim = c(0.25,0.5)
)

axis(mgp = c(3.5,3,0),
     side = 1,
     at = 1:3,
     pos = NA,
     labels = c('Preserved Xi',
                'Partial\nreactivation\nof Xi',
                'Two Xa')
     
)

stripchart(methy_somatic ~ group_2, 
           subset(somatic,
                  histological_type == 'Endometrioid adenocarcinoma'),
           vertical=T,pch=1, method="jitter",cex=1,add=T)


x <- c(2.02, 2.02, 3, 3)
y <- c(7.2,7.3,7.3,7.2)
lines (x,y)

x <- c(1, 1, 1.98, 1.98)
y <- c(7.2, 7.3, 7.3, 7.2)
lines (x,y)

x <- c(1, 1, 3, 3)
y <- c(7.8, 7.9, 7.9, 7.8)
lines (x,y)

text (1.5, 7.45, 'NS')
text (2.5, 7.45, 'NS')
text (2, 8, '**')


x <- c(2.02, 2.02, 3, 3)
y <- (c(7.9, 8.0, 8.0, 7.9)-0.9)/10-0.26
lines (x,y)

x <- c(1, 1, 1.98, 1.98)
y <- (c(7.9, 8.0, 8.0, 7.9)-0.9)/10-0.26
lines (x,y)

x <- c(1, 1, 3, 3)
y <- (c(7.9, 8.0, 8.0, 7.9)-0.9)/10-0.22
lines (x,y)

text (1.5, 0.46, '***')
text (2.5, 0.46, 'NS')
text (2, 0.49, '***')


par (mar=c(5.1, 1, 1, 1.1))
bp <- boxplot(methy_somatic ~ group_2, 
              subset(somatic,
                     histological_type == 'Serous cystadenocarcinoma, NOS'),
              #ylab = 'mean RNA expression of CTAs (log RSEM)',
              frame = FALSE,
              xaxt = 'n',
              yaxt = 'n',
              outpch = NA,
              ylim = c(0.25,0.5)
)

axis(mgp = c(3.5,3,0),
     side = 1,
     at = 1:3,
     pos = NA,
     labels = c('Preserved Xi',
                'Partial\nreactivation\nof Xi',
                'Two Xa')
     
)

stripchart(methy_somatic ~ group_2, 
           subset(somatic,
                  histological_type == 'Serous cystadenocarcinoma, NOS'),
           vertical=T,pch=1, method="jitter",cex=1,add=T)


x <- c(2, 2, 3, 3)
y <- (c(7.9, 8.0, 8.0, 7.9)-0.9)/10-0.26
lines (x,y)

x <- c(1, 1, 2.5, 2.5)
y <- (c(7.9, 8.0, 8.0, 7.9)-0.9)/10-0.23
lines (x,y)

text (1.75, 0.49, '***')
text (2.5, 0.46, 'NS')

dev.off()


##################################################################
###################################################################

layout(matrix(1))

a1 <- aov(lognxct ~ subtype, data_n_XCTs)
posthoc <- TukeyHSD(x=a1, 'subtype', conf.level=0.95)
summary(a1)
posthoc

t.test(lognxct ~ data_n_XCTs$subtype == 'Copy-number high (Serous-like)', data_n_XCTs)
t.test(lognxct~histological_type, data_n_XCTs)
a2 <- aov(lognxct ~ group_2, 
          subset(data_n_XCTs,
                 histological_type == 'Serous cystadenocarcinoma, NOS'))
posthoc2 <- TukeyHSD(x=a2, 'group_2', conf.level=0.95)
summary(a2)
posthoc2
t.test(lognxct ~ subset(data_n_XCTs,
                       histological_type == 'Serous cystadenocarcinoma, NOS')$group_2 == 
         'Preserved Xi', 
       subset(data_n_XCTs,
              histological_type == 'Serous cystadenocarcinoma, NOS'))


a3 <- aov(lognxct ~ group_2, 
          subset(data_n_XCTs,
                 histological_type == 'Endometrioid adenocarcinoma'))
posthoc3 <- TukeyHSD(x=a3, 'group_2', conf.level=0.95)
summary(a3)
posthoc3
t.test(lognxct ~ subset(data_n_XCTs,
                        histological_type == 'Endometrioid adenocarcinoma')$group_2 == 
         'Preserved Xi', 
       subset(data_n_XCTs,
              histological_type == 'Endometrioid adenocarcinoma'))


load('/home/jun/XIST_UT/methy.somatic.mean.RData')
methy_somatic <- data.frame(id = names(methy_somatic), methy_somatic)

somatic <- merge(data_table_group, methy_somatic, by ='id')



##################################################################
library(NMF)

svg(file = '/home/jun/XIST_UT/Figures/bp_somatic.svg',
    width =9, height = 4, pointsize = 12)
layout(matrix(c(1:3), ncol = 3, byrow = TRUE),
       widths = c(1,1,1), respect = FALSE)

par (mar=c(5.1, 5, 1, 1.1))
bp <- boxplot(methy_somatic~subtype_2, somatic,
              ylab = 'methylation of\nsomatic CpG island (beta value)',
              lab = '',
              #border = "white",
              frame = FALSE,
              names = rep('',4),
              xaxt='n',
              outpch = NA,
              ylim = c(0.25,0.5)
              
)

axis(mgp=c(3.5, 2, 0),
     side = 1,
     at = 1:4,
     pos = NA,
     tck = -0.01,
     labels=c('MSI',
              "POLE",
              "CN low",
              'CN high'))

stripchart(methy_somatic~subtype_2, somatic,
           vertical=T,pch=1, method="jitter",cex=1,add=T)

x <- c(1, 1, 3, 3)
y <- (c(7.9, 8.0, 8.0, 7.9)-0.9)/10-0.26
lines (x,y)
x <- c(2, 2, 4, 4)
y <- (c(7.9, 8.0, 8.0, 7.9)-0.9)/10-0.23
lines (x,y)

text (3, 0.49, '***')
text (2, 0.46, 'NS')

par (mar=c(5.1, 0, 1, 3.1))
bp <- boxplot(methy_somatic~histological_type, somatic,
              xaxt = 'n',
              yaxt = 'n',
              #ylab = 'mean RNA expression of CTAs (log RSEM)',
              outpch = NA,
              #border = "white",
              frame = FALSE,
              outpch = NA,
              ylim = c(0.25,0.5)
              
)

axis(mgp=c(3.5, 2, 0),
     side = 1,
     at = 1:2,
     pos = NA,
     tck = -0.01,
     labels=c('Enodmetrioid\nadenocarcinoma',
              "Serous\nadenocarcinoma"))

stripchart(methy_somatic~histological_type, somatic,
           vertical=T,pch=1, method="jitter",cex=1,add=T)

x <- c(1, 1, 2, 2)
y <- (c(7.9, 8.0, 8.0, 7.9)-0.9)/10-0.26
lines (x,y)
text (1.5, 0.46, '***')

x <- c(1.5, 1.5, 2.5, 2.5, 1.5)
y <- c(0.27, 0.42, 0.42, 0.27, 0.27)
lines (x,y)

par (mar=c(5.1, 1, 1, 1.1))
bp <- boxplot(methy_somatic ~ group_2, 
              subset(somatic,
                     histological_type == 'Serous cystadenocarcinoma, NOS'),
              #ylab = 'mean RNA expression of CTAs (log RSEM)',
              frame = FALSE,
              xaxt = 'n',
              yaxt = 'n',
              outpch = NA,
              ylim = c(0.25,0.5)
              )

axis(mgp = c(3.5,3,0),
     side = 1,
     at = 1:3,
     pos = NA,
     labels = c('Preserved Xi',
                'Partial\nreactivation\nof Xi',
                'Two Xa')
     
)

stripchart(methy_somatic ~ group_2, 
           subset(somatic,
                  histological_type == 'Serous cystadenocarcinoma, NOS'),
           vertical=T,pch=1, method="jitter",cex=1,add=T)


x <- c(2, 2, 3, 3)
y <- (c(7.9, 8.0, 8.0, 7.9)-0.9)/10-0.26
lines (x,y)

x <- c(1, 1, 2.5, 2.5)
y <- (c(7.9, 8.0, 8.0, 7.9)-0.9)/10-0.23
lines (x,y)

text (1.75, 0.49, '***')
text (2.5, 0.46, 'NS')



dev.off()


a1 <- aov(methy_somatic~subtype_2, somatic)
posthoc <- TukeyHSD(x=a1, 'subtype_2', conf.level=0.95)
summary(a1)
posthoc

t.test(methy_somatic ~ somatic$subtype_2 == 'Copy-number high (Serous-like)', somatic)
t.test(methy_somatic~histological_type, somatic)
a2 <- aov(methy_somatic ~ group_2, 
          subset(somatic,
                 histological_type == 'Serous cystadenocarcinoma, NOS'))
posthoc2 <- TukeyHSD(x=a2, 'group_2', conf.level=0.95)
summary(a2)
posthoc2


a3 <- aov(methy_somatic ~ group_2, 
          subset(somatic,
                 histological_type == 'Endometrioid adenocarcinoma'))
posthoc3 <- TukeyHSD(x=a3, 'group_2', conf.level=0.95)
summary(a3)
posthoc3


boxplot(methy_somatic~histological_type, somatic)

library(xtable)
prop <- round(prop.table(table(data_table_group$group_2, data_table_group$subtype), 2), 3)*100
write(print(xtable(prop), type = 'html'), file = './Figures/prop.html')
xtable(prop)
