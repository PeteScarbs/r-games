### Meta-analysis for systematic review of self-monitoring and tailored feedback in diet ###

# Last updated: 6th June 2017

# install.packages("metafor")
library("metafor")

# Open analysis dataset
setwd("S:\\Pete\\NUTRITION\\OTHER PROJECTS\\myShop\\FUNDING OPPORTUNITIES\\MRC PUBLIC HEALTH INTERVENTIONS\\SYSTEMATIC REVIEW OF SELF MONITORING")
data <- read.csv("Teasdale meta-analysis_06JUN2017.csv")

# Multi-level random effects meta-analysis with outcomes clustered in papers
meta.ml <- rma.mv(effect.size, meta.var, random = ~1|study.id, data = data, slab = paste(study.name, outcome, sep = ", "))

# Forest plot
par(mar=c(4,4,1,2), font = 1, ps = 16) # decrease margins so the full space is used
forest(meta.ml, xlim=c(-16, 4), at=c(-1, -0.5, 0, 0.5, 1), 
       ilab=cbind(data$int.n, round(data$int.effect, 1), round(data$int.sd, 1), data$con.n, round(data$con.effect, 1), round(data$con.sd, 1)),
       ilab.xpos=c(-8,-7,-6,-4, -3, -2), 
       xlab="     Favours control / Favours intervention", mlab="Random effects multi-level model (individuals and outcomes nested in studies)", psize=1)
# add column headings to the plot
par(font = 2, ps = 13)
text(c(-8,-7,-6,-4, -3, -2),         47, c("n", "Mean", "SD", "n", "Mean", "SD"))
par(font = 2, ps = 14)
text(c(-7,-3),                       48.5, c("Intervention", "Control"))
par(font = 2, ps = 13)
text(-16,                            47, "First author and Year",     pos=4)
text(2.88,                           48, "Standardised Mean")
text(4,                              47, "Difference [95% CI]", pos=2)

# Equivalent to funnel plot
# Need to do a multi-level regression, with outcomes nested in studies, with standard.error as a study-level predictor

# install.packages("lme4")
library(lme4)
funnel <- lmer(effect.size ~ standard.error + (1|study.id), data = data)
Vcov <- vcov(funnel, useScale = FALSE)
betas <- fixef(funnel)
se <- sqrt(diag(Vcov))
zval <- betas / se
pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
cbind(betas, se, zval, pval)
# This shows that effect size gets bigger as standard error gets bigger - indicating publication bias.
# But could be heavily influenced by Mummah, which is a massive outlier. So rerun without.

# Funnel results without Mummah, 2016
data2 <- rbind(data[1:27,],data[29:45,])
funnel2 <- lmer(effect.size ~ standard.error + (1|study.id), data = data2)
Vcov2 <- vcov(funnel2, useScale = FALSE)
betas2 <- fixef(funnel2)
se2 <- sqrt(diag(Vcov2))
zval2 <- betas2 / se2
pval2 <- 2 * pnorm(abs(zval2), lower.tail = FALSE)
cbind(betas2, se2, zval2, pval2)
# And no Standard error is no longer significantly associated with effect size!

# 'One study removed' sensitivity analysis.

new.study.id <- c(rep(1,1),  #1
                  rep(2,2),  #2
                  rep(3,3),  #5
                  rep(4,3),  #6
                  rep(5,3),  #7
                  rep(6,1),  #8
                  rep(7,2),  #11
                  rep(8,3),  #12
                  rep(9,1), #13
                  rep(10,1), #14
                  rep(11,4), #15
                  rep(12,2), #16
                  rep(13,1), #17
                  rep(14,1), #18
                  rep(15,3), #20
                  rep(16,1), #21
                  rep(17,4), #23
                  rep(18,3), #24
                  rep(19,1), #25
                  rep(20,1), #26
                  rep(21,4)) #28
data <- cbind(data, new.study.id)

for (i in 1:21){
  assign("data.x", data[data$new.study.id!=i,])
  assign(paste0("meta.ml.",i),rma.mv(effect.size, meta.var, random = ~1|study.id, 
                                     data = data.x, slab = paste(study.name, outcome, sep = ", ")))
}
sens.estimates <- c(meta.ml.1[1],
                    meta.ml.2[1],
                    meta.ml.3[1],
                    meta.ml.4[1],
                    meta.ml.5[1],
                    meta.ml.6[1],
                    meta.ml.7[1],
                    meta.ml.8[1],
                    meta.ml.9[1],
                    meta.ml.10[1],
                    meta.ml.11[1],
                    meta.ml.12[1],
                    meta.ml.13[1],
                    meta.ml.14[1],
                    meta.ml.15[1],
                    meta.ml.16[1],
                    meta.ml.17[1],
                    meta.ml.18[1],
                    meta.ml.19[1],
                    meta.ml.20[1],
                    meta.ml.21[1])
sens.se        <- c(meta.ml.1[2],
                    meta.ml.2[2],
                    meta.ml.3[2],
                    meta.ml.4[2],
                    meta.ml.5[2],
                    meta.ml.6[2],
                    meta.ml.7[2],
                    meta.ml.8[2],
                    meta.ml.9[2],
                    meta.ml.10[2],
                    meta.ml.11[2],
                    meta.ml.12[2],
                    meta.ml.13[2],
                    meta.ml.14[2],
                    meta.ml.15[2],
                    meta.ml.16[2],
                    meta.ml.17[2],
                    meta.ml.18[2],
                    meta.ml.19[2],
                    meta.ml.20[2],
                    meta.ml.21[2])
sens.lowci     <- c(meta.ml.1[5],
                    meta.ml.2[5],
                    meta.ml.3[5],
                    meta.ml.4[5],
                    meta.ml.5[5],
                    meta.ml.6[5],
                    meta.ml.7[5],
                    meta.ml.8[5],
                    meta.ml.9[5],
                    meta.ml.10[5],
                    meta.ml.11[5],
                    meta.ml.12[5],
                    meta.ml.13[5],
                    meta.ml.14[5],
                    meta.ml.15[5],
                    meta.ml.16[5],
                    meta.ml.17[5],
                    meta.ml.18[5],
                    meta.ml.19[5],
                    meta.ml.20[5],
                    meta.ml.21[5])
sens.hici      <- c(meta.ml.1[6],
                    meta.ml.2[6],
                    meta.ml.3[6],
                    meta.ml.4[6],
                    meta.ml.5[6],
                    meta.ml.6[6],
                    meta.ml.7[6],
                    meta.ml.8[6],
                    meta.ml.9[6],
                    meta.ml.10[6],
                    meta.ml.11[6],
                    meta.ml.12[6],
                    meta.ml.13[6],
                    meta.ml.14[6],
                    meta.ml.15[6],
                    meta.ml.16[6],
                    meta.ml.17[6],
                    meta.ml.18[6],
                    meta.ml.19[6],
                    meta.ml.20[6],
                    meta.ml.21[6])
sens.p         <- c(meta.ml.1[4],
                    meta.ml.2[4],
                    meta.ml.3[4],
                    meta.ml.4[4],
                    meta.ml.5[4],
                    meta.ml.6[4],
                    meta.ml.7[4],
                    meta.ml.8[4],
                    meta.ml.9[4],
                    meta.ml.10[4],
                    meta.ml.11[4],
                    meta.ml.12[4],
                    meta.ml.13[4],
                    meta.ml.14[4],
                    meta.ml.15[4],
                    meta.ml.16[4],
                    meta.ml.17[4],
                    meta.ml.18[4],
                    meta.ml.19[4],
                    meta.ml.20[4],
                    meta.ml.21[4])        
sens.analysis <- cbind(sens.estimates, sens.se, sens.lowci, sens.hici, sens.p) # Basically no difference to the results

# I2 heterogeneity
# code taken from http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
W <- diag(1/data$meta.var)
X <- model.matrix(meta.ml)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(meta.ml$sigma2) / (sum(meta.ml$sigma2) + (meta.ml$k-meta.ml$p)/sum(diag(P)))

# I2 heterogeneity with outlier removed
meta.ml2 <- rma.mv(effect.size, meta.var, random = ~1|study.id, data = data2, slab = paste(study.name, outcome, sep = ", "))
W <- diag(1/data2$meta.var)
X <- model.matrix(meta.ml2)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(meta.ml2$sigma2) / (sum(meta.ml2$sigma2) + (meta.ml2$k-meta.ml2$p)/sum(diag(P)))


# Meta-analyses stratified by risk of bias (0 = high risk of bias; 1 = low risk of bias)
meta.ml.lowbias <- rma.mv(effect.size, meta.var, random = ~1|study.id, data = data[data$bias==1,], slab = paste(study.name, outcome, sep = ", "))
meta.ml.hibias  <- rma.mv(effect.size, meta.var, random = ~1|study.id, data = data[data$bias==0,], slab = paste(study.name, outcome, sep = ", "))

# Forest plot - low bias
par(mar=c(4,4,1,2), font = 1, ps = 16) # decrease margins so the full space is used
forest(meta.ml.lowbias, xlim=c(-16, 4), at=c(-1, -0.5, 0, 0.5, 1), 
       ilab=cbind(data$int.n, round(data$int.effect, 1), round(data$int.sd, 1), data$con.n, round(data$con.effect, 1), round(data$con.sd, 1)),
       ilab.xpos=c(-8,-7,-6,-4, -3, -2), 
       xlab="     Favours control / Favours intervention", mlab="Random effects multi-level model (individuals and outcomes nested in studies)", psize=1)
# add column headings to the plot
par(font = 2, ps = 13)
text(c(-8,-7,-6,-4, -3, -2),         12, c("n", "Mean", "SD", "n", "Mean", "SD"))
par(font = 2, ps = 14)
text(c(-7,-3),                       13, c("Intervention", "Control"))
par(font = 2, ps = 13)
text(-16,                            12, "First author and Year",     pos=4)
text(2.88,                           13, "Standardised Mean")
text(4,                              12, "Difference [95% CI]", pos=2)

# Forest plot - high risk of bias
par(mar=c(4,4,1,2), font = 1, ps = 16) # decrease margins so the full space is used
forest(meta.ml.hibias, xlim=c(-16, 4), at=c(-1, -0.5, 0, 0.5, 1), 
       ilab=cbind(data$int.n, round(data$int.effect, 1), round(data$int.sd, 1), data$con.n, round(data$con.effect, 1), round(data$con.sd, 1)),
       ilab.xpos=c(-8,-7,-6,-4, -3, -2), 
       xlab="     Favours control / Favours intervention", mlab="Random effects multi-level model (individuals and outcomes nested in studies)", psize=1)
# add column headings to the plot
par(font = 2, ps = 13)
text(c(-8,-7,-6,-4, -3, -2),         37, c("n", "Mean", "SD", "n", "Mean", "SD"))
par(font = 2, ps = 14)
text(c(-7,-3),                       38, c("Intervention", "Control"))
par(font = 2, ps = 13)
text(-16,                            37, "First author and Year",     pos=4)
text(2.88,                           38, "Standardised Mean")
text(4,                              37, "Difference [95% CI]", pos=2)
