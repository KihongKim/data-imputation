# MICE with either CART or Random Forests

setwd("/home/kihong/Imputation") # on the 'sapporo' server
rm(list=ls())

library(RPostgreSQL)
library(plyr)
library(lattice)
library(mice)
library(CALIBERrfimpute)
library(missForest)
library(reshape)

### connect to the PostgreSQL database on the sapporo server
conn <- dbConnect(PostgreSQL(), host="sapporo.usp.pdx.edu", user="smartdata", password="xxx", dbname="portland")

### read taxlots20XX
taxlots2011 <- dbReadTable(conn, c("rlis","taxlots2011"))
#
tl <- taxlots2011
tl$rowname <- as.numeric(row.names(tl))
tl$district <- as.factor(tl$district)
tl$landuse <- as.factor(tl$landuse)
# delete tax lots that are not defined in block groups or districts
tl <- subset(tl, !is.na(district) & !is.na(bg_fips))

### define missing values
# 'area' and 'district' do not have missing values
# 'bldgsqft', 'bldgval', and 'landval' may have missing values coded zero
t1 <- ddply(tl, .(landuse), summarize, total=length(bldgsqft))
t2 <- ddply(tl[tl$bldgsqft==0,], .(landuse), summarize, zero.bldgsqft=length(bldgsqft))
t3 <- ddply(tl[tl$bldgval==0,], .(landuse), summarize, zero.bldgval=length(bldgval))
t4 <- ddply(tl[tl$landval==0,], .(landuse), summarize, zero.landval=length(landval))
t5 <- cbind(t1,t2[-1],t3[-1],t4[-1])
t5
#    landuse  total zero.bldgsqft zero.bldgval zero.landval
# 1      AGR  15418          7041         5067         2579
# 2      COM  21752          9458         1976         1232
# 3      FOR  13763          7828         6439         2625
# 4      IND   5661          3202          705          435
# 5      MFR  58280         18517         1336        53089
# 6      PUB   3793          3587         3470         1768
# 7      RUR  15345          3277          937           47
# 8      SFR 415716          8264         2662          296
# 9      VAC  43087         42655        42011         5520
# 10    <NA>   6774          5389         5235         5196

ddply(tl[!is.na(tl$far),], .(district), summarize, freq=length(far),min=min(far),mean=mean(far),max=max(far))
ddply(tl[!is.na(tl$bvpa),], .(district), summarize, freq=length(bvpa),min=min(bvpa),mean=mean(bvpa),max=max(bvpa))
ddply(tl[!is.na(tl$lvpa),], .(district), summarize, freq=length(lvpa),min=min(lvpa),mean=mean(lvpa),max=max(lvpa))

# define missing 
tl$bldgsqft[which(tl$landuse %in% c('COM','IND','MFR','SFR') & tl$bldgsqft==0)] <- NA
tl$bldgval[which(tl$landuse %in% c('COM','IND','MFR','SFR') & tl$bldgval==0)] <- NA
tl$landval[which(tl$landuse %in% c('COM','IND','MFR','SFR') & tl$landval==0)] <- NA

### transform variables of interest
tl$far <- tl$bldgsqft / tl$area		# building SQFT per parcel area
tl$bvpa <- tl$bldgval / tl$area		# building value per parcel area
tl$lvpa <- tl$landval / tl$area		# land value per parcel area

#
save(tl, file="tl.RData")

### Multnohma County
mult <- subset(tl, county=='051')
mult <- mult[c("area","district","bg_fips","pop10","job11","landuse","bldgsqft","bldgval","landval","far","bvpa","lvpa")]
save(mult, file="mult.RData")
length(which(is.na(mult$landuse))) / nrow(mult)
length(which(is.na(mult$bldgsqft))) / nrow(mult)
length(which(is.na(mult$bldgval))) / nrow(mult)
length(which(is.na(mult$landval))) / nrow(mult)

vars <- c("area","district","pop10","job11","landuse","far","bvpa","lvpa")
mult.mice <- mult[vars]
mult.mice$district <- factor(mult.mice$district)
mult.mice$landuse <- factor(mult.mice$landuse)
save(mult.mice, file="mult.mice.RData")
#
# vars2 <- c("area","district","pop10","job11","landuse","bldgsqft","bldgval","landval")
# mult.mice2 <- mult[vars2]
# mult.mice2$district <- factor(mult.mice2$district)
# mult.mice2$landuse <- factor(mult.mice2$landuse)
# save(mult.mice2, file="mult.mice2.RData")
#
md.pattern(mult.mice)
#
validateMICE <- function(x,y) {
	tl.mice$rowname <- as.numeric(row.names(tl.mice))
	imputed <- x[["imp"]][[y]]
	imputed$rowname <- as.numeric(row.names(imputed))
	merged <- merge(imputed, tl.mice[c(y,"rowname")], by="rowname")
	rmse.1 <- sqrt(mean((merged$'1' - merged[[y]])^2))
	rmse.2 <- sqrt(mean((merged$'2' - merged[[y]])^2))
	rmse.3 <- sqrt(mean((merged$'3' - merged[[y]])^2))
	rmse.4 <- sqrt(mean((merged$'4' - merged[[y]])^2))
	rmse.5 <- sqrt(mean((merged$'5' - merged[[y]])^2))
	r2.1 <- summary(lm(merged$'1' ~ merged[[y]]))$r.squared
	r2.2 <- summary(lm(merged$'2' ~ merged[[y]]))$r.squared
	r2.3 <- summary(lm(merged$'3' ~ merged[[y]]))$r.squared
	r2.4 <- summary(lm(merged$'4' ~ merged[[y]]))$r.squared
	r2.5 <- summary(lm(merged$'5' ~ merged[[y]]))$r.squared
	t <- matrix(c(rmse.1, rmse.2, rmse.3, rmse.4, rmse.5,
				  r2.1, r2.2, r2.3, r2.4, r2.5), 2, 5, byrow=TRUE)
	dimnames(t) = list(c("rmse", "r2"), c("mids_1", "mids_2", "mids_3", "mids_4", "mids_5")) 
	return(t)
}
validateMICE2 <- function(x,y) {
	tl.mice$rowname <- as.numeric(row.names(tl.mice))
	imputed <- x[["imp"]][[y]]
	imputed$rowname <- as.numeric(row.names(imputed))
	merged <- merge(imputed, tl.mice[c(y,"rowname")], by="rowname")
	# classification table (ct) for each imputation
	ct.1 <- as.matrix(table(factor(merged$'1', levels(merged[[y]])), merged[[y]]))
	ct.2 <- as.matrix(table(factor(merged$'2', levels(merged[[y]])), merged[[y]]))
	ct.3 <- as.matrix(table(factor(merged$'3', levels(merged[[y]])), merged[[y]]))
	ct.4 <- as.matrix(table(factor(merged$'4', levels(merged[[y]])), merged[[y]]))
	ct.5 <- as.matrix(table(factor(merged$'5', levels(merged[[y]])), merged[[y]]))
	# accuracy rate (ar) for each imputation
	ar.1 <- sum(diag(ct.1))/nrow(merged)
	ar.2 <- sum(diag(ct.2))/nrow(merged)
	ar.3 <- sum(diag(ct.3))/nrow(merged)
	ar.4 <- sum(diag(ct.4))/nrow(merged)
	ar.5 <- sum(diag(ct.5))/nrow(merged)
	#
	table.val <- matrix(c(ar.1, ar.2, ar.3, ar.4, ar.5), 1, 5, byrow=TRUE)
	dimnames(table.val) = list("accuracy rate", c("mice_1", "mice_2", "mice_3", "mice_4", "mice_5"))
	# print(ct.1)
	# print(ct.2)
	# print(ct.3)
	# print(ct.4)
	# print(ct.5)	
	return(table.val)
}
#
set.seed(1234)
mult.mice.val <- subset(mult.mice, !is.na(landuse) & !is.na(far) & !is.na(bvpa) & !is.na(lvpa))
mult.mice.val.missing <- mult.mice.val[c("landuse","far","bvpa","lvpa")]
mult.mice.val.missing <- prodNA(mult.mice.val.missing, noNA=0.05)
mult.mice.val <- cbind(mult.mice.val[,c(1,2,3,4)],mult.mice.val.missing)
save(mult.mice.val, file="mult.mice.val.RData")
#set.seed(1234)
#mult.mice2.val <- subset(mult.mice2, !is.na(landuse) & !is.na(bldgsqft) & !is.na(bldgval) & !is.na(landval))
#mult.mice2.val.missing <- mult.mice2.val[c("landuse","bldgsqft","bldgval","landval")]
#mult.mice2.val.missing <- prodNA(mult.mice2.val.missing, noNA=0.05)
#mult.mice2.val <- cbind(mult.mice2.val[,c(1,2,3,4)],mult.mice2.val.missing)
#save(mult.mice2.val, file="mult.mice2.val.RData")
# validation: pmm
load("mult.mice.val.RData")
system.time(imp.mult.mice.val.pmm <- mice(mult.mice.val, seed=1234, method="pmm"))
print(imp.mult.mice.val.pmm)
save(imp.mult.mice.val.pmm, file="imp.mult.mice.val.pmm.RData")
validateMICE(imp.mult.mice.val.pmm, "far")
validateMICE(imp.mult.mice.val.pmm, "bvpa")
validateMICE(imp.mult.mice.val.pmm, "lvpa")
validateMICE2(imp.mult.mice.val.pmm, "landuse")
# validation: cart
load("mult.mice.val.RData")
system.time(imp.mult.mice.val.cart <- mice(mult.mice.val, seed=1234, method="cart"))
print(imp.mult.mice.val.cart)
save(imp.mult.mice.val.cart, file="imp.mult.mice.val.cart.RData")
validateMICE(imp.mult.mice.val.cart, "far")
validateMICE(imp.mult.mice.val.cart, "bvpa")
validateMICE(imp.mult.mice.val.cart, "lvpa")
validateMICE2(imp.mult.mice.val.cart, "landuse")
# validation: rf
load("mult.mice.val.RData")
system.time(imp.mult.mice.val.rf <- mice(mult.mice.val, seed=1234, method=c("","","","","rfcat","rfcont","rfcont","rfcont")))
print(imp.mult.mice.val.rf)
save(imp.mult.mice.val.rf, file="imp.mult.mice.val.rf.RData")
validateMICE(imp.mult.mice.val.rf, "far")
validateMICE(imp.mult.mice.val.rf, "bvpa")
validateMICE(imp.mult.mice.val.rf, "lvpa")
validateMICE2(imp.mult.mice.val.rf, "landuse")
#


# pmm
load("mult.mice.RData")
system.time(imp.mult.mice.pmm <- mice(mult.mice, seed=1234, method="pmm"))
print(imp.mult.mice.pmm)
save(imp.mult.mice.pmm, file="imp.mult.mice.pmm.RData")
# cart
load("mult.mice.RData")
system.time(imp.mult.mice.cart <- mice(mult.mice, seed=1234, method="cart"))
print(imp.mult.mice.cart)
save(imp.mult.mice.cart, file="imp.mult.mice.cart.RData")
# rf
load("mult.mice.RData")
system.time(imp.mult.mice.rf <- mice(mult.mice, seed=1234, method=c("","","","","rfcat","rfcont","rfcont","rfcont")))
print(imp.mult.mice.rf)
save(imp.mult.mice.rf, file="imp.mult.mice.rf.RData")

#### diagnostic
# stripplot
png("stripplot.png")
stripplot(imp.mult.mice.cart, far+bvpa+lvpa+landuse~.imp, pch=20, cex=1.2)
dev.off()

png("stripplot.far.png")
stripplot(imp.mult.mice.cart, far~.imp, ylab="FAR", pch=20, cex=1.2)
dev.off()

png("stripplot.bvpa.png")
stripplot(imp.mult.mice.cart, bvpa~.imp, ylab="BVPA", pch=20, cex=1.2)
dev.off()

png("stripplot.lvpa.png")
stripplot(imp.mult.mice.cart, lvpa~.imp, ylab="LVPA", pch=20, cex=1.2)
dev.off()

png("stripplot.png")
par(mfrow=c(1,3))
stripplot(imp.mult.mice.cart, far~.imp, ylab="FAR", pch=20, cex=1.2)
stripplot(imp.mult.mice.cart, bvpa~.imp, ylab="BVPA", pch=20, cex=1.2)
stripplot(imp.mult.mice.cart, lvpa~.imp, ylab="LVPA", pch=20, cex=1.2)
dev.off()

stripplot(imp.mult.mice.cart, landuse~.imp, ylab="LU", pch=20, cex=1.2)


# scatterplot
png("sp.far.bvpa.png"); xyplot(imp.mult.mice.cart, far ~ bvpa | .imp, pch = 20, cex = 1.4); dev.off()
png("sp.far.lvpa.png"); xyplot(imp.mult.mice.cart, far ~ lvpa | .imp, pch = 20, cex = 1.4); dev.off()
png("sp.bvpa.lvpa.png"); xyplot(imp.mult.mice.cart, bvpa ~ lvpa | .imp, pch = 20, cex = 1.4); dev.off()
# densityplot
png("dp.far.png"); densityplot(imp.mult.mice.cart, ~far | .imp); dev.off()
png("dp.bvpa.png"); densityplot(imp.mult.mice.cart, ~bvpa | .imp); dev.off()
png("dp.lvpa.png"); densityplot(imp.mult.mice.cart, ~lvpa | .imp); dev.off()



##### DQIs
#
load("mult.RData")
#mult.res <- subset(mult, landuse %in% c("SFR","MFR"))
mult.nonres <- subset(mult, !(landuse %in% c("SFR","MFR")))
mult.nonres.bg <- ddply(mult.nonres, .(bg_fips), summarise, bldgsqft_bg=sum(bldgsqft, na.rm=TRUE), job11_bg=min(job11))
r2.0 <- summary(lm(bldgsqft_bg ~ job11_bg, data=mult.nonres.bg))$r.squared
round(r2.0,3)
#
load("imp.mult.mice.pmm.RData")
load("imp.mult.mice.cart.RData")
load("imp.mult.mice.rf.RData")
computeDQI <- function(x,m) {
	for (i in 1:m) {
		comp <- complete(x, i)
		comp <- cbind(comp, mult[c("bg_fips")])
		comp$bldgsqft <- with(comp, area*far)
		comp.nonres <- subset(comp, !(landuse %in% c("SFR","MFR")))
		comp.nonres.bg <- ddply(comp.nonres, .(bg_fips), summarise, bldgsqft_bg=sum(bldgsqft, na.rm=TRUE), job11_bg=min(job11))
		fit <- lm(bldgsqft_bg ~ job11_bg, data=comp.nonres.bg)
		print(round(summary(fit)$r.squared,3))
	}
}
computeDQI(imp.mult.mice.pmm, 5)
computeDQI(imp.mult.mice.cart, 5)
computeDQI(imp.mult.mice.rf, 5)

#
DQIs <- data.frame(R2.0,R2.1,R2.2,R2.3,R2.4,R2.5)


##### regression analysis of building value
load("mult.mice.RData")
load("imp.mult.mice.pmm.RData")
load("imp.mult.mice.cart.RData")
load("imp.mult.mice.rf.RData")
library(MASS)
#
fm <- as.formula(bvpa ~ far + lvpa + I(district=="1 ") + I(landuse=="COM"))
#
fit.orig <- lm(bvpa ~ far + lvpa + I(district=="1 ")+I(landuse=="MFR"), data=mult.mice)
summary(fit.orig)
fit.orig2 <- lm(bvpa ~ far + lvpa + I(district=="1 ")*I(landuse=="MFR"), data=mult.mice)
summary(fit.orig2)
fit.orig3 <- lm(bvpa ~ far + lvpa + district + landuse, data=mult.mice)
summary(fit.orig3)
#
fit.cart <- with(imp.mult.mice.cart, lm(bvpa ~ far + lvpa + I(district=="1 ")+I(landuse=="MFR")))
round(summary(pool(fit.cart)), 3)
fit.cart2 <- with(imp.mult.mice.cart, lm(bvpa ~ far + lvpa + I(district=="1 ")*I(landuse=="MFR")))
round(summary(pool(fit.cart2)), 3)
#
fit.pmm <- with(imp.mult.mice.pmm, lm(bvpa ~ far + lvpa + I(district=="1 ") + I(landuse=="COM")))
round(summary(pool(fit.pmm)),3)
#
fit.rf <- with(imp.mult.mice.rf, lm(bvpa ~ far + lvpa + I(district=="1 ") + I(landuse=="COM")))
round(summary(pool(fit.rf)), 3)
#
robustfit.cart <- with(imp.mult.mice.cart, rlm(bvpa ~ far + lvpa + I(district=="1 ") + I(landuse=="COM"), maxit=100))
round(summary(pool(robustfit.cart)),3)



#summary(fit.cart$ana[[3]])


fit.cart <- with(imp.mult.mice.cart, lm(lvpa ~ I(district=="1 "):bvpa)) 

### create a data set for MICE
vars <- c("area","district","pop10","job11","landuse","far","bvpa","lvpa")
tl.mice <- tl[vars]
save(tl.mice, file="tl.mice.RData")

### inspect the pattern of missing values across variables
md.pattern(tl.mice)

### validation

## create a validation data set
set.seed(1234)
# draw a 5% random sample for each variable and replace them with NA
# tl.mice.val <- subset(tl.mice, !is.na(landuse) & !is.na(far) & !is.na(bvpa) & !is.na(lvpa))
# tl.mice.val$landuse[sample(nrow(tl.mice.val), 0.05*nrow(tl.mice.val))] <- NA
# tl.mice.val$far[sample(nrow(tl.mice.val), 0.05*nrow(tl.mice.val))] <- NA
# tl.mice.val$bvpa[sample(nrow(tl.mice.val), 0.05*nrow(tl.mice.val))] <- NA
# tl.mice.val$lvpa[sample(nrow(tl.mice.val), 0.05*nrow(tl.mice.val))] <- NA
# draw a 5% randome sample for all missing variables and replace them with NA using the prodNA function of the missForest package
tl.mice.val <- subset(tl.mice, !is.na(landuse) & !is.na(far) & !is.na(bvpa) & !is.na(lvpa))
tl.mice.val.missing <- tl.mice.val[c("landuse","far","bvpa","lvpa")]
tl.mice.val.missing <- prodNA(tl.mice.val.missing, noNA=0.05) # missing at random (MAR)
tl.mice.val <- cbind(tl.mice.val[,c(1,2,3,4)],tl.mice.val.missing)

# SFR
SFR <- subset(tl.mice, landuse=="SFR")
SFR <- SFR[c(-5)]
SFR.val <- subset(SFR, !is.na(far) & !is.na(bvpa) & !is.na(lvpa))
SFR.val.missing <- SFR.val[c("far","bvpa","lvpa")]
SFR.val.missing <- prodNA(SFR.val.missing, noNA=0.05)
SFR.val <- cbind(SFR.val[,c(1,2,3,4)], SFR.val.missing)
save(SFR.val, file="SFR.val.RData")
system.time(imp.SFR.val.rf <- mice(SFR.val, seed=1234, method="rfcont", visitSequence="monotone"))
print(imp.SFR.val.rf)
save(imp.SFR.val.rf, file="imp.SFR.val.rf.RData")
system.time(imp.SFR.val.cart <- mice(SFR.val, seed=1234, method="cart", visitSequence="monotone"))
print(imp.SFR.val.cart)
save(imp.SFR.val.cart, file="imp.SFR.val.cart.RData")

# MFR
MFR <- subset(tl.mice, landuse=="MFR")
MFR <- MFR[c(-5)]
MFR.val <- subset(MFR, !is.na(far) & !is.na(bvpa) & !is.na(lvpa))
MFR.val.missing <- MFR.val[c("far","bvpa","lvpa")]
MFR.val.missing <- prodNA(MFR.val.missing, noNA=0.05)
MFR.val <- cbind(MFR.val[,c(1,2,3,4)], MFR.val.missing)
save(MFR.val, file="MFR.val.RData")
system.time(imp.MFR.val.rf <- mice(MFR.val, seed=1234, method="rfcont", visitSequence="monotone"))
print(imp.MFR.val.rf)
save(imp.MFR.val.rf, file="imp.MFR.val.rf.RData")
system.time(imp.MFR.val.cart <- mice(MFR.val, seed=1234, method="cart", visitSequence="monotone"))
print(imp.MFR.val.cart)
save(imp.MFR.val.cart, file="imp.MFR.val.cart.RData")

# COM
COM <- subset(tl.mice, landuse=="COM")
COM <- COM[c(-5)]
COM.val <- subset(COM, !is.na(far) & !is.na(bvpa) & !is.na(lvpa))
COM.val.missing <- COM.val[c("far","bvpa","lvpa")]
COM.val.missing <- prodNA(COM.val.missing, noNA=0.05)
COM.val <- cbind(COM.val[,c(1,2,3,4)], COM.val.missing)
save(COM.val, file="COM.val.RData")
system.time(imp.COM.val.rf <- mice(COM.val, seed=1234, method="rfcont", visitSequence="monotone"))
#print(imp.COM.val.rf)
save(imp.COM.val.rf, file="imp.COM.val.rf.RData")
system.time(imp.COM.val.cart <- mice(COM.val, seed=1234, method="cart", visitSequence="monotone"))
#print(imp.COM.val.cart)
save(imp.COM.val.cart, file="imp.COM.val.cart.RData")
#
test<-validateMICE(imp.COM.val.cart, "far")
validateMICE(imp.COM.val.cart, "bvpa")
validateMICE(imp.COM.val.cart, "lvpa")
validateMICE(imp.COM.val.rf, "far")
validateMICE(imp.COM.val.rf, "bvpa")
validateMICE(imp.COM.val.rf, "lvpa")

### IND
IND <- subset(tl.mice, landuse=="IND")
IND <- IND[c(-5)]
IND.val <- subset(IND, !is.na(far) & !is.na(bvpa) & !is.na(lvpa))
IND.val.missing <- IND.val[c("far","bvpa","lvpa")]
IND.val.missing <- prodNA(IND.val.missing, noNA=0.05)
IND.val <- cbind(IND.val[,c(1,2,3,4)], IND.val.missing)
save(IND.val, file="IND.val.RData")
# rf
system.time(imp.IND.val.rf <- mice(IND.val, seed=1234, method="rfcont", visitSequence="monotone", ntree=1000, mtry=ncol(IND.val)))
print(imp.IND.val.rf)
save(imp.IND.val.rf, file="imp.IND.val.rf.RData")
# cart
system.time(imp.IND.val.cart <- mice(IND.val, seed=1234, method="cart", visitSequence="monotone"))
print(imp.IND.val.cart)
save(imp.IND.val.cart, file="imp.IND.val.cart.RData")
# pmm
system.time(imp.IND.val.pmm <- mice(IND.val, seed=1234, method="pmm", visitSequence="monotone"))
print(imp.IND.val.pmm)
save(imp.IND.val.pmm, file="imp.IND.val.pmm.RData")
#
validateMICE(imp.IND.val.pmm, "far")
validateMICE(imp.IND.val.pmm, "bvpa")
validateMICE(imp.IND.val.pmm, "lvpa")
validateMICE(imp.IND.val.cart, "far")
validateMICE(imp.IND.val.cart, "bvpa")
validateMICE(imp.IND.val.cart, "lvpa")
validateMICE(imp.IND.val.rf, "far")
validateMICE(imp.IND.val.rf, "bvpa")
validateMICE(imp.IND.val.rf, "lvpa")




#
save(tl.mice.val, file="tl.mice.val.RData")



## inspect the pattern of missing values across variables in the validation data set
md.pattern(tl.mice.val)

## diagnostics
ddply(tl.mice.val2[!is.na(tl.mice.val2$far),], .(district), summarize, freq=length(far),min=min(far),mean=mean(far),max=max(far))
ddply(tl.mice.val2[is.na(tl.mice.val2$far),], .(district), summarize, freq=length(far))
ddply(tl.mice.val2[!is.na(tl.mice.val2$bvpa),], .(district), summarize, freq=length(bvpa),min=min(bvpa),mean=mean(bvpa),max=max(bvpa))
ddply(tl.mice.val2[is.na(tl.mice.val2$bvpa),], .(district), summarize, freq=length(bvpa))
ddply(tl.mice.val2[!is.na(tl.mice.val2$lvpa),], .(district), summarize, freq=length(lvpa),min=min(lvpa),mean=mean(lvpa),max=max(lvpa))
ddply(tl.mice.val2[is.na(tl.mice.val2$lvpa),], .(district), summarize, freq=length(lvpa))


## conduct MICE

# MICE via CART using the mice package with the validation data set 1
system.time(imp.val.cart <- mice(tl.mice.val, seed=1234, method="cart"))
print(imp.val.cart)
save(imp.val.cart, file="imp.val.cart.RData")
# MICE via RF using the mice package with the validation data set 1
system.time(imp.val.rf <- mice(tl.mice.val, seed=1234, method="rf", ntree=5))
print(imp.val.rf)
save(imp.val.rf, file="imp.val.rf.RData")
# MICE via RF using the CALIBERrfimpute package with the validation data set 1
system.time(imp.val.CALIBERrf <- mice(tl.mice.val, seed=1234, method=c("","","","rfcat","rfcont","rfcont","rfcont")))
print(imp.val.CALIBERrf)
save(imp.val.CALIBERrf, file="imp.val.CALIBERrf.RData")

# MICE via CART using the mice package with the validation data set 2
system.time(imp.val2.cart <- mice(tl.mice.val2, seed=1234, method="cart"))
print(imp.val2.cart)
save(imp.val2.cart, file="imp.val2.cart.RData")
# MICE via RF using the CALIBERrfimpute package with the validation data set 2
system.time(imp.val2.CALIBERrf <- mice(tl.mice.val2, seed=1234, method=c("","","","","","rfcat","rfcont","rfcont","rfcont")))
print(imp.val2.CALIBERrf)
save(imp.val2.CALIBERrf, file="imp.val2.CALIBERrf.RData")

## create validation results

# create a function to merge imputed and observed values for each variable in the validation data set
merge.imp.obs <- function(x, y) {
	# x: a data set that includes imputed values in the validation data set
	# y: a missing variable
	tl.mice$rowname <- as.numeric(row.names(tl.mice))
	imputed <- x[["imp"]][[y]]
	imputed$rowname <- as.numeric(row.names(imputed))
	imp.obs <- merge(imputed, tl.mice[c(y,"rowname")], by="rowname")
	return(imp.obs)
}
# create a function to validate MICE results for a continuous variable
validateMICE <- function(x,y) {
	tl.mice$rowname <- as.numeric(row.names(tl.mice))
	imputed <- x[["imp"]][[y]]
	imputed$rowname <- as.numeric(row.names(imputed))
	merged <- merge(imputed, tl.mice[c(y,"rowname")], by="rowname")
	#
	observed <- x[deparse(substitute(x))]
	names(observed) <- "observed"
	x <- cbind(x, observed)
	#
	rmse.1 <- sqrt(mean((x$'1' - x$observed)^2))
	rmse.2 <- sqrt(mean((x$'2' - x$observed)^2))
	rmse.3 <- sqrt(mean((x$'3' - x$observed)^2))
	rmse.4 <- sqrt(mean((x$'4' - x$observed)^2))
	rmse.5 <- sqrt(mean((x$'5' - x$observed)^2))
	r2.1 <- summary(lm(x$'1' ~ x$observed))$r.squared
	r2.2 <- summary(lm(x$'2' ~ x$observed))$r.squared
	r2.3 <- summary(lm(x$'3' ~ x$observed))$r.squared
	r2.4 <- summary(lm(x$'4' ~ x$observed))$r.squared
	r2.5 <- summary(lm(x$'5' ~ x$observed))$r.squared
	table.val <- matrix(c(rmse.1, rmse.2, rmse.3, rmse.4, rmse.5,
						  r2.1, r2.2, r2.3, r2.4, r2.5), 2, 5, byrow=TRUE)
	dimnames(table.val) = list(c("rmse", "r2"), c("mids_1", "mids_2", "mids_3", "mids_4", "mids_5")) 
	return(table.val)
}
# create a function to validate MICE results for a categorical variable
validateMICE2 <- function(x) {
	#
	observed <- x[deparse(substitute(x))]
	names(observed) <- "observed"
	x <- cbind(x, observed)
	# classification table (ct) for each imputation
	ct.1 <- as.matrix(table(as.factor(x$'1'), x$observed))
	ct.2 <- as.matrix(table(as.factor(x$'2'), x$observed))
	ct.3 <- as.matrix(table(as.factor(x$'3'), x$observed))
	ct.4 <- as.matrix(table(as.factor(x$'4'), x$observed))
	ct.5 <- as.matrix(table(as.factor(x$'5'), x$observed))
	# accuracy rate (ar) for each imputation
	ar.1 <- sum(diag(ct.1))/nrow(x)
	ar.2 <- sum(diag(ct.2))/nrow(x)
	ar.3 <- sum(diag(ct.3))/nrow(x)
	ar.4 <- sum(diag(ct.4))/nrow(x)
	ar.5 <- sum(diag(ct.5))/nrow(x)
	#
	table.val <- matrix(c(ar.1, ar.2, ar.3, ar.4, ar.5), 1, 5, byrow=TRUE)
	dimnames(table.val) = list("accuracy rate", c("mice_1", "mice_2", "mice_3", "mice_4", "mice_5"))
	print(ct.1)
	print(ct.2)
	print(ct.3)
	print(ct.4)
	print(ct.5)	
	return(table.val)
}

# apply the function of merge.imp.obs
landuse <- merge.imp.obs(imp.tl.mice.val.cart, "landuse")
far <- merge.imp.obs(imp.tl.mice.val.cart, "far")
bvpa <- merge.imp.obs(imp.tl.mice.val.cart, "bvpa")
lvpa <- merge.imp.obs(imp.tl.mice.val.cart, "lvpa")

landuse <- merge.imp.obs(imp.val.CALIBERrf, "landuse")
far <- merge.imp.obs(imp.val.CALIBERrf, "far")
bvpa <- merge.imp.obs(imp.val.CALIBERrf, "bvpa")
lvpa <- merge.imp.obs(imp.val.CALIBERrf, "lvpa")

landuse <- merge.imp.obs(imp.val2.cart, "landuse")
far <- merge.imp.obs(imp.val2.cart, "far")
bvpa <- merge.imp.obs(imp.val2.cart, "bvpa")
lvpa <- merge.imp.obs(imp.val2.cart, "lvpa")

landuse <- merge.imp.obs(imp.val2.CALIBERrf, "landuse")
far <- merge.imp.obs(imp.val2.CALIBERrf, "far")
bvpa <- merge.imp.obs(imp.val2.CALIBERrf, "bvpa")
lvpa <- merge.imp.obs(imp.val2.CALIBERrf, "lvpa")

# apply the function of validateMICE1/2
validateMICE2(landuse)
validateMICE1(far)
validateMICE1(bvpa)
validateMICE1(lvpa)

## validate graphically
# scatterplot between imputed and observed for each missing variable
png("scatter.val.pmm.png")
par(mfrow=c(3,5)


# far
png("scatter.far.png")
par(mfrow=c(3,2))
plot(far$far,far$"1")
plot(far$far,far$"2")
plot(far$far,far$"3")
plot(far$far,far$"4")
plot(far$far,far$"5")
dev.off()
# bvpa
png("scatter.bvpa.png")
par(mfrow=c(3,2))
plot(bvpa$bvpa,bvpa$"1")
plot(bvpa$bvpa,bvpa$"2")
plot(bvpa$bvpa,bvpa$"3")
plot(bvpa$bvpa,bvpa$"4")
plot(bvpa$bvpa,bvpa$"5")
dev.off()
# lvpa
png("scatter.lvpa.png")
par(mfrow=c(3,2))
plot(lvpa$lvpa,lvpa$"1")
plot(lvpa$lvpa,lvpa$"2")
plot(lvpa$lvpa,lvpa$"3")
plot(lvpa$lvpa,lvpa$"4")
plot(lvpa$lvpa,lvpa$"5")
dev.off()

# histogram
png("histograms.png")
par(mfrow=c(1,3))
hist(tl.mice.val2$far)
hist(tl.mice.val2$bvpa)
hist(tl.mice.val2$lvpa)
dev.off()

# stripplot to see the distribution of imputed values over observed values
png("stripplot.imp.val.cart.png"); stripplot(imp.val.cart, pch=20, cex=1.2); dev.off()
png("stripplot.imp.val2.cart.png"); stripplot(imp.val2.cart, pch=20, cex=1.2); dev.off()


png("stripplot.imp.val.cart.png"); stripplot(imp.val.cart, pch=20, cex=1.2); dev.off()

png("test1.png");plot(far$far,far$"1"); dev.off()
png("test2.png");plot(far$far,far$"2"); dev.off()
png("test3.png");plot(far$far,far$"3"); dev.off()
png("test4.png");plot(far$far,far$"4"); dev.off()
png("test5.png");plot(far$far,far$"5"); dev.off()

# xyplot
png("xyplot.imp.val.far.1.png"); xyplot(imp.val, far~bvpa | .imp, pch=20, cex=1.4); dev.off()
png("xyplot.imp.val.far.2.png"); xyplot(imp.val, far~lvpa | .imp, pch=20, cex=1.4); dev.off()
# png("xyplot.imp.val.far.3.png"); xyplot(imp.val, far~juris_city | .imp, pch=20, cex=1.4); dev.off()
png("xyplot.imp.val.far.4.png"); xyplot(imp.val, far~landuse | .imp, pch=20, cex=1.4); dev.off()

png("xyplot.imp.val.bvpa.1.png"); xyplot(imp.val, bvpa~bvpa | .imp, pch=20, cex=1.4); dev.off()
png("xyplot.imp.val.bvpa.2.png"); xyplot(imp.val, bvpa~lvpa | .imp, pch=20, cex=1.4); dev.off()
# png("xyplot.imp.val.bvpa.3.png"); xyplot(imp.val, bvpa~juris_city | .imp, pch=20, cex=1.4); dev.off()
png("xyplot.imp.val.bvpa.4.png"); xyplot(imp.val, bvpa~landuse | .imp, pch=20, cex=1.4); dev.off()

png("xyplot.imp.val.lvpa.1.png"); xyplot(imp.val, lvpa~bvpa | .imp, pch=20, cex=1.4); dev.off()
png("xyplot.imp.val.lvpa.2.png"); xyplot(imp.val, lvpa~lvpa | .imp, pch=20, cex=1.4); dev.off()
# png("xyplot.imp.val.lvpa.3.png"); xyplot(imp.val, lvpa~juris_city | .imp, pch=20, cex=1.4); dev.off()
png("xyplot.imp.val.lvpa.4.png"); xyplot(imp.val, lvpa~landuse | .imp, pch=20, cex=1.4); dev.off()

# densityplot
png("densityplot.imp.val.far.png"); densityplot(imp.val, ~far | .imp); dev.off()
png("densityplot.imp.val.bvpa.png"); densityplot(imp.val, ~bvpa | .imp); dev.off()
png("densityplot.imp.val.lvpa.png"); densityplot(imp.val, ~lvpa | .imp); dev.off()


########### conductin MICE with CART in the whole data set
setwd("/home/kihong/Imputation")
#load("tl.RData")
load("tl.mice.RData")
load("tl.mice.val.RData")
library(mice)
library(CALIBERrfimpute)
#
system.time(imp.tl.mice.cart <- mice(tl.mice, method="cart", cp=0.01, visitSequence='monotone', seed=1234))
#  
print(imp.tl.mice.cart)
save(imp.tl.mice.cart, file="imp.tl.mice.cart.RData")
#
system.time(imp.tl.mice.val.cart <- mice(tl.mice.val, method="cart", visitSequence='monotone', seed=1234))
# minbucket=100, cp=0.001, 
print(imp.tl.mice.val.cart)
save(imp.tl.mice.val.cart, file="imp.tl.mice.val.cart.RData")

# minbucket: the minimum number of observations in any terminal node; minsplit = minbucket*3
# cp: complexity parameter. Any split that does not decrease the overall lack of ﬁt by a factor of cp is not attempted.

#
system.time(imp.tl.mice.val.pmm <- mice(tl.mice.val, method="pmm", seed=1234))
print(imp.tl.mice.val.pmm)
save(imp.tl.mice.val.pmm, file="imp.tl.mice.val.pmm.RData")

# with a subset including tax lots of COM, IND, SFR, and MFR
tl.mice2 <- subset(tl.mice, landuse %in% c('COM','IND','MFR','SFR'))
tl.mice2$landuse <- factor(tl.mice2$landuse)
tl.mice2 <- tl.mice2[c(-3,-4)]
system.time(imp.tl.mice2.cart <- mice(tl.mice2, method="cart", visitSequence='monotone', seed=1234))
print(imp.tl.mice2.cart)
save(imp.tl.mice2.cart, file="imp.tl.mice2.cart.RData")


tl.mice$cbd <- factor(ifelse(as.numeric(tl.mice$district)==1, 1, 0))
tl.mice3 <- tl.mice[c(-2,-3,-4)]

system.time(imp.tl.mice3.cart <- mice(tl.mice3, method="cart", visitSequence='monotone', seed=1234))
print(imp.tl.mice3.cart)
save(imp.tl.mice3.cart, file="imp.tl.mice3.cart.RData")

SFR <- subset(tl.mice, landuse=="SFR")
SFR <- SFR[c(-5)]
system.time(imp.SFR.cart <- mice(SFR, method="cart", visitSequence='monotone', seed=1234))
print(imp.SFR.cart)
save(imp.SFR.cart, file="imp.SFR.cart.RData")
system.time(imp.SFR.rf <- mice(SFR, method="rfcont", visitSequence='monotone', seed=1234))
print(imp.SFR.rf)
save(imp.SFR.rf, file="imp.SFR.rf.RData")


print(imp.cart)
save(imp.cart, file="imp.cart.RData")

system.time(imp.tl.mice6 <- mice(tl.mice6, seed=1234, method="cart"))
print(imp.tl.mice6.cart)
save(imp.tl.mice6.cart, file="imp.tl.mice6.cart.RData")




#
# tl$d_VAC <- with(tl, ifelse(landuse=='VAC', 'vacant', 'occupied'))
# tl$d_bldgsqft <- with(tl, ifelse(bldgsqft==0, 'zero', 'non-zero'))
# tl$d_bldgval <- with(tl, ifelse(bldgval==0, 'zero', 'non-zero'))
# table(tl$d_VAC,tl$d_bldgsqft)
# table(tl$d_VAC,tl$d_bldgval)


### to Prof. van Buuren
parcel <- subset(tl.mice, !is.na(landuse) & !is.na(far) & !is.na(bvpa) & !is.na(lvpa))
set.seed(1234)
parcel$landuse[sample(nrow(tl.mice.val), 0.05*nrow(tl.mice.val))] <- NA
parcel$far[sample(nrow(tl.mice.val), 0.05*nrow(tl.mice.val))] <- NA
parcel$bvpa[sample(nrow(tl.mice.val), 0.05*nrow(tl.mice.val))] <- NA
parcel$lvpa[sample(nrow(tl.mice.val), 0.05*nrow(tl.mice.val))] <- NA

parcel.ext <- cbind(parcel, landuse.1.area = NA)
ini <- mice(parcel.ext, maxit=0, print=F)
meth <- ini$method
meth[c("landuse","far","bvpa","lvpa")] <- "cart"
meth["landuse.1.area"] <- "~I(landuse.1*area)"
pred <- ini$predictorMatrix
pred[c("far","bvpa","lvpa"), "landuse.1.area"] <- 1
pred[, "landuse"] <- 0
pred[, "area"] <- 0
pred["landuse","area"] <- 1
imp <- mice(parcel.ext, seed=1234, method=meth, predictorMatrix=pred)


# a dry-run
ini.val <- mice(tl.mice.val, maxit=0, print=F)
ini.val$method
ini.val$predictorMatrix
ini.val$visitSequence








### add interaction terms
tl.mice.val$COM.area <- with(tl.mice.val, ifelse(landuse=='COM', area, 0))
tl.mice.val$IND.area <- with(tl.mice.val, ifelse(landuse=='IND', area, 0))
tl.mice.val$MFR.area <- with(tl.mice.val, ifelse(landuse=='MFR', area, 0))
tl.mice.val$SFR.area <- with(tl.mice.val, ifelse(landuse=='SFR', area, 0))

tl.mice.val <- cbind(tl.mice.val,
					COM.area = NA,
					IND.area = NA,
					MFR.area = NA,
					SFR.area = NA)


### a dry-run
ini <- mice(tl.mice.val, maxit=0, print=F)

meth <- ini$method
meth[c("landuse","far","bvpa","lvpa")] <- "cart"
meth["COM.area"] <- "~I(landuse.1*area)"
meth["IND.area"] <- "~I(landuse.3*area)"
meth["MFR.area"] <- "~I(landuse.4*area)"
meth["SFR.area"] <- "~I(landuse.7*area)"

pred <- ini$predictorMatrix
pred[c("far","bvpa","lvpa"), c("COM.area","IND.area","MFR.area","SFR.area")] <- 1
pred[, "landuse"] <- 0
pred[, "area"] <- 0
pred["landuse","area"] <- 1
# pred[c("COM.area","IND.area","MFR.area","SFR.area"), c("far","bvpa","lvpa")] <- 1

#system.time(imp <- mice(tl.mice.val, seed=1234, method=meth, predictorMatrix=pred))
imp <- mice(tl.mice.val, seed=1234, method=meth, predictorMatrix=pred)
print(imp)






tl.mice.val <- cbind(tl.mice.val,
					landuse.1.area = NA,
					landuse.2.area = NA,
					landuse.3.area = NA,
					landuse.4.area = NA,
					landuse.5.area = NA,
					landuse.6.area = NA,
					landuse.7.area = NA,
					landuse.8.area = NA)
					

meth <- ini$method
meth[c("landuse","far","bvpa","lvpa")] <- "cart"
meth["landuse.1.area"] <- "˜I(landuse.1*area)"
meth["landuse.2.area"] <- "˜I(landuse.2*area)"
meth["landuse.3.area"] <- "˜I(landuse.3*area)"
meth["landuse.4.area"] <- "˜I(landuse.4*area)"
meth["landuse.5.area"] <- "˜I(landuse.5*area)"
meth["landuse.6.area"] <- "˜I(landuse.6*area)"
meth["landuse.7.area"] <- "˜I(landuse.7*area)"
meth["landuse.8.area"] <- "˜I(landuse.8*area)"





pred <- ini$predictorMatrix
pred[, c("area","landuse")] <- 0
pred["landuse", "area"] <- 1
pred[c("far","bvpa","lvpa"),
	 c("landuse.1.area",
	   "landuse.3.area",
	   "landuse.4.area",
	   "landuse.7.area")] <- 1

system.time(imp.tl.mice.val <- mice(tl.mice.val, seed=1234, meth="meth", pred="pred", minbucket=5))
print(imp.tl.mice.val)

				 COM.area = NA,
				 IND.area = NA,
				 MFR.area = NA,
				 SFR.area = NA,
				 COM.far = NA,
				 IND.far = NA,
				 MFR.far = NA,
				 SFR.far = NA,
				 COM.bvpa = NA,
				 IND.bvpa = NA,
				 MFR.bvpa = NA,
				 SFR.bvpa = NA,
				 COM.lvpa = NA,
				 IND.lvpa = NA,
				 MFR.lvpa = NA,
				 SFR.lvpa = NA)

### conduct a dry run to create the mids object called ini containing the default settings
ini <- mice(tl.mice.val, maxit=0, print=FALSE)
# create the matrices of some default settings
pred <- ini$predictorMatrix
meth <- ini$method

### redefine 'meth' for each variable and each interaction term
meth["landuse"] <- "cart"
meth["far"] <- "cart"
meth["bvpa"] <- "cart"
meth["lvpa"] <- "cart"
meth["COM.area"] <- "~I(landuse.1*area)"
meth["IND.area"] <- "~I(landuse.3*area)"
meth["MFR.area"] <- "~I(landuse.4*area)"
meth["SFR.area"] <- "~I(landuse.7*area)"
meth["COM.far"] <- "~I(landuse.1*far)"
meth["IND.far"] <- "~I(landuse.3*far)"
meth["MFR.far"] <- "~I(landuse.4*far)"
meth["SFR.far"] <- "~I(landuse.7*far)"
meth["COM.bvpa"] <- "~I(landuse.1*bvpa)"
meth["IND.bvpa"] <- "~I(landuse.3*bvpa)"
meth["MFR.bvpa"] <- "~I(landuse.4*bvpa)"
meth["SFR.bvpa"] <- "~I(landuse.7*bvpa)"
meth["COM.lvpa"] <- "~I(landuse.1*lvpa)"
meth["IND.lvpa"] <- "~I(landuse.3*lvpa)"
meth["MFR.lvpa"] <- "~I(landuse.4*lvpa)"
meth["SFR.lvpa"] <- "~I(landuse.7*lvpa)"

### redefine 'pred'

#test
pred[, c("area","landuse")] <- 0
pred["landuse", "area"] <- 1
pred[c("far","bvpa","lvpa"), c("COM.area",
			"IND.area",
			"MFR.area",
			"SFR.area")] <- 1



pred[, c("area","landuse","far","bvpa","lvpa")] <- 0
pred["landuse", c("area","far","bvpa","lvpa")] <- 1
pred["far", c("COM.area",
			"IND.area",
			"MFR.area",
			"SFR.area",
			"COM.bvpa",
			"IND.bvpa",
			"MFR.bvpa",
			"SFR.bvpa",
			"COM.lvpa",
			"IND.lvpa",
			"MFR.lvpa",
			"SFR.lvpa")] <- 1
pred["bvpa", c("COM.area",
			"IND.area",
			"MFR.area",
			"SFR.area",
			"COM.far",
			"IND.far",
			"MFR.far",
			"SFR.far",
			"COM.lvpa",
			"IND.lvpa",
			"MFR.lvpa",
			"SFR.lvpa")] <- 1
pred["lvpa", c("COM.area",
			"IND.area",
			"MFR.area",
			"SFR.area",
			"COM.far",
			"IND.far",
			"MFR.far",
			"SFR.far",
			"COM.bvpa",
			"IND.bvpa",
			"MFR.bvpa",
			"SFR.bvpa")] <- 1

### conduct MICE via CART in the validation data set
system.time(imp.tl.mice.val <- mice(tl.mice.val, seed=1234, method="meth", predictorMatrix="pred", minbucket=5))
print(imp.tl.mice.val)











# inspect if the landuse variable is MAR or MNAR
missing.landuse <- subset(tl, is.na(landuse))
missing.landuse <- missing.landuse[c("area","juris_city","bldgsqft","bldgval","landval","landuse")]
summary(missing.landuse)




### validation
# create a validation data set
# in which we draw a 5% random sample for each variable and replace them with NA
imp.val <- subset(imp, !is.na(far) & !is.na(bvpa) & !is.na(lvpa) & !is.na(landuse))
set.seed(1234)
imp.val$far[sample(nrow(imp.val), 0.05*nrow(imp.val))] <- NA
imp.val$bvpa[sample(nrow(imp.val), 0.05*nrow(imp.val))] <- NA
imp.val$lvpa[sample(nrow(imp.val), 0.05*nrow(imp.val))] <- NA
imp.val$landuse[sample(nrow(imp.val), 0.05*nrow(imp.val))] <- NA

# inspect missing values across the variables from the validation data set
md.pattern(imp.val)

# conduct MICE via CART using the validation data set
system.time(imp.val.mice <- mice(imp.val, seed=1234, method="cart", minbucket=5))
print(imp.val.mice)


return(far)
	return(bvpa)
	return(lvpa)
	return(landuse)
}
	far.rmse.1 <- sqrt(mean((far$'1' - far$far)^2))	
	far.rmse.2 <- sqrt(mean((far$'2' - far$far)^2))	
	far.rmse.3 <- sqrt(mean((far$'3' - far$far)^2))	
	far.rmse.4 <- sqrt(mean((far$'4' - far$far)^2))	
	far.rmse.5 <- sqrt(mean((far$'5' - far$far)^2))	
	#
	far.r2.1 <- summary(lm(far$'1' ~ $far))$r.squared
	far.r2.2 <- summary(lm(far$'2' ~ $far))$r.squared
	far.r2.3 <- summary(lm(far$'3' ~ $far))$r.squared
	far.r2.4 <- summary(lm(far$'4' ~ $far))$r.squared
	far.r2.5 <- summary(lm(far$'5' ~ $far))$r.squared
	#
	far.table <- matrix(
}

rmse.1 <- sqrt(mean((imp.far$'1'-imp.far$far)^2))
rmse.2 <- sqrt(mean((imp.far$'2'-imp.far$far)^2))
rmse.3 <- sqrt(mean((imp.far$'3'-imp.far$far)^2))
rmse.4 <- sqrt(mean((imp.far$'4'-imp.far$far)^2))
rmse.5 <- sqrt(mean((imp.far$'5'-imp.far$far)^2))
r2.1 <- summary(lm(imp.far$'1' ~ imp.far$far))$r.squared
r2.2 <- summary(lm(imp.far$'2' ~ imp.far$far))$r.squared
r2.3 <- summary(lm(imp.far$'3' ~ imp.far$far))$r.squared
r2.4 <- summary(lm(imp.far$'4' ~ imp.far$far))$r.squared
r2.5 <- summary(lm(imp.far$'5' ~ imp.far$far))$r.squared

table.far <- matrix(c(rmse.1, rmse.2, rmse.3, rmse.4, rmse.5, r2.1, r2.2, r2.3, r2.4, r2.5), 2, 5, byrow=TRUE)
dimnames(table.far) = list(c("rmse", "r2"), c("mi_1", "mi_2", "mi_3", "mi_4", "mi_5")) 
table.far


ds$juris_city <- as.factor(ds$juris_city)		# convert data type of juris_city from character to factor

ds$d_VAC <- with(ds, ifelse(landuse=='VAC', 'vacant', 'occupied'))
ds$d_bldgsqft <- with(ds, ifelse(bldgsqft==0, 'zero', 'non-zero'))
ds$d_bldgval <- with(ds, ifelse(bldgval==0, 'zero', 'non-zero'))
table(ds$d_VAC,ds$d_bldgsqft)
table(ds$d_VAC,ds$d_bldgval)


### create a data set for imputation
ds <- taxlots2011[c("area","bldgsqft","landval","bldgval","landuse","juris_city")]

SFR <- subset(ds, landuse=='SFR')				# select tax lots with SFR
SFR$bldgsqft[which(SFR$bldgsqft==0)] <- NA		# recode bldgsqft with zero to NA
SFR$landval[which(SFR$landval==0)] <- NA		# recode landval with zero to NA
SFR$bldgval[which(SFR$bldgval==0)] <- NA		# recode bldgval with zero to NA
SFR$juris_city <- as.factor(SFR$juris_city)		# convert data type of juris_city from character to factor

SFR$far <- SFR$bldgsqft / SFR$area				# floor area ratio
SFR$lvpa <- SFR$landval / SFR$area				# land value per area
SFR$bvpa <- SFR$bldgval / SFR$area				# building value per area

SFR <- SFR[c("area","juris_city","far","lvpa","bvpa")]

### validation 1

# create a validation data set
# in which we draw a 5% random sample for each variable and replace them with NA
SFR.val <- subset(SFR, !is.na(far) & !is.na(lvpa) & !is.na(bvpa))
set.seed(1234)
SFR.val$far[sample(nrow(SFR.val),0.05*nrow(SFR.val))] <- NA
SFR.val$lvpa[sample(nrow(SFR.val),0.05*nrow(SFR.val))] <- NA
SFR.val$bvpa[sample(nrow(SFR.val),0.05*nrow(SFR.val))] <- NA

# inspect missing values across the variables from the validation data set
md.pattern(SFR.val)

# conduct MICE via CART using the validation data set
system.time(imp.cart.SFR.val <- mice(SFR.val, seed=1234, method="cart", minbucket=5))
print(imp.cart.SFR.val)

# 

SFR$rowname <- as.numeric(row.names(SFR))



# validation 1
imp.far <- imp.cart.SFR1$imp$far
imp.far$rowname <- as.numeric(row.names(imp.far))
imp.far <- merge(imp.far, SFR, by="rowname")

imp.lvpa <- imp.cart.SFR1$imp$lvpa
imp.lvpa$rowname <- as.numeric(row.names(imp.lvpa))
imp.lvpa <- merge(imp.lvpa, SFR, by="rowname")

imp.bvpa <- imp.cart.SFR1$imp$bvpa
imp.bvpa$rowname <- as.numeric(row.names(imp.bvpa))
imp.bvpa <- merge(imp.bvpa, SFR, by="rowname")

# validation 2
imp.far <- imp.cart.SFR.val2$imp$far
imp.far$rowname <- as.numeric(row.names(imp.far))
imp.far <- merge(imp.far, SFR, by="rowname")

imp.lvpa <- imp.cart.SFR.val2$imp$lvpa
imp.lvpa$rowname <- as.numeric(row.names(imp.lvpa))
imp.lvpa <- merge(imp.lvpa, SFR, by="rowname")

imp.bvpa <- imp.cart.SFR.val2$imp$bvpa
imp.bvpa$rowname <- as.numeric(row.names(imp.bvpa))
imp.bvpa <- merge(imp.bvpa, SFR, by="rowname")

# validation 2
imp.far <- imp.cart.SFR.val3$imp$far
imp.far$rowname <- as.numeric(row.names(imp.far))
imp.far <- merge(imp.far, SFR, by="rowname")

imp.lvpa <- imp.cart.SFR.val3$imp$lvpa
imp.lvpa$rowname <- as.numeric(row.names(imp.lvpa))
imp.lvpa <- merge(imp.lvpa, SFR, by="rowname")

imp.bvpa <- imp.cart.SFR.val3$imp$bvpa
imp.bvpa$rowname <- as.numeric(row.names(imp.bvpa))
imp.bvpa <- merge(imp.bvpa, SFR, by="rowname")

rmse.1 <- sqrt(mean((imp.far$'1'-imp.far$far)^2))
rmse.2 <- sqrt(mean((imp.far$'2'-imp.far$far)^2))
rmse.3 <- sqrt(mean((imp.far$'3'-imp.far$far)^2))
rmse.4 <- sqrt(mean((imp.far$'4'-imp.far$far)^2))
rmse.5 <- sqrt(mean((imp.far$'5'-imp.far$far)^2))
r2.1 <- summary(lm(imp.far$'1' ~ imp.far$far))$r.squared
r2.2 <- summary(lm(imp.far$'2' ~ imp.far$far))$r.squared
r2.3 <- summary(lm(imp.far$'3' ~ imp.far$far))$r.squared
r2.4 <- summary(lm(imp.far$'4' ~ imp.far$far))$r.squared
r2.5 <- summary(lm(imp.far$'5' ~ imp.far$far))$r.squared

table.far <- matrix(c(rmse.1, rmse.2, rmse.3, rmse.4, rmse.5, r2.1, r2.2, r2.3, r2.4, r2.5), 2, 5, byrow=TRUE)
dimnames(table.far) = list(c("rmse", "r2"), c("mi_1", "mi_2", "mi_3", "mi_4", "mi_5")) 
table.far

rmse.1 <- sqrt(mean((imp.lvpa$'1'-imp.lvpa$lvpa)^2))
rmse.2 <- sqrt(mean((imp.lvpa$'2'-imp.lvpa$lvpa)^2))
rmse.3 <- sqrt(mean((imp.lvpa$'3'-imp.lvpa$lvpa)^2))
rmse.4 <- sqrt(mean((imp.lvpa$'4'-imp.lvpa$lvpa)^2))
rmse.5 <- sqrt(mean((imp.lvpa$'5'-imp.lvpa$lvpa)^2))
r2.1 <- summary(lm(imp.lvpa$'1' ~ imp.lvpa$lvpa))$r.squared
r2.2 <- summary(lm(imp.lvpa$'2' ~ imp.lvpa$lvpa))$r.squared
r2.3 <- summary(lm(imp.lvpa$'3' ~ imp.lvpa$lvpa))$r.squared
r2.4 <- summary(lm(imp.lvpa$'4' ~ imp.lvpa$lvpa))$r.squared
r2.5 <- summary(lm(imp.lvpa$'5' ~ imp.lvpa$lvpa))$r.squared

table.lvpa <- matrix(c(rmse.1, rmse.2, rmse.3, rmse.4, rmse.5, r2.1, r2.2, r2.3, r2.4, r2.5), 2, 5, byrow=TRUE)
dimnames(table.lvpa) = list(c("rmse", "r2"), c("mi_1", "mi_2", "mi_3", "mi_4", "mi_5")) 
table.lvpa

rmse.1 <- sqrt(mean((imp.bvpa$'1'-imp.bvpa$bvpa)^2))
rmse.2 <- sqrt(mean((imp.bvpa$'2'-imp.bvpa$bvpa)^2))
rmse.3 <- sqrt(mean((imp.bvpa$'3'-imp.bvpa$bvpa)^2))
rmse.4 <- sqrt(mean((imp.bvpa$'4'-imp.bvpa$bvpa)^2))
rmse.5 <- sqrt(mean((imp.bvpa$'5'-imp.bvpa$bvpa)^2))
r2.1 <- summary(lm(imp.bvpa$'1' ~ imp.bvpa$bvpa))$r.squared
r2.2 <- summary(lm(imp.bvpa$'2' ~ imp.bvpa$bvpa))$r.squared
r2.3 <- summary(lm(imp.bvpa$'3' ~ imp.bvpa$bvpa))$r.squared
r2.4 <- summary(lm(imp.bvpa$'4' ~ imp.bvpa$bvpa))$r.squared
r2.5 <- summary(lm(imp.bvpa$'5' ~ imp.bvpa$bvpa))$r.squared

table.bvpa <- matrix(c(rmse.1, rmse.2, rmse.3, rmse.4, rmse.5, r2.1, r2.2, r2.3, r2.4, r2.5), 2, 5, byrow=TRUE)
dimnames(table.bvpa) = list(c("rmse", "r2"), c("mi_1", "mi_2", "mi_3", "mi_4", "mi_5")) 
table.bvpa



# validation 2
# the first record of juris_city is replaced with NA
SFR.val2 <- SFR.val
SFR.val2[1,2] <- NA

md.pattern(SFR.val2)

system.time(imp.cart.SFR.val2 <- mice(SFR.val2, seed=1234, method="cart", minbucket=5))
print(imp.cart.SFR.val2)


# validation 3
# the first record of area is replaced with NA
SFR.val3 <- SFR.val
SFR.val3[1,1] <- NA

md.pattern(SFR.val3)

system.time(imp.cart.SFR.val3 <- mice(SFR.val3, seed=1234, method="cart", minbucket=5))
print(imp.cart.SFR.val3)











system.time(imp.pmm <- mice(SFR, seed=1234))



# stripplot
png("stripplot.pmm.png"); stripplot(imp.pmm, pch=20, cex=1.2); dev.off()
png("stripplot.cart.png"); stripplot(imp.cart, pch=20, cex=1.2); dev.off()

# xyplot
png("xyplot.pmm1.png"); xyplot(imp.pmm, bldgsqft~landval | .imp, pch=20, cex=1.4); dev.off()
png("xyplot.pmm2.png"); xyplot(imp.pmm, bldgsqft~bldgval | .imp, pch=20, cex=1.4); dev.off()
png("xyplot.pmm3.png"); xyplot(imp.pmm, bldgsqft~juris_city | .imp, pch=20, cex=1.4); dev.off()

png("xyplot.cart1.png"); xyplot(imp.cart, bldgsqft~landval | .imp, pch=20, cex=1.4); dev.off()
png("xyplot.cart2.png"); xyplot(imp.cart, bldgsqft~bldgval | .imp, pch=20, cex=1.4); dev.off()
png("xyplot.cart3.png"); xyplot(imp.cart, bldgsqft~juris_city | .imp, pch=20, cex=1.4); dev.off()

png("xyplot.cart1.png"); xyplot(imp.cart, far~lvpa | .imp, pch=20, cex=1.4); dev.off()
png("xyplot.cart2.png"); xyplot(imp.cart, far~bvpa | .imp, pch=20, cex=1.4); dev.off()
png("xyplot.cart3.png"); xyplot(imp.cart, far~juris_city | .imp, pch=20, cex=1.4); dev.off()

# densityplot
png("densityplot.pmm.bldgsqft.png"); densityplot(imp.pmm, ~bldgsqft | .imp); dev.off()
png("densityplot.pmm.landval.png"); densityplot(imp.pmm, ~landval | .imp); dev.off()
png("densityplot.pmm.bldgval.png"); densityplot(imp.pmm, ~bldgval | .imp); dev.off()

png("densityplot.cart.bldgsqft.png"); densityplot(imp.cart, ~bldgsqft | .imp); dev.off()
png("densityplot.cart.landval.png"); densityplot(imp.cart, ~landval | .imp); dev.off()
png("densityplot.cart.bldgval.png"); densityplot(imp.cart, ~bldgval | .imp); dev.off()

png("densityplot.cart.far.png"); densityplot(imp.cart, ~far | .imp); dev.off()
png("densityplot.cart.lvpa.png"); densityplot(imp.cart, ~lvpa | .imp); dev.off()
png("densityplot.cart.bvpa.png"); densityplot(imp.cart, ~bvpa | .imp); dev.off()


# after MICE
fit <- with(imp.cart, lm(far~lvpa+bvpa))
summary(fit)
est <- pool(fit)
est






