# Impute Missing Values on Parcel Data
# by MICE (Multiple Imputation by Chained Equations) with Recursive Partitioning
################################### 
# KiHong Kim and Liming Wang, Ph.D.
# Portland State University
###################################
# the relevant paper was submitted for the 2015 TRB meeting
 
setwd("/home/kihong/Imputation") # on the sapporo server

library(mice)
library(CALIBERrfimpute)
library(missForest)
library(RPostgreSQL)
library(lattice)

### create a data set to conduct MICE
# function
createDS <- function(taxlots) {
	tl <- taxlots
	# define missingness
	tl$bldgsqft[which(tl$landuse %in% c("COM","IND","MFR","SFR") & tl$bldgsqft==0)] <- NA
	tl$bldgval[which(tl$landuse %in% c('COM','IND','MFR','SFR') & tl$bldgval==0)] <- NA
	tl$landval[which(tl$landuse %in% c('COM','IND','MFR','SFR') & tl$landval==0)] <- NA
	# transform missing variables
	tl$far <- tl$bldgsqft / tl$area		# building SQFT per parcel area
	tl$bvpa <- tl$bldgval / tl$area		# building value per parcel area
	tl$lvpa <- tl$landval / tl$area		# land value per parcel area
	# select parcels in Multnohma County
	mult <- subset(tl,county=="051")
	# select independent variables used for an MICE imputer
	# non-missing variables: area, district, pop10, job11
	# missing variables: landuse, far, bvpa, lvpa
	vars <- c("area","district","pop10","job11","landuse","far","bvpa","lvpa")
	# create a data set for MICE
	ds <- mult[vars]
	ds$district <- factor(ds$district) # from numeric to factor
	ds$landuse <- factor(ds$landuse) # from numeric to factor
	return(ds)
}
# connect to the PostgreSQL database where taxlot tables are stored
conn <- dbConnect(PostgreSQL(), host="sapporo.usp.pdx.edu", user="smartdata", password="Smartaa00", dbname="portland")
taxlots2011 <- dbReadTable(conn, c("rlis","taxlots2011"))
# create a data set to perform MICE
createDS(taxlots2011)
#
save(ds, file="ds.RData")

### inspect missing pattern
md.pattern(ds)

### conduct the CART-based MICE
# according to the validation results shown below,
# the CART-based MICE performed better than MICE with PMM and Random Forests
imp.ds.cart <- mice(ds, seed=1234, method="cart")
#
print(imp.ds.cart)
save(imp.ds.cart, file="imp.ds.cart.RData")

### graphically diagnose the results of CART-based MICE
png("stripplot.far.png")
stripplot(imp.mult.mice.cart, far~.imp, ylab="FAR", pch=20, cex=1.2)
dev.off()
png("stripplot.bvpa.png")
stripplot(imp.mult.mice.cart, bvpa~.imp, ylab="BVPA", pch=20, cex=1.2)
dev.off()
png("stripplot.lvpa.png")
stripplot(imp.mult.mice.cart, lvpa~.imp, ylab="LVPA", pch=20, cex=1.2)
dev.off()

##############
# validation #
##############

### create a data set for validation
# function
createDS.val <- function(ds) {
	set.seed(1234)
	ds.val <- subset(ds, !is.na(landuse) & !is.na(far) & !is.na(bvpa) & !is.na(lvpa))
	ds.val.missVars <- ds.val[c("landuse","far","bvpa","lvpa")]
	ds.val.missVars <- prodNA(ds.val.missVars, noNA=0.05) # generate Missing At Random (MAR) using the missForest package
	ds.val <- cbind(ds.val[,c(1,2,3,4)], ds.val.missVars)
	return(ds.val)
}
#
createDS.val(ds)
#
save(ds.val, file="ds.val.RData")

### inspect missing pattern
md.pattern(ds.val)

### conduct MICE via PMM, CART, and Random Forests
# load("ds.val.RData")
#
imp.ds.val.pmm <- mice(ds.val, seed=1234, method="pmm")
print(imp.ds.val.pmm)
save(imp.ds.val.pmm, file="imp.ds.val.pmm.RData")
#
imp.ds.val.cart <- mice(ds.val, seed=1234, method="cart")
print(imp.ds.val.cart)
save(imp.ds.val.cart, file="imp.ds.val.cart.RData")
#
imp.ds.val.rf <- mice(ds.val, seed=1234, method=method=c("","","","","rfcat","rfcont","rfcont","rfcont"))
print(imp.ds.val.rf)
save(imp.ds.val.rf, file="imp.ds.val.rf.RData")

### validations
# validation function for a continuous missing variable
validateMICE <- function(x,y) {
	# x: a imputed data set
	# y: a missing variable for validation
	ds$rowname <- as.numeric(row.names(ds))
	imputed <- x[["imp"]][[y]]
	imputed$rowname <- as.numeric(row.names(imputed))
	merged <- merge(imputed, ds[c(y,"rowname")], by="rowname")
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
# validation function for a categorical missing variable
validateMICE2 <- function(x,y) {
	# x: a imputed data set
	# y: a missing variable for validation
	ds$rowname <- as.numeric(row.names(ds))
	imputed <- x[["imp"]][[y]]
	imputed$rowname <- as.numeric(row.names(imputed))
	merged <- merge(imputed, ds[c(y,"rowname")], by="rowname")
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
# validation: PMM
validateMICE(imp.ds.val.pmm, "far")
validateMICE(imp.ds.val.pmm, "bvpa")
validateMICE(imp.ds.val.pmm, "lvpa")
validateMICE2(imp.ds.val.pmm, "landuse")
# validation: cart
validateMICE(imp.ds.val.cart, "far")
validateMICE(imp.ds.val.cart, "bvpa")
validateMICE(imp.ds.val.cart, "lvpa")
validateMICE2(imp.ds.val.cart, "landuse")
# validation: rf
validateMICE(imp.ds.val.rf, "far")
validateMICE(imp.ds.val.rf, "bvpa")
validateMICE(imp.ds.val.rf, "lvpa")
validateMICE2(imp.ds.val.rf, "landuse")

########################################
# regression analysis after imputation #
########################################



