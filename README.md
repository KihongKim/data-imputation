data-imputation
===============

Imputing Missing Land Use Data: A Multiple Imputation by Chained Equations (MICE) Approach Based on Recursive Partitioning

Kihong Kim and Liming Wang, Ph.D.

Nohad A. Toulan School of Urban Studies and Planning
Portland State University

Submitted for the 94th Annual Meeting of the Transportation Research Board on January 11-15, 2015

ABSTRACT

Missing data has always been a challenging problem in transportation and land use modeling and becomes more prominent as models move to use disaggregated and higher spatial resolution data, such as parcel data. So far missing land use data are largely handled in an ad-hoc basis without assessment of imputation quality. In this paper, we introduce a sophisticated multiple imputation process for missing values and test it with real data. Multiple imputation by chained equations (MICE) is a flexible and practical method to deal with multivariate missing data. The key of success in a MICE application depends on how well an imputation model is constructed for each missing variable. Recently, the use of recursive partitioning is suggested as imputation models to implement MICE because of the potential to capture complex interaction effects on missing variables with minimal effort to set up the models. In this study, we apply the non-parametric MICE approach to impute missing values in parcel data. We evaluate the performance of MICE with two recursive partitioning tools (i.e., CART and Random Forests) using a cross-validation technique. We find that, under limited computing power, the CART-based MICE performs better, especially for continuous missing variables that are severely skewed.
