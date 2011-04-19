# +-----------------------------------------------------------------------------------------------+
# |									LONGITUDINAL DATA											  |
# +-----------------------------------------------------------------------------------------------+
# ####################################################
# ----------------------------------------------------
# 				Loading data
# ----------------------------------------------------
tbc <- read.table('data/tbc.csv', header=FALSE, sep=',')		# Loading comma-separated values
names(tbc) <- c('1996', '1997', '1998', 'Counts')				# Naming columns

# ----------------------------------------------------
# 	Generating primer data from the above crosstab
# ----------------------------------------------------
data <- tbc
r <- rep(data[,1], data[,4])								#
for (i in 2:3) {											# Transforming the crosstab to primer
	r <- cbind(r, rep(data[,i], data[,4]))					# data for all cases for later analyse
}															#
colnames(r) <- names(tbc)[1:3]								# Naming columns

# ----------------------------------------------------
# 			Generating Vebb diagram of counts
# ----------------------------------------------------
# How to install limma?
# See: http://www.ats.ucla.edu/stat/r/faq/venn.htm 
library(limma)								# Loading limma
c <- vennCounts(r)							# Counting the number of persons in each intersect
c[1,4] <- '?'								# Number of uncounted/unmarked homeless people
pdf('img/venn1.pdf')						# To save image
vennDiagram(c, counts.col='grey')			# Plot diagram
dev.off()									# Closing pdf device

# ----------------------------------------------------
# 			Descriptives of the CMR samples
# ----------------------------------------------------
library(Rcapture)							# Loading Rcapture
(desc <- descriptive(tbc, dfreq=TRUE))		# Printing descriptives
plot(desc)									# Plotting descriptives

# ----------------------------------------------------
# 	Buliding log-linear models of the CMR samples
# ----------------------------------------------------
closedp(tbc, dfreq=TRUE)							# All possible models
closedpCI.t(tbc, dfreq=TRUE, m='Mth', h='Chao')		# Compute confidence interval for given model

# +-----------------------------------------------------------------------------------------------+
# |									DATA of 1996												  |
# +-----------------------------------------------------------------------------------------------+
# ####################################################
# ----------------------------------------------------
# 				Loading data
# ----------------------------------------------------
d1996 <- read.table('1996.csv', header=FALSE, sep=',')

# ----------------------------------------------------
# 			Descriptives of the CMR samples
# ----------------------------------------------------
desc <- descriptive(d1996, dfreq=TRUE)				# Printing descriptives
plot(desc)											# Ploting descriptives

# ----------------------------------------------------
# 	Buliding log-linear models of the CMR samples
# ----------------------------------------------------
closedp.0(d1996, dfreq=TRUE)							# All possible models
closedpCI.0(d1996, dfreq=TRUE, m='M0')					# Compute confidence interval for given model
closedpCI.0(d1996, dfreq=TRUE, m='Mh', h='Chao')		#
closedpCI.0(d1996, dfreq=TRUE, m='Mh', h='Poisson')		#
closedpCI.0(d1996, dfreq=TRUE, m='Mh', h='Darroch')		#
closedpCI.0(d1996, dfreq=TRUE, m='Mh', h='Gamma')		#


# +-----------------------------------------------------------------------------------------------+
# |							Comparison of results												  |
# +-----------------------------------------------------------------------------------------------+
# ####################################################
# ----------------------------------------------------
# 		Get confidence intervals in a data frame
# ----------------------------------------------------
results <- data.frame(matrix(0, 4, 4))
names(results) <- c('method', 'mean', 'low', 'high')
results[1,] <- c('Dávid-Snijders (1996-1998)', c(12345, 8654, 17610))
results[2,] <- c('Mth Darroch (1996-1998)',
		as.numeric(closedpCI.t(tbc, dfreq=TRUE, m='Mth', h='Darroch')$CI))
results[3,] <- c('Dávid-Snijders (1996)', c(3913, 1605, 9545))
results[4,] <- c('Mh Chao (1996)', 
		as.numeric(closedpCI.0(d1996, dfreq=TRUE, m='Mh', h='Chao')$CI))
# ----------------------------------------------------
# 		Rearrange method's factor leveles
# ----------------------------------------------------
results$method <- factor(results$method,
		levels=c('Dávid-Snijders (1996-1998)',
				'Mth Darroch (1996-1998)',
				'Dávid-Snijders (1996)',
				'Mh Chao (1996)'))
results[,2] <- as.numeric(results[,2])
results[,3] <- as.numeric(results[,3])
results[,4] <- as.numeric(results[,4])
# ----------------------------------------------------
# 					Draw plot!
# ----------------------------------------------------
library(ggplot2)													# Loading ggplot2
p <- ggplot(results, aes(x=method, y=mean)) + geom_point() +		# Draw points of estimates
		geom_errorbar(aes(ymin=low, ymax=high)) + 					# with errorbars of given CI
		coord_flip() + 												# Flipping coordinates
		theme_bw() + 												# Setting theme
		xlab('') + 													# Whith no labels on x and
		ylab('')													# y axis
print(p)															# Show (or save) plot