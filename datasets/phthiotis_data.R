#####################################################
# Section 4.2: Phthiotis Population Dataset
#####################################################

data1 <- read.table('fthiotida.txt')
xnames <- data1[,1]
covs <- scale(cbind(data1[,23], log(data1[,26])))
xnames <- data1[,1]
# multinomial count data responses
y <- as.matrix(data1[,22:2])
# desing matrix
X <- as.matrix(cbind(1, covs))

head(y)

head(X)



