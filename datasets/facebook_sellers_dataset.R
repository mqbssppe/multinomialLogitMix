#############################################################
# Section 4.3: Facebook Live Sellers in Thailand Data Set
#############################################################

df <- read.csv("Live_20210128.csv")
df$status_published <- as.Date(df$status_published,format = "%m/%d/%Y")
ind <- which(df$status_published >= "2017-01-01")
df <- df[ind,]
ind <- which(df$status_published <= "2018-12-31")
df <- df[ind,]
ind <- which(df$num_reactions < 41)
df <- df[-ind,]
ind <- which(df$status_type == "link")
df <- df[-ind,]


set.seed(1)

ind1 <- sample(which(df$status_type == "video"),  100)
ind2 <- sample(which(df$status_type == "photo"),  100)                  
ind3 <- sample(which(df$status_type == "status"), 100)  
df <- df[c(ind1, ind2, ind3), ]


# multinomial count data (response)
y <- as.matrix(df[,12:7])
apply(y/rowSums(y), 2, quantile, probs = c(0.99))
x1 <- ifelse(df$status_type == 'status', 1, 0)
x2 <- ifelse(df$status_type == 'photo', 1, 0)
xType <- cbind(x1, x2)
x <- xType
# desing matrix
X <- as.matrix(cbind(1, xType))


head(y)

head(X)
