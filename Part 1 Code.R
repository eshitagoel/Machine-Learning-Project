########### ST443 - Final Project ############
####Dataset Link:https://www.kaggle.com/harlfoxem/housesalesprediction

setwd()
data = read.csv("kc_house_data.csv")

################################# Part 1 Data Preparation ################################
###1.1 Overview of the dataset
head(data)
summary(data)

###1.2 Missing values checking

table(is.na(data))

###1.3 Drop irrelevant columns (based on common sense)
data <- subset(data, select = c(-id,-zipcode))

###1.4 Change the format for variable 'date'
library(lubridate) 
data$date<-(substr(data$date, 1, 8))
data$date<- ymd(data$date)
data$date<-as.numeric(data$date)

###1.5 change numerical variables that majority of the observations equal to 0 to categorical variables
#most of houses are not renovated
data$yr_renovated[data$yr_renovated>0] <- 1 
#most of houses do not have basements 
data$sqft_basement[data$sqft_basement>0] <- 1 

###1.6 Round the values of specific variables
data$bathrooms=round(data$bathrooms)
data$floors<-round(data$floors)
#----------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------#
################################# Part 2 Exploratory Data Analysis ################################

### 2.1 Correlation between variables
library(reshape2)
library(corrplot)
cor(data)
corelation_matrix <- round(cor(data),2)
corrplot(corelation_matrix, method="shade",title = "Figure 1. Correlation between variables")

# pairs(data, cex=0.1)

### 2.2 Response variable-price
library(ggplot2)
summary(data$price)
boxplot(data$price,main="Figure 2. Boxplot for Response Variable - Housing Price")######log y linear####
hist(data$price, main="Figure 3. Histogram for Response Variable - Housing Price", xlab="Price")
hist(log(data$price), main="Figure 4. Histogram for Response Variable (log)- Housing Price", xlab="Price")

### 2.3 Relationship between Housing Prices and some significant numerical variables (correlation>0.2) 
numercial<-c(5,12,18)
par(mfrow=c(2,2))

for (i in numercial){
  plot(data[,i], data$price, main=paste("Scatterplot for",colnames(data)[i] ,"and housing price"),
  xlab=paste(colnames(data)[i]), ylab="price", pch=19)
  abline(lm(data$price~data[,i]), col="red") 
  lines(lowess(data[,i],data$price), col="blue")
}

### 2.4 Relationship between Housing Prices and some significant categorical variables (correlation>0.2) 
categorical<-c(3,4,7,8,9,11) 
par(mfrow=c(3,2))
for (i in categorical){
  boxplot(data[, 2] ~ data[, i], main=paste('Boxplot for',colnames(data)[i]),xlab=colnames(data)[i],ylab='price')
}

### 2.5 longitude and latitude map
library(magrittr)
ggplot(data,aes(x=long,y=lat,color=price)) + geom_point(alpha=0.02) + 
  ggtitle( "Spatial Layout of Housing Data colored by price" )

################################# Part 3 Feature Engineering ################################

####3.1 bedrooms
data[data$bedrooms >10,] #only two houses
data <- data[!data$bedrooms>10,]

###3.2 Convert data type to categorical
data$waterfront <-as.factor(data$waterfront)
data$floors <-as.factor(data$floors)
data$view <-as.factor(data$view)
data$condition <-as.factor(data$condition)
data$yr_renovated <-as.factor(data$yr_renovated)
data$sqft_basement <-as.factor(data$sqft_basement)

### 3.3 Other changes
#Convert building year into 3 categories
data$yr_built<-cut(data$yr_built,3)
data$yr_built<-factor(data$yr_built,levels = c("(1.9e+03,1.94e+03]","(1.94e+03,1.98e+03]","(1.98e+03,2.02e+03]"),labels = c("1900_1940","1940_1980","1980_2020"))
data$yr_built<-as.factor(data$yr_built)
par(mfrow=c(1,1))
boxplot(data$price ~ data$yr_built, main='Boxplot for built year',xlab='Built Year',ylab='Price')
#----------------------------------------------------------------------------------#


################################# Part 4 Regression Part ################################

#----------------------------------------------------------------------------------#
##We first need to split our data into training and validation sets so that we can test the accuracy of our models. 
set.seed(20)
train = sample(nrow(data), size = 0.8*nrow(data)) #20% validation data
train_set <- data[train,]
#----------------------------------------------------------------------------------#

##Part 4.1 linear Regression##
old <- Sys.time()
lm.fit = lm(price ~ ., data = data, subset = train)
summary(lm.fit)
lm.mse = mean((data$price - predict(object = lm.fit, newdata = data))[-train] ^ 2)
lm.mse #34613676333
(time.lm <- Sys.time() - old)
#----------------------------------------------------------------------------------#

##Part 4.2 log-linear Regression (as response variable is extremely right-skewed) 
old <- Sys.time()
loglinear_fit = lm(log(price) ~., data = data, subset = train)
summary(loglinear_fit)
# MSE OF testing data
ll.mse =mean((data$price - exp(predict(object = loglinear_fit, newdata = data)))[-train] ^ 2)
ll.mse #36074236094
(time.ll <- Sys.time() - old)
#----------------------------------------------------------------------------------#

## Part 4.3 Best Subset Selection 
old <- Sys.time()
library(leaps)
K =10
set.seed(20)
folds = sample(rep(1:10, length = nrow(train_set)))
predict_regsubsets = function(object, newdata, id, ...) {
  form = as.formula(object$call[[2]])
  mat = model.matrix(form, newdata)
  coef_i = coef(object, id = id)
  mat[, names(coef_i)] %*% coef_i
}

cv_errors =matrix(0, 10, 27)
for(k in 1:10){
  fit_bss = regsubsets(price ~ ., data=train_set[folds!=k,], nvmax=27)
  for(i in 1:27){
    pred =predict_regsubsets(fit_bss, train_set[folds==k,], id=i)
    cv_errors[k,i] =mean((train_set$price[folds==k]-pred)^2)
  }
}
mse_cv =apply(cv_errors,2,mean)
plot(mse_cv, ylab="MSE", xlab="Model Size", pch=19, type="b")
which.min(mse_cv)
points(which.min(mse_cv), mse_cv[which.min(mse_cv)], col="red", cex=2, pch=20)
bss_fit = regsubsets(price ~ ., data = train_set,nvmax = 27)
coef(bss_fit, id = 21)

pred =predict_regsubsets(bss_fit, data[-train,], id=21)
bss.mse =mean(((data[-train,]$price)-pred)^2)
bss.mse #34646016404
(time.bss <- Sys.time() - old)
#----------------------------------------------------------------------------------#

## Part 4.4 Ridge Regression 
old <- Sys.time()
library(glmnet)
x <-model.matrix(price~., data=train_set)[,-1]
y <-train_set$price
grid <- 10^seq(-2,10, length=100)
set.seed(20)
cv.ridge <-cv.glmnet(x, y, alpha=0,lambda=grid) #k=10 by default
plot(cv.ridge)
cv.ridge
coef(cv.ridge) # Coefficient vector corresponding to the mse which is within one standard error of the lowest mse using the best lambda.
test_mat <- model.matrix(price ~ ., data = data[-train,])
test_y <- data[-train,]$price
pred <- test_mat %*% coef(cv.ridge)
ridge.mse <- mean((test_y-pred)^2)
ridge.mse #34500445644
(time.ridge <- Sys.time() - old)
#----------------------------------------------------------------------------------#

## Part 4.5 Lasso Regression 
old <- Sys.time()
library('glmnet')
grid <- 10^seq(-2,10, length=100)
set.seed(20)
cv.lasso <-cv.glmnet(x, y,lamda=grid)
plot(cv.lasso)
cv.lasso
coef(cv.lasso) # coefficent vector corresponding to the mse which is within one standard error of the lowest mse using the best lambda.
x <-model.matrix(price~., data=train_set)[,-1]
y <-train_set$price
train_lasso <- subset(x, select = -c(sqft_lot,floors2,floors3,floors4,condition2,condition4,sqft_basement1,sqft_lot15))
train_lasso <- cbind(train_lasso,y)
colnames(train_lasso)[20] <- "price"
train_lasso <- as.data.frame(train_lasso)
lmlasso.fit = lm(price ~., data = train_lasso)
summary(lmlasso.fit)

test_mat = model.matrix(price ~ ., data = data[-train,])
pred = test_mat[, names(coef(lmlasso.fit))] %*% coef(lmlasso.fit)
test_y <- data[-train,]$price
lasso.mse <- mean((test_y-pred)^2)
lasso.mse #34591063624
(time.lasso <- Sys.time() - old)
#----------------------------------------------------------------------------------#

## Part 4.6 Principal Component Regression 
old <- Sys.time()
library(pls)
set.seed(20)
fit.pcr = pcr(price ~ ., data = data[train,], validation = "CV", scale = T)
summary(fit.pcr)
selectNcomp(fit.pcr, method = "onesigma", plot = TRUE)
fit.pcr = pcr(price ~ ., data = data[train,],scale = T,ncomp = 25)
bhat.pcr = fit.pcr$coefficients[,1,25]
sd = apply(x, 2, sd)
bhat.pcr.rec = bhat.pcr/sd
bhat.pcr.rec
bhat.lm = coef(lm(price ~ ., data = data[train,]))[-1]
cbind(bhat.pcr.rec, bhat.lm)

#see the results of test MSE of all the number of components 
mse = rep(NA, 27)
fit.pcr = pcr(price ~ ., data = data[train,], validation = "CV", scale = T)
for(j in 1:27){
  yhat = predict(fit.pcr, ncomp = j, newdata = data[-train,])
  mse[j] = mean((yhat - data$price[-train])^2)
}
plot(mse, type = "o", xlab = "Number of Components", ylab = "Test MSE" )
points(which.min(mse), mse[which.min(mse)], col="red", cex=2, pch=20)

#see the result of selected number of components(25)
yhat = predict(fit.pcr, ncomp = 25, newdata = data[-train,])
pcr.mse = mean((yhat - data$price[-train])^2)
pcr.mse #34744802347
(time.pcr <- Sys.time() - old)
#----------------------------------------------------------------------------------#

## Part 4.7 Generalized Additive Models (GAM)
old <- Sys.time()
library(gam)
gam.fit <-gam(price~ s(date) + s(bedrooms) + s(bathrooms) + s(sqft_living) + s(sqft_lot)
              + floors+ waterfront+ view + condition +s(grade)+s(sqft_above)+ sqft_basement + yr_built
              + yr_renovated +s(lat)+s(long)+s(sqft_living15)+s(sqft_lot15),data=data[train,]) 
summary(gam.fit) #all the variables are significant
par(mfrow=c(3,4))
plot(gam.fit, se=TRUE, col="blue",pages=1) #need to use 'arrow'to see the remaining plots
gam.mse <- mean((data$price - predict(object = gam.fit, newdata = data))[-train] ^ 2)
gam.mse #25632462458
(time.gam <- Sys.time() - old)
#----------------------------------------------------------------------------------#

## Part 4.8 Random Forest 
old <- Sys.time()
library(randomForest)
library(mlbench)
library(caret)
control <- trainControl(method="cv", number=5, search="grid") 
set.seed(20)
rfgrid <- expand.grid(.mtry=c(6:12))
#since random forest is not so sensitive to the number of trees we selected,
#for the sake of computational efficiency,we simply choose B=200 rather than tuning this hyperparameter
rf.fit <- train(price~., data=train_set, method="rf", metric='RMSE', importance=TRUE,tuneGrid=rfgrid, trControl=control,ntree=200)
print(rf.fit) #mtry =12
set.seed(20)
rf.fit <-randomForest(price~., data=train_set, mtry=12, importance=TRUE, ntree = 200)
plot(rf.fit)
importance(rf.fit)
varImpPlot(rf.fit) #Second: RSS
rf_pred = predict(rf.fit, data[-train,])
y_test = data[-train,][, 2]
rf.mse <- mean((y_test - rf_pred)^2)
rf.mse #14729682237
(time.rf <- Sys.time() - old)
#----------------------------------------------------------------------------------#

## Part 4.9 Gradient Boosting
old <- Sys.time()
library(gbm)

# try different combinations before, this gb_grid is the modified version
gb_grid = expand.grid(
  shrinkage = c(0.01, 0.05, 0.1), 
  interaction.depth = c(3,4,5), 
  n.minobsinnode = c(10,15),
  bag.fraction = c(0.8,0.9,1),
  n.trees = c(1000,2000,4000,5000), 
  MSE = 0)
for(i in 1:nrow(gb_grid)){
  set.seed(20)
  gb_tuning = gbm(
    price ~., 
    data = data[train,], 
    distribution = 'gaussian', 
    shrinkage = gb_grid$shrinkage[i],
    interaction.depth = gb_grid$interaction.depth[i],
    n.minobsinnode = gb_grid$n.minobsinnode[i],
    bag.fraction = gb_grid$bag.fraction[i],
    n.trees = gb_grid$n.trees[i],
    train.fraction=0.8)
  gb_grid$MSE[i] = gb_tuning$valid.error
}

gb_grid
which.min(gb_grid$MSE) #0.10 5 10 0.9 1000 

gb.fit <- gbm(price ~., data=data[train,],distribution = 'gaussian', shrinkage =0.1 ,
              interaction.depth = 5, n.minobsinnode = 10, bag.fraction = 0.9, n.trees =1000)
summary(gb.fit)
gb.pred <- predict(gb.fit,data[-train,])
y_test = data[-train,][, 2]
gb.mse <- mean((y_test-gb.pred)^2)
gb.mse #14515768968
(time.gb <- Sys.time() - old)
#----------------------------------------------------------------------------------#

## Part 4.10 XG Boosting
old <- Sys.time()
library(xgboost)
library(caret)
trControl = trainControl(method = "cv", number = 5, search = "grid")
set.seed(20)
xgbGrid <-  expand.grid(max_depth = c(3, 5, 7), 
                        nrounds = c(100,150,200,300,400,500),
                        subsample =c(0.5,0.8,1),
                        # below are default values
                        eta = 0.3,
                        gamma = 0,
                        min_child_weight = 1,
                        colsample_bytree = 0.6)
xgboosting.fit = train(price~., data = train_set, method = "xgbTree", trControl = trControl, tuneGrid = xgbGrid)
print(xgboosting.fit) # nrounds = 150, max_depth = 5, subsample = 1
xg_pred = predict(xgboosting.fit, data[-train,])
y_test = data[-train,][, 2]
xgb.mse= mean((y_test - xg_pred)^2)
xgb.mse  #15561512842
(time.xgb <- Sys.time() - old)
#----------------------------------------------------------------------------------#

## Part 4.11 Neutral Network
old <- Sys.time()
NNdata <- model.matrix(price~., data = data)
NNdata<-subset(NNdata, select = -c(1) )
NNdata<-cbind(NNdata,data$price)
colnames(NNdata)[28] <- 'price'

library(neuralnet)
set.seed(20)
maxs <- apply(NNdata, 2, max) 
mins <- apply(NNdata, 2, min)
scaled <- as.data.frame(scale(NNdata, center = mins, scale = maxs - mins))
NN.train <- scaled[train,]
NN.test<- scaled[-train,]
nn <- neuralnet(price~ date+bedrooms+bathrooms+sqft_living+sqft_lot+floors2
                +floors3+floors4+waterfront1+view1+view2+view3+view4+condition2
                +condition3+condition4+condition5+grade+sqft_above+sqft_basement1
                +yr_built1940_1980+yr_built1980_2020+yr_renovated1+lat+long+sqft_living15
                +sqft_lot15,data=NN.train,hidden=c(5,3),linear.output=T)
plot(nn)
pred.nn <- compute(nn,NN.test[,1:27])
pred.nn <- pred.nn$net.result*(max(data$price)-min(data$price))+min(data$price)
nnresult <- (NN.test$price)*(max(data$price)-min(data$price))+min(data$price)
NN.MSE<-mean((nnresult - pred.nn)^2)
NN.MSE #17482355982
(time.nn <- Sys.time() - old)
#----------------------------------------------------------------------------------#

################################# Part 5 Compare between different models ################################

model <- c('Linear Regression','Log-Linear','Best Subset Selection','Ridge Regression',
           'Lasso Regression','Principal Component Regression','Generalized Additive Models',
           'Random Forest','Gradient Boosting','XG Boosting','Neutral Network')
MSE <- c(lm.mse,ll.mse,bss.mse,ridge.mse,lasso.mse,pcr.mse,gam.mse,rf.mse,gb.mse,xgb.mse,NN.MSE)
Elapsed.Time <- c(time.lm,time.ll,time.bss,time.ridge,time.lasso,time.pcr,time.gam,time.rf,time.gb,time.xgb,time.nn)
models <- cbind(model,MSE,Elapsed.Time)
models
#----------------------------------------------------------------------------------#


################################# Part 6 Classification Analysis ###############################
#### 6.1 Convert price to categorical. 321725 is the 1st quantile; 645000 is the 3rd quantile
a=data$price >= 645000
b=data$price < 645000 & 321725 <= data$price
c=data$price < 321725
data$price[a] = "High"
data$price[b] = "Medium"
data$price[c] = "Low"
data$price <-as.factor(data$price)
attach(data)
#----------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------#
### 6.2 Multinomial Logistic Regression
old <- Sys.time()
library(nnet)

# Multinomial Logistic Regression
glm_fit = multinom(price ~.,data = data, subset = train)
# Predicted probabilities for the testing data set
glm_pred = predict(glm_fit, data[-train,], type = "class")
# Misclassfication error rate
error_rate_mlr =mean(glm_pred!=data$price[-train])
error_rate_mlr
(time_mlr <- Sys.time() - old)
#----------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------#
### 6.3 Linear Discriminant Analysis
old <- Sys.time()
library(MASS) 

# Linear Discriminant Analysis
lda_fit = lda(price ~., data=data, subset = train)
# Predicted probabilities for the testing data set
lda_pred = predict(lda_fit, data[-train,])$class
# Misclassfication error rate
error_rate_lda =mean(lda_pred!= data$price[-train])
error_rate_lda
(time_lda <- Sys.time() - old)
#----------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------#
### 6.4 Quadratic Discriminant Analysis
old <- Sys.time()
qda_fit = qda(price ~. , data=data, subset = train)
# Error in qda.default(x, grouping, ...) : rank deficiency in group Low
qda_pred = predict(qda_fit, data[-train,])$class
error_rate_qda =mean(qda_pred != data$price[-train])
(time_qda <- Sys.time() - old)
#----------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------#
### 6.5 k Nearest Neighbors
old <- Sys.time()
library(class)

# Create training data for X
train_X = cbind(data$date, data$bedrooms, data$bathrooms, data$sqft_living, data$sqft_lot, data$floors, 
                data$waterfront, data$view, data$condition, data$grade, data$sqft_above, data$sqft_basement, 
                data$yr_built, data$yr_renovated, data$lat, data$long, data$sqft_living15, data$sqft_lot15)[train,]
# Create testing data for X
test_X = cbind(data$date, data$bedrooms, data$bathrooms, data$sqft_living, data$sqft_lot, data$floors, 
               data$waterfront, data$view, data$condition, data$grade, data$sqft_above, data$sqft_basement, 
               data$yr_built, data$yr_renovated, data$lat, data$long, data$sqft_living15, data$sqft_lot15)[-train,]
# Create training data for Y
train_level = data$price[train]
# Create testing data for Y
level_test = data[-train,]$price

#### 6.5.1 Find the optimal k
erro =matrix(0, 200 , 1)
for(k in 1:200){
  knn_pred = knn(train_X, test_X, train_level, k)
  erro[k]=mean(knn_pred!=level_test)
}
plot(erro, xlab="k", ylab="Misclassfication error rate", type="l")
(best_model = which.min(erro)) 
points(best_model, erro[best_model], col="red", cex=2, pch=20)  

## optimal k is 39 (odd number)
knn_pred = knn(train_X, test_X, train_level, k=39)
error_rate_knn = mean(knn_pred!=level_test)
error_rate_knn
(time_knn <- Sys.time() - old)
#----------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------#
### 6.6 Classification Trees
old <- Sys.time()
library(tree)

price.test <- data[-train, "price"]
tree.price <- tree(price ~ ., data , subset=train)
plot(tree.price)
text(tree.price, pretty = 0)

## pruning the tree
cv.price <-cv.tree(tree.price, FUN = prune.misclass)
plot(cv.price$size, cv.price$dev, type="b")
# prune the tree to get the 4-node tree
prune.price <-prune.misclass(tree.price, best=4)
plot(prune.price)
text(prune.price, pretty = 0)

tree.pred <-predict(prune.price, price.test, type="class")[-train]
table(tree.pred, price.test)
error_rate_tree = mean(tree.pred!=price.test)
error_rate_tree
(time_tree <- Sys.time() - old)
#----------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------#
#### 6.6.1 Random Forest
old <- Sys.time()
library(randomForest)

rf.price <- randomForest(price~., data=data, subset=train, mtry=sqrt(18), importance=TRUE)
varImpPlot(rf.price)
plot(rf.price)

yhat.rd <-predict(rf.price, newdata=data[-train,])
table(yhat.rd, price.test)
error_rate_rm = mean(yhat.rd!=price.test)
error_rate_rm
(time_rm <- Sys.time() - old)
#----------------------------------------------------------------------------------#


############# Part 7 Compare between different Classification models ################
error_rate=c(error_rate_mlr, error_rate_lda, error_rate_knn, error_rate_tree, error_rate_rm)
model <- c('Multinomial Logistic Regression','Linear Discriminant Analysis','k Nearest Neighbors',
           'Pruned Classification Tree','Random Forest')
Elapsed_Time <- c(time_mlr,time_lda,time_knn,time_tree,time_rm)
(models <- cbind(model,error_rate,Elapsed_Time))
plot(error_rate, xlab="different Classification models", ylab="Misclassfication error rate", type="l")
best_model = which.min(error_rate)
points(best_model, error_rate[best_model], col="red", cex=2, pch=20)
#----------------------------------------------------------------------------------#