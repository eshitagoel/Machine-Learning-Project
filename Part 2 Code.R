## PART 2

##0.1 We first need to find the matrix B 

p <- 100
set.seed(7)
B <- matrix(sample(c(0,0.5), p^2, replace=TRUE, prob=c(0.9,0.1)), p)
i <- lower.tri(B)
B[i] <- t(B)[i]
diag(B) <- rep(0, p)
B

#0.2 To calculate delta, we apply the definition of a matrix being positive definite

delta <- -min(eigen(B, symmetric=TRUE, only.values=TRUE)$values)
delta 

# We can choose any value of delta that is more than this, so let us put delta = 4

delta <- 4

# We now calculate theta

theta = B + delta*diag(p)
theta

#0.3 We check if the matrix is positive definite: 

# install.packages("matrixcalc")
library(matrixcalc)

is.positive.definite(theta)

#0.4 Now we convert the covariance matrix to the correlation matrix using the cov2cov() function

theta <- cov2cor(theta)
theta

#0.5 The Edge set are all pairs (i,j) such that the ij^th entry on the matrix is non zero

df = data.frame(paste(rep(1:p,each =p),1:p,sep=","),rep(0,p*p))
names(df)<- c("Nodes","Edges")
df

for (i in 1:p){
  for (j in 1:p){
    if(theta[i,j]!=0){
      df["Edges"][df['Nodes'] == paste(i,j,sep=",")] = 1
    }
  }
}
df

node.lasso <- data.frame(df)
graph.lasso <- data.frame(df)


##1 We now generate our samples x1, ..., xn from multivariate normal distribution 
# with zero mean and covariance matrix sigma = inverse of theta

sigma = solve(theta)
set.seed(7)

library(MASS)
X <- mvrnorm(n=200, mu=rep(0,100), Sigma = sigma)

##2.1 Node-Wise Lasso Approach

#MSE Measurement (not the best measurement)
library(glmnet)

node.lasso['Edges'] = rep(0,p)
grid <- 10^seq(-6,-2, length=100) 
K =5
set.seed(0)
folds = sample(rep(1:5, length = 200))
cv_error = matrix(0,100,80)


for (i in 1:p){
  
  if(i==1){
    x <- X[,2:p]
    y <- X[,1]
  } else if(i==p){
    x <- X[,1:p-1]
    y <- X[,p]
  }else{
    x <- X[,c(1:(i-1),(i+1):p)]
    y <- X[,i]
  }
  for (g in 1:80){
    for(k in 1:5){
      error <- rep(0,5)
      lasso = glmnet(x,y, data=X[folds!=k,],alpha=1,lambda=grid[g])
      pred =predict(lasso, x[folds==k,])
      error[k]=(y[folds==k]-pred)^2
      cv_error[i,g] = mean(error)
    }
  }
}

mse_cv =apply(cv_error,2,sum)
plot(mse_cv, ylab="MSE", xlab="lambda", pch=19, type="b")
which.min(mse_cv)
points(which.min(mse_cv), mse_cv[which.min(mse_cv)], col="red", cex=2, pch=20)


###True Positive Rate Measurement
library(MASS)
set.seed(7)
X <- mvrnorm(n=200, mu=rep(0,100), Sigma = sigma)
node.lasso['Edges'] = rep(0,p)
grid <- 10^seq(-2,2, length=10) 
TP_rate = matrix(0,p,10)
NP = sum(df$Edges == 1) 

for (g in 1:10){
  for (i in 1:p){
    
    if(i==1){
      x <- X[,2:p]
      y <- X[,1]
    } else if(i==p){
      x <- X[,1:p-1]
      y <- X[,p]
    }else{
      x <- X[,c(1:(i-1),(i+1):p)]
      y <- X[,i]
    }
    for (j in 2:p){
      if (i==1){
        if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
      } else if(i==p){
        if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
      } else{
        if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
          if(j<=i){
            node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
          } else{
            node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
          }
        }
      }
      
      TP_rate[i,g] = sum(df$Edges == 1 & node.lasso$Edges == 1) / NP 
    }
  }
}

TP_lambda =apply(TP_rate,2,sum)
plot(TP_lambda, ylab="sun of TP", xlab="lambda", pch=19, type="b")
points(which.max(TP_lambda), TP_lambda[which.max(TP_lambda)], col="red", cex=2, pch=20)
lambda.max =grid[which.max(TP_lambda)] #0.2154435

## Apply the optimal lambda

node.lasso <- data.frame(df)
node.lasso['Edges'] = rep(0,p)

for (i in 1:p){
  
  if(i==1){
    x <- X[,2:p]
    y <- X[,1]
  } else if(i==p){
    x <- X[,1:p-1]
    y <- X[,p]
  }else{
    x <- X[,c(1:(i-1),(i+1):p)]
    y <- X[,i]
  }
  
  for (j in 2:p){
    if (i==1){
      if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
        node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
    } else if(i==p){
      if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
        node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
    } else{
      if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
        if(j<=i){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
        } else{
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
        }
      }
    }
  }
}

#2.1.1 Implementing "Joint" and "Or" Rule to retrive edges

node.lasso1 <- data.frame(df)
node.lasso2 <- data.frame(df)

node.lasso1['Edges'] = rep(0,p)
node.lasso2['Edges'] = rep(0,p)

for (i in 1:p){
  for (j in 1:i){
    if ((node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] == 1) & (node.lasso["Edges"][node.lasso['Nodes'] == paste(j,i,sep=",")] == 1) ){
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 1
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 1
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 1
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 1
    } else if((node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] == 1) | (node.lasso["Edges"][node.lasso['Nodes'] == paste(j,i,sep=",")] == 1) ){
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 1
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 1
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 0
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 0
    } else{
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 0
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 0
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 0
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 0
    }
  }
}

#2.1.2 Checking if the edges are symmetric 

for (i in 1:p){
  for (j in 1:p){
    if(node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] != node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")]){
      print('Not Symmetric')
      break
    }
    
    if(node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] != node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")]){
      print('Not Symmetric')
      break
    }
  }
}

# Proportion of correct classification

nwlasso1_mean_200 = mean(df['Edges'] == node.lasso1['Edges'])
nwlasso2_mean_200 = mean(df['Edges'] == node.lasso2['Edges'])

##2.2 Graphical Lasso Approach

# create the covariance matrix for X
s <- var(X)
s

# Using glasso

library(glasso)
g.lasso <- glasso(s, rho=0.1)

summary(g.lasso)

# We want the estimated inverse covariance matrix : 

omega <- g.lasso$wi

graph.lasso['Edges'] = rep(0,p)

for (i in 1:p){
  for (j in 1:p){
    if(omega[i,j]!=0){
      graph.lasso["Edges"][graph.lasso['Nodes'] == paste(i,j,sep=",")] = 1
    }
  }
}

# Proportion of correct classification

mean(df['Edges'] == graph.lasso['Edges'])


##3
NP = sum(df$Edges == 1) # No. of Positive
NN = sum(df$Edges == 0) # No. of Negative

#3.1 Node-Wise Lasso Approach: best values for tuning parameteres

(TP1.1 = sum(df$Edges == 1 & node.lasso1$Edges == 1) / NP) # True Positive Rate
(TP1.2 = sum(df$Edges == 1 & node.lasso2$Edges == 1) / NP) # True Positive Rate

(FP1.1 = sum(df$Edges == 0 & node.lasso1$Edges == 1) / NN) # False Positive Rate
(FP1.2 = sum(df$Edges == 0 & node.lasso2$Edges == 1) / NN) # False Positive Rate

#3.2 Graphical Lasso Approach: Applying Cross Validation to get the best values for tuning parameteres

library(CVglasso)
cv.glasso <-CVglasso(X)
plot(cv.glasso)

cv.glasso$Tuning
cv.omega = cv.glasso$Omega

graph.lasso['Edges'] = rep(0,p)

for (i in 1:p){
  for (j in 1:p){
    if(cv.omega[i,j]!=0){
      graph.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
    }
  }
}

# Proportion of correct classification

Glasso_mean_200 = mean(df['Edges'] == graph.lasso['Edges'])

(TP2 = sum(df$Edges == 1 & graph.lasso$Edges == 1) / NP) # True Positive Rate
(FP2 = sum(df$Edges == 0 & graph.lasso$Edges == 1) / NN) # False Positive Rate

#3.3 Plotting TPRλ vs FPRλ over a fine grid of values of λ produces a ROC curve.

#3.3.1 Node - wise Lasso ROC curve
roc.lasso <- data.frame(df)

value <- 10^seq(-2,0,length.out=10)
threshold = sort(unique(c(0, value)))
nth = length(threshold)
FP1 = TP1 = rep(0, nth)

for (k in 1:nth){
  roc.lasso['Edges'] = rep(0,p)
  
  for (i in 1:p){
    if(i==1){
      x <- X[,2:p]
      y <- X[,1]
    } else if(i==p){
      x <- X[,1:p-1]
      y <- X[,p]
    }else{
      x <- X[,c(1:(i-1),(i+1):p)]
      y <- X[,i]
    }
    fit.lasso <-glmnet(x,y, lambda = threshold[k])
    for (j in 2:p){
      if (i==1){
        if (coef(fit.lasso)[j]!=0){
          roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
      } else if(i==p){
        if (coef(fit.lasso)[j]!=0){
          roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
      } else{
        if (coef(fit.lasso)[j]!=0){
          if(j<=i){
            roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
          } else{
            roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j,sep=",")] = 1
          }
        }
      }
    }
  }
  
  TP1[k] = sum(df$Edges == 1 & roc.lasso$Edges == 1) / NP # True Positive Rate
  FP1[k] = sum(df$Edges == 0 & roc.lasso$Edges == 1) / NN # False Positive Rate
}

plot(
  x = c(0,1),
  y = c(0,1),
  type = "n",
  main = "ROC",
  xlab = "False Positive Rate",
  ylab = "True Positive Rate"
)
lines(FP1, TP1, col = "black", lwd = 3)

#3.4.1 Calculating Area Under Curve

library(pROC)
auc.lasso = auc(df$Edges, roc.lasso$Edges)

#3.3.2 Graphical Lasso ROC curve

roc.graph <- data.frame(df)

threshold = sort(unique(c(0, cv.glasso$Lambdas)))
nth = length(threshold)
FP.g.200 = TP.g.200 = rep(0, nth)

for (k in 1:nth){
  roc.graph['Edges'] = rep(0,p)
  fit.graph <- glasso(s, rho= threshold[k])
  
  # We want the estimated inverse covariance matrix : 
  
  omega <- fit.graph$wi
  
  for (i in 1:p){
    for (j in 1:p){
      if(omega[i,j]!=0){
        roc.graph["Edges"][roc.graph['Nodes'] == paste(i,j,sep=",")] = 1
      }
    }
  }
  
  
  TP.g.200[k] = sum(df$Edges == 1 & roc.graph$Edges == 1) / NP # True Positive Rate
  FP.g.200[k] = sum(df$Edges == 0 & roc.graph$Edges == 1) / NN # False Positive Rate
}

plot(
  x = c(0,1),
  y = c(0,1),
  type = "n",
  main = "ROC",
  xlab = "False Positive Rate",
  ylab = "True Positive Rate"
)
lines(FP.g.200, TP.g.200, col = "black", lwd = 3)

#3.4.2 Calculating Area Under Curve

library(pROC)
auc.graph = auc(df$Edges, roc.graph$Edges)

# Using roc() function to plot ROC Curve
plot(roc(df$Edges,roc.graph$Edges))


##4 compare the sample performance: False Positives, False Negatives.

#4.1 Node - wise Lasso

FP1.1
FP1.2
FN1.1 = sum(df$Edges == 1 & node.lasso1$Edges == 0) / NP # False Negatives Rate
FN1.2 = sum(df$Edges == 1 & node.lasso2$Edges == 0) / NP # False Negatives Rate
FN1.1
FN1.2

#4.2 Graphical Lasso

FP2
FN2 = sum(df$Edges == 1 & graph.lasso$Edges == 0) / NP # False Negatives Rate
FN2


#5 We procedure 50 times

##5.1  Replicate the Node-Wise Lasso Approach 50 times

mean.lasso1.1_50 = rep(0,50)
mean.lasso1.2_50 = rep(0,50)
FP1.1_50 = rep(0,50)
FP1.2_50 = rep(0,50)
FN1.1_50 = rep(0,50)
FN1.2_50 = rep(0,50)

for (n in 1:50){
  X <- mvrnorm(n=200, mu=rep(0,100), Sigma = sigma)

  node.lasso['Edges'] = rep(0,p)
  grid <- 10^seq(-2,2, length=10) 
  TP_rate = matrix(0,p,10)
  NP = sum(df$Edges == 1) 
  
  for (g in 1:10){
    for (i in 1:p){
      
      if(i==1){
        x <- X[,2:p]
        y <- X[,1]
      } else if(i==p){
        x <- X[,1:p-1]
        y <- X[,p]
      }else{
        x <- X[,c(1:(i-1),(i+1):p)]
        y <- X[,i]
      }
      for (j in 2:p){
        if (i==1){
          if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
            node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
        } else if(i==p){
          if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
            node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
        } else{
          if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
            if(j<=i){
              node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
            } else{
              node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
            }
          }
        }
        
        TP_rate[i,g] = sum(df$Edges == 1 & node.lasso$Edges == 1) / NP 
      }
    }
  }
  
  TP_lambda =apply(TP_rate,2,sum)
  plot(TP_lambda, ylab="sun of TP", xlab="lambda", pch=19, type="b")
  points(which.max(TP_lambda), TP_lambda[which.max(TP_lambda)], col="red", cex=2, pch=20)
  lambda.max =grid[which.max(TP_lambda)] 
  
  ## Apply the optimal lambda
  
  node.lasso <- data.frame(df)
  node.lasso['Edges'] = rep(0,p)
  
  for (i in 1:p){
    
    if(i==1){
      x <- X[,2:p]
      y <- X[,1]
    } else if(i==p){
      x <- X[,1:p-1]
      y <- X[,p]
    }else{
      x <- X[,c(1:(i-1),(i+1):p)]
      y <- X[,i]
    }
    
    for (j in 2:p){
      if (i==1){
        if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
      } else if(i==p){
        if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
      } else{
        if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
          if(j<=i){
            node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
          } else{
            node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
          }
        }
      }
    }
  }
  
  #Implementing "Joint" and "Or" Rule to retrive edges
  
  node.lasso1 <- data.frame(df)
  node.lasso2 <- data.frame(df)
  
  node.lasso1['Edges'] = rep(0,p)
  node.lasso2['Edges'] = rep(0,p)
  
  for (i in 1:p){
    for (j in 1:i){
      if ((node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] == 1) & (node.lasso["Edges"][node.lasso['Nodes'] == paste(j,i,sep=",")] == 1) ){
        node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 1
        node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 1
        node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 1
        node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 1
      } else if((node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] == 1) | (node.lasso["Edges"][node.lasso['Nodes'] == paste(j,i,sep=",")] == 1) ){
        node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 1
        node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 1
        node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 0
        node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 0
      } else{
        node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 0
        node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 0
        node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 0
        node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 0
      }
    }
  }
  mean.lasso1.1_50[n] = mean(df['Edges'] == node.lasso1['Edges'])
  mean.lasso1.2_50[n] = mean(df['Edges'] == node.lasso2['Edges'])
  FP1.1_50[n] = sum(df$Edges == 0 & node.lasso1$Edges == 1) / NN # False Positive Rate
  FP1.2_50[n] = sum(df$Edges == 0 & node.lasso2$Edges == 1) / NN # False Positive Rate
  FN1.1_50[n] = sum(df$Edges == 1 & node.lasso1$Edges == 0) / NP # False Negatives Rate
  FN1.2_50[n] = sum(df$Edges == 1 & node.lasso2$Edges == 0) / NP # False Negatives Rate
  }

NWlasso_mean1.1_50 = mean(mean.lasso1.1_50)
NWlasso_mean1.2_50 = mean(mean.lasso1.2_50)

# Mean of False Positives and False Negatives rate
mean.FP1.1_50 = mean(FP1.1_50)
mean.FP1.2_50 = mean(FP1.2_50)
mean.FN1.1_50 = mean(FN1.1_50)
mean.FN1.2_50 = mean(FN1.2_50)

# standard error of False Positives and False Negatives rate
sd.FP1.1_50 = sd(FP1.1_50)
sd.FP1.2_50 = sd(FP1.2_50)
sd.FN1.1_50 = sd(FN1.1_50)
sd.FN1.2_50 = sd(FN1.2_50)

x = c(FP1.1_50, FP1.2_50)
gp = c(rep(1,50), rep(2,50))
boxplot(x ~ gp, col=c("blue","red"), pch=19)

x = c(FN1.1_50, FN1.2_50)
gp = c(rep(1,50), rep(2,50))
boxplot(x ~ gp, col=c("blue","red"), pch=19)


#5.2 Replicate the Graphical Lasso Approach 50 times
mean.lasso2_50 = rep(0,50)
FP2_50 = rep(0,50)
FN2_50 = rep(0,50)

# Replicate 50 times
for (n in 1:50){
  X <- mvrnorm(n=200, mu=rep(0,100), Sigma = sigma)

# Cross Validation
cv.glasso <-CVglasso(X)
cv.omega = cv.glasso$Omega
graph.lasso['Edges'] = rep(0,p)

for (i in 1:p){
  for (j in 1:p){
    if(cv.omega[i,j]!=0){
      graph.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
    }
  }
}

# Proportion of correct classification
mean.lasso2_50[n]=mean(df['Edges'] == graph.lasso['Edges'])

# False Positives, False Negatives
FP2_50[n] = sum(df$Edges == 0 & graph.lasso$Edges == 1) / NN # False Positive Rate
FN2_50[n] = sum(df$Edges == 1 & graph.lasso$Edges == 0) / NP # False Negatives Rate

}
Glasso_mean = mean(mean.lasso2_50)

# Mean of False Positives and False Negatives rate
mean.FP2_50 = mean(FP2_50)
mean.FN2_50 = mean(FN2_50)

# standard error of False Positives and False Negatives rate
sd.FP2_50 = sd(FP2_50)
sd.FN2_50 = sd(FN2_50)

boxplot(FP2_50)
boxplot(FN2_50)


##6 various simulation settings for n and p
#6.1 n<p (p=100, n=50)
#6.1.1 Node-Wise Lasso Approach
set.seed(7)
X <- mvrnorm(n=50, mu=rep(0,100), Sigma = sigma)
node.lasso['Edges'] = rep(0,p)
grid <- 10^seq(-2,2, length=10) 
TP_rate = matrix(0,p,10)
NP = sum(df$Edges == 1) 

for (g in 1:10){
  for (i in 1:p){
    
    if(i==1){
      x <- X[,2:p]
      y <- X[,1]
    } else if(i==p){
      x <- X[,1:p-1]
      y <- X[,p]
    }else{
      x <- X[,c(1:(i-1),(i+1):p)]
      y <- X[,i]
    }
    for (j in 2:p){
      if (i==1){
        if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
      } else if(i==p){
        if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
      } else{
        if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
          if(j<=i){
            node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
          } else{
            node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
          }
        }
      }
      
      TP_rate[i,g] = sum(df$Edges == 1 & node.lasso$Edges == 1) / NP 
    }
  }
}

TP_lambda =apply(TP_rate,2,sum)
plot(TP_lambda, ylab="sun of TP", xlab="lambda", pch=19, type="b")
points(which.max(TP_lambda), TP_lambda[which.max(TP_lambda)], col="red", cex=2, pch=20)
lambda.max =grid[which.max(TP_lambda)] # 0.5994843
node.lasso <- data.frame(df)
node.lasso['Edges'] = rep(0,p)

for (i in 1:p){
  
  if(i==1){
    x <- X[,2:p]
    y <- X[,1]
  } else if(i==p){
    x <- X[,1:p-1]
    y <- X[,p]
  }else{
    x <- X[,c(1:(i-1),(i+1):p)]
    y <- X[,i]
  }
  
  for (j in 2:p){
    if (i==1){
      if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
        node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
    } else if(i==p){
      if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
        node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
    } else{
      if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
        if(j<=i){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
        } else{
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
        }
      }
    }
  }
}

###Implementing "Joint" and "Or" Rule to retrive edges

node.lasso1 <- data.frame(df)
node.lasso2 <- data.frame(df)

node.lasso1['Edges'] = rep(0,p)
node.lasso2['Edges'] = rep(0,p)

for (i in 1:p){
  for (j in 1:i){
    if ((node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] == 1) & (node.lasso["Edges"][node.lasso['Nodes'] == paste(j,i,sep=",")] == 1) ){
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 1
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 1
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 1
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 1
    } else if((node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] == 1) | (node.lasso["Edges"][node.lasso['Nodes'] == paste(j,i,sep=",")] == 1) ){
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 1
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 1
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 0
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 0
    } else{
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 0
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 0
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 0
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 0
    }
  }
}

# Proportion of correct classification

nwlasso1_mean_50 = mean(df['Edges'] == node.lasso1['Edges'])
nwlasso2_mean_50 = mean(df['Edges'] == node.lasso2['Edges'])

#ROC
roc.lasso <- data.frame(df)
value <- 10^seq(-2,2,length.out=10)
threshold = sort(unique(c(0, value)))
nth = length(threshold)
FP1_50_50 = TP1_50_50 = rep(0, nth)

for (k in 1:nth){
  roc.lasso['Edges'] = rep(0,p)
  for (i in 1:p){
    if(i==1){
      x <- X[,2:p]
      y <- X[,1]
    } else if(i==p){
      x <- X[,1:p-1]
      y <- X[,p]
    }else{
      x <- X[,c(1:(i-1),(i+1):p)]
      y <- X[,i]
    }
    fit.lasso <-glmnet(x,y, lambda = threshold[k])
    for (j in 2:p){
      if (i==1){
        if (coef(fit.lasso)[j]!=0){
          roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
      } else if(i==p){
        if (coef(fit.lasso)[j]!=0){
          roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
      } else{
        if (coef(fit.lasso)[j]!=0){
          if(j<=i){
            roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
          } else{
            roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j,sep=",")] = 1
          }
        }
      }
    }
  }
  TP1_50_50[k] = sum(df$Edges == 1 & roc.lasso$Edges == 1) / NP # True Positive Rate
  FP1_50_50[k] = sum(df$Edges == 0 & roc.lasso$Edges == 1) / NN # False Positive Rate
}


#6.1.2 Graphical Lasso Approach
  cv.glasso <-CVglasso(X)
  cv.glasso$Tuning
  cv.omega = cv.glasso$Omega
  graph.lasso['Edges'] = rep(0,p)
  for (i in 1:p){
    for (j in 1:p){
      if(cv.omega[i,j]!=0){
        graph.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
      }
    }
  }
Glasso_mean_50 = mean(df['Edges'] == graph.lasso['Edges'])

# ROC
roc.graph <- data.frame(df)
threshold = sort(unique(c(0, cv.glasso$Lambdas)))
nth = length(threshold)
FP.g.50 = TP.g.50 = rep(0, nth)

for (k in 1:nth){
  roc.graph['Edges'] = rep(0,p)
  fit.graph <- glasso(s, rho= threshold[k])
  omega <- fit.graph$wi
  for (i in 1:p){
    for (j in 1:p){
      if(omega[i,j]!=0){
        roc.graph["Edges"][roc.graph['Nodes'] == paste(i,j,sep=",")] = 1
      }
    }
  }
  
  
  TP.g.50[k] = sum(df$Edges == 1 & roc.graph$Edges == 1) / NP # True Positive Rate
  FP.g.50[k] = sum(df$Edges == 0 & roc.graph$Edges == 1) / NN # False Positive Rate
}


##6.2 n=p (p=100, n=100)
#6.2.1 Node-Wise Lasso Approach
set.seed(7)
X <- mvrnorm(n=100, mu=rep(0,100), Sigma = sigma)
node.lasso['Edges'] = rep(0,p)
grid <- 10^seq(-2,2, length=10) 
TP_rate = matrix(0,p,10)
NP = sum(df$Edges == 1) 

for (g in 1:10){
  for (i in 1:p){
    if(i==1){
      x <- X[,2:p]
      y <- X[,1]
    } else if(i==p){
      x <- X[,1:p-1]
      y <- X[,p]
    }else{
      x <- X[,c(1:(i-1),(i+1):p)]
      y <- X[,i]
    }
    for (j in 2:p){
      if (i==1){
        if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
      } else if(i==p){
        if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
      } else{
        if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
          if(j<=i){
            node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
          } else{
            node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
          }
        }
      }
      
      TP_rate[i,g] = sum(df$Edges == 1 & node.lasso$Edges == 1) / NP 
    }
  }
}

TP_lambda =apply(TP_rate,2,sum)
plot(TP_lambda, ylab="sun of TP", xlab="lambda", pch=19, type="b")
points(which.max(TP_lambda), TP_lambda[which.max(TP_lambda)], col="red", cex=2, pch=20)
lambda.max =grid[which.max(TP_lambda)] 
node.lasso <- data.frame(df)
node.lasso['Edges'] = rep(0,p)

for (i in 1:p){
  
  if(i==1){
    x <- X[,2:p]
    y <- X[,1]
  } else if(i==p){
    x <- X[,1:p-1]
    y <- X[,p]
  }else{
    x <- X[,c(1:(i-1),(i+1):p)]
    y <- X[,i]
  }
  
  for (j in 2:p){
    if (i==1){
      if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
        node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
    } else if(i==p){
      if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
        node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
    } else{
      if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
        if(j<=i){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
        } else{
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
        }
      }
    }
  }
}

###Implementing "Joint" and "Or" Rule to retrive edges

node.lasso1 <- data.frame(df)
node.lasso2 <- data.frame(df)

node.lasso1['Edges'] = rep(0,p)
node.lasso2['Edges'] = rep(0,p)

for (i in 1:p){
  for (j in 1:i){
    if ((node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] == 1) & (node.lasso["Edges"][node.lasso['Nodes'] == paste(j,i,sep=",")] == 1) ){
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 1
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 1
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 1
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 1
    } else if((node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] == 1) | (node.lasso["Edges"][node.lasso['Nodes'] == paste(j,i,sep=",")] == 1) ){
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 1
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 1
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 0
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 0
    } else{
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 0
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 0
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 0
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 0
    }
  }
}



# Proportion of correct classification

nwlasso1_mean_100 = mean(df['Edges'] == node.lasso1['Edges'])
nwlasso2_mean_100 = mean(df['Edges'] == node.lasso2['Edges'])

#ROC
roc.lasso <- data.frame(df)
value <- 10^seq(-2,2,length.out=10)
threshold = sort(unique(c(0, value)))
nth = length(threshold)
FP1_100 = TP1_100 = rep(0, nth)


for (k in 1:nth){
  roc.lasso['Edges'] = rep(0,p)
  
  for (i in 1:p){
    if(i==1){
      x <- X[,2:p]
      y <- X[,1]
    } else if(i==p){
      x <- X[,1:p-1]
      y <- X[,p]
    }else{
      x <- X[,c(1:(i-1),(i+1):p)]
      y <- X[,i]
    }
    fit.lasso <-glmnet(x,y, lambda = threshold[k])
    for (j in 2:p){
      if (i==1){
        if (coef(fit.lasso)[j]!=0){
          roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
      } else if(i==p){
        if (coef(fit.lasso)[j]!=0){
          roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
      } else{
        if (coef(fit.lasso)[j]!=0){
          if(j<=i){
            roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
          } else{
            roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j,sep=",")] = 1
          }
        }
      }
    }
  }
  TP1_100[k] = sum(df$Edges == 1 & roc.lasso$Edges == 1) / NP # True Positive Rate
  FP1_100[k] = sum(df$Edges == 0 & roc.lasso$Edges == 1) / NN # False Positive Rate
}


#6.2.2 Graphical Lasso Approach
cv.glasso <-CVglasso(X)
cv.glasso$Tuning
cv.omega = cv.glasso$Omega
graph.lasso['Edges'] = rep(0,p)
for (i in 1:p){
  for (j in 1:p){
    if(cv.omega[i,j]!=0){
      graph.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
    }
  }
}
 
Glasso_mean_100 = mean(df['Edges'] == graph.lasso['Edges'])

# ROC
roc.graph <- data.frame(df)
threshold = sort(unique(c(0, cv.glasso$Lambdas)))
nth = length(threshold)
FP.g.100 = TP.g.100 = rep(0, nth)

for (k in 1:nth){
  roc.graph['Edges'] = rep(0,p)
  fit.graph <- glasso(s, rho= threshold[k])
  omega <- fit.graph$wi
  for (i in 1:p){
    for (j in 1:p){
      if(omega[i,j]!=0){
        roc.graph["Edges"][roc.graph['Nodes'] == paste(i,j,sep=",")] = 1
      }
    }
  }
  TP.g.100[k] = sum(df$Edges == 1 & roc.graph$Edges == 1) / NP # True Positive Rate
  FP.g.100[k] = sum(df$Edges == 0 & roc.graph$Edges == 1) / NN # False Positive Rate
}


##6.3 n>p (p=100, n=200)
#6.2.1 Node-Wise Lasso Approach
# The orginal case
nwlasso1_mean_200
nwlasso2_mean_200
#6.3.2 Graphical Lasso Approach
# The orginal case
Glasso_mean_200


##6.4 n>>p (p=100, n=1000)
#6.4.1 Node-Wise Lasso Approach
set.seed(7)
X <- mvrnorm(n=1000, mu=rep(0,100), Sigma = sigma)
node.lasso['Edges'] = rep(0,p)
grid <- 10^seq(-2,2, length=10) 
TP_rate = matrix(0,p,10)
NP = sum(df$Edges == 1) 

for (g in 1:10){
  for (i in 1:p){
    
    if(i==1){
      x <- X[,2:p]
      y <- X[,1]
    } else if(i==p){
      x <- X[,1:p-1]
      y <- X[,p]
    }else{
      x <- X[,c(1:(i-1),(i+1):p)]
      y <- X[,i]
    }
    for (j in 2:p){
      if (i==1){
        if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
      } else if(i==p){
        if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
      } else{
        if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
          if(j<=i){
            node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
          } else{
            node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
          }
        }
      }
      
      TP_rate[i,g] = sum(df$Edges == 1 & node.lasso$Edges == 1) / NP 
    }
  }
}

TP_lambda =apply(TP_rate,2,sum)
plot(TP_lambda, ylab="sun of TP", xlab="lambda", pch=19, type="b")
points(which.max(TP_lambda), TP_lambda[which.max(TP_lambda)], col="red", cex=2, pch=20)
lambda.max =grid[which.max(TP_lambda)] 
node.lasso <- data.frame(df)
node.lasso['Edges'] = rep(0,p)

for (i in 1:p){
  
  if(i==1){
    x <- X[,2:p]
    y <- X[,1]
  } else if(i==p){
    x <- X[,1:p-1]
    y <- X[,p]
  }else{
    x <- X[,c(1:(i-1),(i+1):p)]
    y <- X[,i]
  }
  
  for (j in 2:p){
    if (i==1){
      if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
        node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
    } else if(i==p){
      if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
        node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
    } else{
      if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
        if(j<=i){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
        } else{
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
        }
      }
    }
  }
}

###Implementing "Joint" and "Or" Rule to retrive edges

node.lasso1 <- data.frame(df)
node.lasso2 <- data.frame(df)

node.lasso1['Edges'] = rep(0,p)
node.lasso2['Edges'] = rep(0,p)

for (i in 1:p){
  for (j in 1:i){
    if ((node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] == 1) & (node.lasso["Edges"][node.lasso['Nodes'] == paste(j,i,sep=",")] == 1) ){
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 1
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 1
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 1
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 1
    } else if((node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] == 1) | (node.lasso["Edges"][node.lasso['Nodes'] == paste(j,i,sep=",")] == 1) ){
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 1
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 1
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 0
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 0
    } else{
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 0
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 0
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 0
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 0
    }
  }
}

# Proportion of correct classification
nwlasso1_mean_1000 = mean(df['Edges'] == node.lasso1['Edges'])
nwlasso1_mean_1000 = mean(df['Edges'] == node.lasso2['Edges'])

#ROC
roc.lasso <- data.frame(df)
value <- 10^seq(-2,2,length.out=10)
threshold = sort(unique(c(0, value)))
nth = length(threshold)
FP1_1000 = TP1_1000 = rep(0, nth)

for (k in 1:nth){
  roc.lasso['Edges'] = rep(0,p)
  
  for (i in 1:p){
    if(i==1){
      x <- X[,2:p]
      y <- X[,1]
    } else if(i==p){
      x <- X[,1:p-1]
      y <- X[,p]
    }else{
      x <- X[,c(1:(i-1),(i+1):p)]
      y <- X[,i]
    }
    fit.lasso <-glmnet(x,y, lambda = threshold[k])
    for (j in 2:p){
      if (i==1){
        if (coef(fit.lasso)[j]!=0){
          roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
      } else if(i==p){
        if (coef(fit.lasso)[j]!=0){
          roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
      } else{
        if (coef(fit.lasso)[j]!=0){
          if(j<=i){
            roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
          } else{
            roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j,sep=",")] = 1
          }
        }
      }
    }
  }
  TP1_1000[k] = sum(df$Edges == 1 & roc.lasso$Edges == 1) / NP # True Positive Rate
  FP1_1000[k] = sum(df$Edges == 0 & roc.lasso$Edges == 1) / NN # False Positive Rate
}


#6.4.2 Graphical Lasso Approach
cv.glasso <-CVglasso(X)
cv.glasso$Tuning
cv.omega = cv.glasso$Omega
graph.lasso['Edges'] = rep(0,p)
for (i in 1:p){
  for (j in 1:p){
    if(cv.omega[i,j]!=0){
      graph.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
    }
  }
}

Glasso_mean_1000 = mean(df['Edges'] == graph.lasso['Edges'])

# ROC
roc.graph <- data.frame(df)
threshold = sort(unique(c(0, cv.glasso$Lambdas)))
nth = length(threshold)
FP.g.1000 = TP.g.1000 = rep(0, nth)

for (k in 1:nth){
  roc.graph['Edges'] = rep(0,p)
  fit.graph <- glasso(s, rho= threshold[k])
  omega <- fit.graph$wi
  for (i in 1:p){
    for (j in 1:p){
      if(omega[i,j]!=0){
        roc.graph["Edges"][roc.graph['Nodes'] == paste(i,j,sep=",")] = 1
      }
    }
  }
  
  
  TP.g.1000[k] = sum(df$Edges == 1 & roc.graph$Edges == 1) / NP # True Positive Rate
  FP.g.1000[k] = sum(df$Edges == 0 & roc.graph$Edges == 1) / NN # False Positive Rate
}


#6.5 0.5 with probability 0.2 or 0 with probability 0.8.
p <- 100
set.seed(7)
B <- matrix(sample(c(0,0.5), p^2, replace=TRUE, prob=c(0.8,0.2)), p)
i <- lower.tri(B)
B[i] <- t(B)[i]
diag(B) <- rep(0, p)
delta <- -min(eigen(B, symmetric=TRUE, only.values=TRUE)$values)
delta 
delta <- 4

theta = B + delta*diag(p)
is.positive.definite(theta)
theta <- cov2cor(theta)

df = data.frame(paste(rep(1:p,each =p),1:p,sep=","),rep(0,p*p))
names(df)<- c("Nodes","Edges")
for (i in 1:p){
  for (j in 1:p){
    if(theta[i,j]!=0){
      df["Edges"][df['Nodes'] == paste(i,j,sep=",")] = 1
    }
  }
}
node.lasso <- data.frame(df)
graph.lasso <- data.frame(df)
sigma = solve(theta)

#6.5.1 Node-Wise Lasso Approach
set.seed(7)
X <- mvrnorm(n=200, mu=rep(0,100), Sigma = sigma)
node.lasso['Edges'] = rep(0,p)
grid <- 10^seq(-2,2, length=10) 
TP_rate = matrix(0,p,10)
NP = sum(df$Edges == 1) 

for (g in 1:10){
  for (i in 1:p){
    
    if(i==1){
      x <- X[,2:p]
      y <- X[,1]
    } else if(i==p){
      x <- X[,1:p-1]
      y <- X[,p]
    }else{
      x <- X[,c(1:(i-1),(i+1):p)]
      y <- X[,i]
    }
    for (j in 2:p){
      if (i==1){
        if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
      } else if(i==p){
        if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
      } else{
        if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
          if(j<=i){
            node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
          } else{
            node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
          }
        }
      }
      
      TP_rate[i,g] = sum(df$Edges == 1 & node.lasso$Edges == 1) / NP 
    }
  }
}

TP_lambda =apply(TP_rate,2,sum)
plot(TP_lambda, ylab="sun of TP", xlab="lambda", pch=19, type="b")
points(which.max(TP_lambda), TP_lambda[which.max(TP_lambda)], col="red", cex=2, pch=20)
lambda.max =grid[which.max(TP_lambda)] 
node.lasso <- data.frame(df)
node.lasso['Edges'] = rep(0,p)

for (i in 1:p){
  
  if(i==1){
    x <- X[,2:p]
    y <- X[,1]
  } else if(i==p){
    x <- X[,1:p-1]
    y <- X[,p]
  }else{
    x <- X[,c(1:(i-1),(i+1):p)]
    y <- X[,i]
  }
  
  for (j in 2:p){
    if (i==1){
      if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
        node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
    } else if(i==p){
      if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
        node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
    } else{
      if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
        if(j<=i){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
        } else{
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
        }
      }
    }
  }
}

###Implementing "Joint" and "Or" Rule to retrive edges

node.lasso1 <- data.frame(df)
node.lasso2 <- data.frame(df)

node.lasso1['Edges'] = rep(0,p)
node.lasso2['Edges'] = rep(0,p)

for (i in 1:p){
  for (j in 1:i){
    if ((node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] == 1) & (node.lasso["Edges"][node.lasso['Nodes'] == paste(j,i,sep=",")] == 1) ){
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 1
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 1
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 1
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 1
    } else if((node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] == 1) | (node.lasso["Edges"][node.lasso['Nodes'] == paste(j,i,sep=",")] == 1) ){
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 1
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 1
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 0
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 0
    } else{
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 0
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 0
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 0
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 0
    }
  }
}

# Proportion of correct classification
nwlasso1_mean_0.8=mean(df['Edges'] == node.lasso1['Edges'])
nwlasso2_mean_0.8=mean(df['Edges'] == node.lasso2['Edges'])

#ROC
roc.lasso <- data.frame(df)
value <- 10^seq(-2,0,length.out=10)
threshold = sort(unique(c(0, value)))
nth = length(threshold)
FP1.0.8 = TP1.0.8 = rep(0, nth)

for (k in 1:nth){
  roc.lasso['Edges'] = rep(0,p)
  
  for (i in 1:p){
    if(i==1){
      x <- X[,2:p]
      y <- X[,1]
    } else if(i==p){
      x <- X[,1:p-1]
      y <- X[,p]
    }else{
      x <- X[,c(1:(i-1),(i+1):p)]
      y <- X[,i]
    }
    fit.lasso <-glmnet(x,y, lambda = threshold[k])
    for (j in 2:p){
      if (i==1){
        if (coef(fit.lasso)[j]!=0){
          roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
      } else if(i==p){
        if (coef(fit.lasso)[j]!=0){
          roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
      } else{
        if (coef(fit.lasso)[j]!=0){
          if(j<=i){
            roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
          } else{
            roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j,sep=",")] = 1
          }
        }
      }
    }
  }
  TP1.0.8[k] = sum(df$Edges == 1 & roc.lasso$Edges == 1) / NP # True Positive Rate
  FP1.0.8[k] = sum(df$Edges == 0 & roc.lasso$Edges == 1) / NN # False Positive Rate
}

#6.5.2 Graphical Lasso Approach
  cv.glasso <-CVglasso(X)
  cv.glasso$Tuning
  cv.omega = cv.glasso$Omega
  graph.lasso['Edges'] = rep(0,p)
  for (i in 1:p){
    for (j in 1:p){
      if(cv.omega[i,j]!=0){
        graph.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
      }
    }
  }
Glasso_mean_0.8 = mean(df['Edges'] == graph.lasso['Edges'])

#plot ROC
roc.graph <- data.frame(df)
threshold = sort(unique(c(0, cv.glasso$Lambdas)))
nth = length(threshold)
FP.g.0.8 = TP.g.0.8 = rep(0, nth)

for (k in 1:nth){
  roc.graph['Edges'] = rep(0,p)
  fit.graph <- glasso(s, rho= threshold[k])
  omega <- fit.graph$wi
  for (i in 1:p){
    for (j in 1:p){
      if(omega[i,j]!=0){
        roc.graph["Edges"][roc.graph['Nodes'] == paste(i,j,sep=",")] = 1
      }
    }
  }
  TP.g.0.8[k] = sum(df$Edges == 1 & roc.graph$Edges == 1) / NP # True Positive Rate
  FP.g.0.8[k] = sum(df$Edges == 0 & roc.graph$Edges == 1) / NN # False Positive Rate
}


#6.6 0.5 with probability 0.05 or 0 with probability 0.95.
p <- 100
set.seed(7)
B <- matrix(sample(c(0,0.5), p^2, replace=TRUE, prob=c(0.95,0.05)), p)
i <- lower.tri(B)
B[i] <- t(B)[i]
diag(B) <- rep(0, p)
delta <- -min(eigen(B, symmetric=TRUE, only.values=TRUE)$values)
delta 
delta <- 4

theta = B + delta*diag(p)
is.positive.definite(theta)
theta <- cov2cor(theta)

df = data.frame(paste(rep(1:p,each =p),1:p,sep=","),rep(0,p*p))
names(df)<- c("Nodes","Edges")
for (i in 1:p){
  for (j in 1:p){
    if(theta[i,j]!=0){
      df["Edges"][df['Nodes'] == paste(i,j,sep=",")] = 1
    }
  }
}
node.lasso <- data.frame(df)
graph.lasso <- data.frame(df)
sigma = solve(theta)

#6.6.1 Node-Wise Lasso Approach
set.seed(7)
X <- mvrnorm(n=200, mu=rep(0,100), Sigma = sigma)
node.lasso['Edges'] = rep(0,p)
grid <- 10^seq(-2,2, length=10) 
TP_rate = matrix(0,p,10)
NP = sum(df$Edges == 1) 

for (g in 1:10){
  for (i in 1:p){
    
    if(i==1){
      x <- X[,2:p]
      y <- X[,1]
    } else if(i==p){
      x <- X[,1:p-1]
      y <- X[,p]
    }else{
      x <- X[,c(1:(i-1),(i+1):p)]
      y <- X[,i]
    }
    for (j in 2:p){
      if (i==1){
        if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
      } else if(i==p){
        if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
      } else{
        if (coef(glmnet(x,y, lambda=grid[g]))[j]!=0){
          if(j<=i){
            node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
          } else{
            node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
          }
        }
      }
      
      TP_rate[i,g] = sum(df$Edges == 1 & node.lasso$Edges == 1) / NP 
    }
  }
}

TP_lambda =apply(TP_rate,2,sum)
plot(TP_lambda, ylab="sun of TP", xlab="lambda", pch=19, type="b")
points(which.max(TP_lambda), TP_lambda[which.max(TP_lambda)], col="red", cex=2, pch=20)
lambda.max =grid[which.max(TP_lambda)] 
node.lasso <- data.frame(df)
node.lasso['Edges'] = rep(0,p)

for (i in 1:p){
  
  if(i==1){
    x <- X[,2:p]
    y <- X[,1]
  } else if(i==p){
    x <- X[,1:p-1]
    y <- X[,p]
  }else{
    x <- X[,c(1:(i-1),(i+1):p)]
    y <- X[,i]
  }
  
  for (j in 2:p){
    if (i==1){
      if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
        node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
    } else if(i==p){
      if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
        node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
    } else{
      if (coef(glmnet(x,y, lambda=lambda.max))[j]!=0){
        if(j<=i){
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
        } else{
          node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
        }
      }
    }
  }
}

###Implementing "Joint" and "Or" Rule to retrive edges

node.lasso1 <- data.frame(df)
node.lasso2 <- data.frame(df)

node.lasso1['Edges'] = rep(0,p)
node.lasso2['Edges'] = rep(0,p)

for (i in 1:p){
  for (j in 1:i){
    if ((node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] == 1) & (node.lasso["Edges"][node.lasso['Nodes'] == paste(j,i,sep=",")] == 1) ){
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 1
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 1
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 1
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 1
    } else if((node.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] == 1) | (node.lasso["Edges"][node.lasso['Nodes'] == paste(j,i,sep=",")] == 1) ){
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 1
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 1
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 0
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 0
    } else{
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(i,j,sep=",")] = 0
      node.lasso1["Edges"][node.lasso1['Nodes'] == paste(j,i,sep=",")] = 0
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(i,j,sep=",")] = 0
      node.lasso2["Edges"][node.lasso2['Nodes'] == paste(j,i,sep=",")] = 0
    }
  }
}

# Proportion of correct classification
nwlasso1_mean_0.95 = mean(df['Edges'] == node.lasso1['Edges'])
nwlasso2_mean_0.95 = mean(df['Edges'] == node.lasso2['Edges'])

# ROC
roc.lasso <- data.frame(df)
value <- 10^seq(-2,0,length.out=10)
threshold = sort(unique(c(0, value)))
nth = length(threshold)
FP1.0.95 = TP1.0.95 = rep(0, nth)

for (k in 1:nth){
  roc.lasso['Edges'] = rep(0,p)
  
  for (i in 1:p){
    if(i==1){
      x <- X[,2:p]
      y <- X[,1]
    } else if(i==p){
      x <- X[,1:p-1]
      y <- X[,p]
    }else{
      x <- X[,c(1:(i-1),(i+1):p)]
      y <- X[,i]
    }
    fit.lasso <-glmnet(x,y, lambda = threshold[k])
    for (j in 2:p){
      if (i==1){
        if (coef(fit.lasso)[j]!=0){
          roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j,sep=",")] = 1 }
      } else if(i==p){
        if (coef(fit.lasso)[j]!=0){
          roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1}
      } else{
        if (coef(fit.lasso)[j]!=0){
          if(j<=i){
            roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j-1,sep=",")] = 1
          } else{
            roc.lasso["Edges"][roc.lasso['Nodes'] == paste(i,j,sep=",")] = 1
          }
        }
      }
    }
  }
  TP1.0.95[k] = sum(df$Edges == 1 & roc.lasso$Edges == 1) / NP # True Positive Rate
  FP1.0.95[k] = sum(df$Edges == 0 & roc.lasso$Edges == 1) / NN # False Positive Rate
}

#6.6.2 Graphical Lasso Approach
cv.glasso <-CVglasso(X)
cv.glasso$Tuning
cv.omega = cv.glasso$Omega
graph.lasso['Edges'] = rep(0,p)
for (i in 1:p){
  for (j in 1:p){
    if(cv.omega[i,j]!=0){
      graph.lasso["Edges"][node.lasso['Nodes'] == paste(i,j,sep=",")] = 1
    }
  }
}

Glasso_mean_0.95 = mean(df['Edges'] == graph.lasso['Edges'])

#plot ROC
roc.graph <- data.frame(df)
threshold = sort(unique(c(0, cv.glasso$Lambdas)))
nth = length(threshold)
FP.g.0.95 = TP.g.0.95 = rep(0, nth)

for (k in 1:nth){
  roc.graph['Edges'] = rep(0,p)
  fit.graph <- glasso(s, rho= threshold[k])
  omega <- fit.graph$wi
  for (i in 1:p){
    for (j in 1:p){
      if(omega[i,j]!=0){
        roc.graph["Edges"][roc.graph['Nodes'] == paste(i,j,sep=",")] = 1
      }
    }
  }
  TP.g.0.95[k] = sum(df$Edges == 1 & roc.graph$Edges == 1) / NP # True Positive Rate
  FP.g.0.95[k] = sum(df$Edges == 0 & roc.graph$Edges == 1) / NN # False Positive Rate
}


#6.5 Compare the results
#6.5.1 Node wise lasso
par(mfrow=c(1,2))
plot(
  x = c(0,1),
  y = c(0,1),
  type = "n",
  main = "Node-wise Lasso ROC",
  xlab = "False Positive Rate",
  ylab = "True Positive Rate"
)
lines(FP1_50_50, TP1_50_50, col = 1, lwd = 3)
lines(FP1_100, TP1_100, col = 2, lwd = 3)
lines(FP1, TP1, col = 3, lwd = 3)
lines(FP1_1000, TP1_1000, col = 4, lwd = 3)
legend("bottomright", legend = c("n=50,", 
                                 "n=100",
                                 "n=200",
                                 "n=1000"), col = c(1:4), fill = c(1:4))

plot(
  x = c(0,1),
  y = c(0,1),
  type = "n",
  main = "Node-wise Lasso ROC",
  xlab = "False Positive Rate",
  ylab = "True Positive Rate"
)
lines(FP1.0.8, TP1.0.8, col = 1, lwd = 3)
lines(FP1, TP1, col = 2, lwd = 3)
lines(FP1.0.95, TP1.0.95, col = 3, lwd = 3)
legend("bottomright", legend = c("p=0.8",
                                 "p=0.9",
                                 "p=0.95"), col = c(1:3), fill = c(1:3))

#6.5.2 Glasso
plot(
  x = c(0,1),
  y = c(0,1),
  type = "n",
  main = "Graphical Lasso ROC",
  xlab = "False Positive Rate",
  ylab = "True Positive Rate"
)

lines(FP.g.50, TP.g.50, col = 1, lwd = 3)
lines(FP.g.100, TP.g.100, col = 2, lwd = 3)
lines(FP.g.200, TP.g.200, col = 3, lwd = 3,)
lines(FP.g.1000, TP.g.1000, col = 4, lwd = 3)
legend("bottomright", legend = c("n=50,", 
                                 "n=100",
                                 "n=200",
                                 "n=1000"), col = c(1:4), fill = c(1:4))
plot(
  x = c(0,1),
  y = c(0,1),
  type = "n",
  main = "Graphical Lasso ROC",
  xlab = "False Positive Rate",
  ylab = "True Positive Rate"
)
lines(FP.g.0.8, TP.g.0.8, col = 1, lwd = 3)
lines(FP.g.200, TP.g.200, col = 2, lwd = 3,)
lines(FP.g.0.95, TP.g.0.95, col = 3, lwd = 3)
legend("bottomright", legend = c("p=0.8",
                                 "p=0.9",
                                 "p=0.95"), col = c(1:3), fill = c(1:3))
