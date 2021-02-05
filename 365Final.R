#downloading data
x <- get.hist.quote(instrument = "FSREX", start = "2014-12-01", end = "2019-12-31",
                    quote = "AdjClose", compression = "m")
fsrex <- as.vector (x)
head (fsrex)

x <- get.hist.quote(instrument = "RYTLX", start = "2014-12-01", end = "2019-12-31",
                    quote = "AdjClose", compression = "m")
rytlx <- as.vector (x)
head (rytlx)

x <- get.hist.quote(instrument = "PRHSX", start = "2014-12-01", end = "2019-12-31",
                    quote = "AdjClose", compression = "m")
prhsx <- as.vector (x)
head (prhsx)

x <- get.hist.quote(instrument = "XBTOX", start = "2014-12-01", end = "2019-12-31",
                    quote = "AdjClose", compression = "m")
xbtox <- as.vector (x)
head (xbtox)

x <- get.hist.quote(instrument = "^GSPC", start = "2014-12-01", end = "2019-12-31",
                    quote = "AdjClose", compression = "m")
sp500 <- as.vector (x)
head (sp500)

#retrieving risk free data
rff <- c (0.03, 0.03, 0.02, 0.03, 0.02, 0.02, 0.02, 0.03, 0.07, 0.02, 0.02, 0.12, 0.23, 0.26, 0.31,
          0.29, 0.23, 0.27, 0.27, 0.30, 0.30, 0.29, 0.33, 0.45, 0.51, 0.51, 0.52, 0.74, 0.80, 0.89, 
          0.98, 1.07, 1.01, 1.03, 1.07, 1.23, 1.32, 1.41, 1.57, 1.70, 1.76, 1.86, 1.90, 1.96, 2.03,
          2.13, 2.25, 2.33, 2.37, 2.37, 2.39, 2.40, 2.38, 2.35, 2.17, 2.10, 1.95, 1.89, 1.65, 1.54,
          1.54)
head (rff)
length (rff)
rfree<-(1 + rff/100)^(1/12)-1
head (rfree)

rfree <- rfree[-1]
head (rfree)

# 1 sentence description
#FSREX - Fidelity Series Real Estate Income Fund - invests in preferred and common stocks of real
#estate investment trusts and debt securities of real estate entities.
#RYTLX - Rydex Telecommunications - invests in equity securities of telecommunications companies and 
#derivatives.
#PRHSX - T. Rowe Price Health Sciences Fund - buys stocks in the medical and health care industries
#XBTOX - John Hancock Financial Opportunities Fund - invests in equity securities of US regional banks

#calculating returns
fsrex.ret <- (fsrex[-1]-fsrex[-61])/fsrex[-61]
head (fsrex.ret)

rytlx.ret <- (rytlx[-1]-rytlx[-61])/rytlx[-61]
head (rytlx.ret)

prhsx.ret <- (prhsx[-1]-prhsx[-61])/prhsx[-61]
head (prhsx.ret)

xbtox.ret <- (xbtox[-1]-xbtox[-61])/xbtox[-61]
head (xbtox.ret)

sp500.ret <- (sp500[-1]-sp500[-61])/sp500[-61]
head (sp500.ret)

#defining excess returns
fsrex.ex <- fsrex.ret - rfree
rytlx.ex <- rytlx.ret - rfree
prhsx.ex <- prhsx.ret - rfree
xbtox.ex <- xbtox.ret - rfree
sp500.ex <- sp500.ret - rfree
mean (sp500.ex)
sd(sp500.ex)

#combining stocks
funds <- cbind(fsrex.ex, rytlx.ex, prhsx.ex, xbtox.ex)
colnames (funds) <- c("FSREX", "RYTLX", "PRHSX", "XBTOX")

#4 calculating mean and sd
apply (funds, 2, mean)
apply (funds, 2, sd)
#prhsx has the highest mean, but also highest sd. Fsrex has lowest sd but low mean. Rytlx has 
#lowest mean and mid sd, arguably the worst

#5 sample covariance matrix
Smat <- cov(funds)
Smat <- as.matrix (Smat)
Corrmat <- as.matrix(cor(funds))

#6 estimate weights of risk-averse portfolio, lambda = 5
funds.mean <- apply (funds, 2, mean)
funds.mean <- as.vector (funds.mean)
funds.sd <- as.vector(apply(funds, 2, sd))

w0 <- solve(Smat, c(1,1,1,1)) 
w_mv <- w0/sum(w0)
w_mv

m <- sum (w_mv*funds.mean)
vbar <- solve (Smat, funds.mean - m*c(1,1,1,1))
vbar
RA5.wgts <- w_mv + vbar/5
RA5.portmean <- sum(RA5.wgts*funds.mean)
RA5.portmean
RA5.portsd <- (RA5.wgts%*%Smat%*%RA5.wgts)^.5
RA5.portsd


#7 non-negative weights
library (quadprog)
A <- cbind(c(1,1,1,1), diag(4))
t(A)
b <- c(1,0,0,0,0)
qpsol <- solve.QP (Dmat = (5)*Smat, dvec = funds.mean, Amat = A, bvec = b, meq = 1)
RA5.nonneg.wghts <- as.vector (qpsol$solution)
RA5.nonneg.mean <- RA5.nonneg.wghts%*%funds.mean
RA5.nonneg.sd <- ((qpsol$solution)%*%Smat%*%(qpsol$solution))^.5

#8 Risk aversion parameter = 20
RA20.wgts <- w_mv + vbar/20
RA20.portmean <- sum(RA20.wgts*funds.mean)
RA20.portmean
RA20.portsd <- (RA20.wgts%*%Smat%*%RA20.wgts)^0.5
RA20.portsd

#9 see notes

#10 finding the tangency portfolio
tan.wgts <- solve(Smat, funds.mean)/sum(solve(Smat, funds.mean))
tan.wgts
tan.mean <- sum (tan.wgts*funds.mean)
tan.mean
tan.sd <- (tan.wgts%*%Smat%*%tan.wgts)^0.5
tan.sd

#11 finding Sharpe Ratios
RA5.SR <- RA5.portmean/RA5.portsd
RA5.nonneg.SR <- RA5.nonneg.mean/RA5.nonneg.sd
RA20.SR <- RA20.portmean/RA20.portsd
tan.SR <- tan.mean/tan.sd
tan.SR

#12 finding market model parameters - see notes
fsrex.mm <- lm(fsrex.ex~sp500.ex)
rytlx.mm <- lm(rytlx.ex~sp500.ex)
prhsx.mm <- lm(prhsx.ex~sp500.ex)
xbtox.mm <- lm(xbtox.ex~sp500.ex)
funds.mm <- c (fsrex.mm, rytlx.mm, prhsx.mm, xbtox.mm)
funds.mm$coefficients
lm(funds ~ sp500.ex)$coefficients

#13 testing if funds are mispriced
f.alphapval<-function(y)
  + {summary(lm(y~sp500.ex))$coefficients[1, 4]}
funds.pval <- apply (funds, 2, f.alphapval) 

#14 finding non-market risk
fsrex.nmr <- summary (fsrex.mm)$sigma
rytlx.nmr <- summary (rytlx.mm)$sigma
prhsx.nmr <- summary (prhsx.mm)$sigma
xbtox.nmr <- summary (xbtox.mm)$sigma
funds.nmr <- c(fsrex.nmr, rytlx.nmr, prhsx.nmr, xbtox.nmr)
funds.nmr

#15 proportion of risk explained by the market
f.rsq <- function(y)
  + {summary(lm(y~sp500.ex))$r.squared}
apply (funds, 2, f.rsq)

#16 finding portfolio with minimum variance and beta = 1
library (quadprog)
A <- cbind(c(1,1,1,1), funds.beta)
t(A)
b <- c(1,1)
qpsol <- solve.QP (Dmat = (2)*Smat, dvec = c(0,0,0,0), Amat = A, bvec = b, meq = 2)
qpsol$solution
B1.portmean <- sum (qpsol$solution * funds.mean)
B1.portmean
B1.portsd <- ((qpsol$solution)%*%Smat%*%(qpsol$solution))^0.5
B1.portsd

qpsol2 <- solve.QP (Dmat = Smat, dvec = c(0,0,0,0), Amat = A, bvec = b, meq = 2)
qpsol2$solution
B1.portmean2 <- sum (qpsol2$solution * funds.mean)
B1.portmean2
B1.portsd2 <- ((qpsol2$solution)%*%Smat%*%(qpsol2$solution))^0.5
B1.portsd2


#17 parameters of market model from portfolio found in 16
B1port <- 0.01122*fsrex.ex + 0.56978*rytlx.ex + 0.13096*prhsx.ex + 0.28804*xbtox.ex
B1port.mm <- lm(B1port~sp500.ex)
summary (B1port.mm)

#18 finding Sharpe, Treynor, and Appraial ratio for each fund
fsrex.SR <- mean(fsrex.ex)/sd(fsrex.ex)
rytlx.SR <- mean(rytlx.ex)/sd(rytlx.ex)
prhsx.SR <- mean(prhsx.ex)/sd(prhsx.ex)
xbtox.SR <- mean(xbtox.ex)/sd(xbtox.ex)
funds.SR <- c(fsrex.SR, rytlx.SR, prhsx.SR, xbtox.SR)
funds.SR
#better way to do this
f.SR <- function(y) {mean(y)/sd(y)}
funds.SR <- apply (funds, 2, f.SR)
#Treynor Ratio
f.TR <- function(y, z) {y*sd(sp500.ex)/cor(z, sp500.ex)}
funds.TR <- f.TR (funds.SR, funds)
funds.TR
funds.mm <- lm (funds~sp500.ex)
funds.alpha<-funds.mm$coefficients[1,] 
funds.beta<-funds.mm$coefficients[2,]
funds.TR <- apply (funds, 2, mean)/ funds.beta
funds.TR


#appraial ratio
f.AR <- function(y) {summary (lm(y~sp500.ex))$coefficients[1,1]/summary(lm(y~sp500.ex))$sigma}
funds.AR <- apply (funds, 2, f.AR)


