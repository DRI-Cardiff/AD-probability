gwas<-read.table(file="GWAS_SNPs.txt", header=T, sep=" ")
M<-nrow(gwas)

OdR<-gwas$OR
log(OdR)
f0<-gwas$AltAlleleFreq
f1<-OdR*f0/(1-f0*(1-OdR))
lor <- log(f1*(1-f0)/f0/(1-f1))  # log odds ratio

Nc <- 10000  # number of cases
Nn <- 10000  # number of controls
K <- Nc/(Nc+Nn)  # prevalence

# create random genotypes based on allele frequencies
gen <- matrix(nrow=Nc+Nn, ncol=M)
for (m in (1:M)) {
  gen[1:Nc, m] <- rbinom(Nc, 1, f1[m])+rbinom(Nc, 1, f1[m]);
  gen[(Nc+1):(Nc+Nn), m] <- rbinom(Nc, 1, f0[m])+rbinom(Nc, 1, f0[m])
}

status <- c(rep(1, Nc), rep(0, Nn))  # case/non-case status
PRS <- array(dim=Nc+Nn)
for (i in (1:(Nc+Nn))) PRS[i] <- sum(lor*gen[i,])  # polygenic risk score


# do logistic regression on the generated PRS from simulation
cdat <- glm(status ~ PRS, family=binomial)$coefficients

# calculate means, variances and st.dev. of PRS in cases and in non-cases
m1 <- mean(PRS[1:Nc])
m0 <- mean(PRS[(Nc+1):(Nc+Nn)])
v1 <- var(PRS[1:Nc])
v0 <- var(PRS[(Nc+1):(Nc+Nn)])
s1 <- sqrt(v1)
s0 <- sqrt(v0)

# hence population data
mp <- K*m1 + (1-K)*m0  # mean PRS in population
vp <- K*s1*s1 + (1-K)*s0*s0 + K*(1-K)*(m1 - m0)*(m1 - m0) # PRS variance in population

# theoretical regression coefficients
r0 <- (s1*s1 + (m0 - m1)*(m0 - m1))/s0/s0
r1 <- (s0*s0 + (m0 - m1)*(m0 - m1))/s1/s1
beta <- (m1 - m0)*(K*(1-K)*((r0+r1)/2 - 1) + K*s1*s1/s0/s0 + (1-K)*s0*s0/s1/s1)/vp
alpha <- log(K*s0/(1-K)/s1) + ((r0-1)*K + (1-r1)*(1-K))/2 - mp*beta


# plots:
grid <- (-40:40)/10
X <- sqrt(vp)*grid+mp
Y <- 1/(1 + exp(-alpha - beta*X))

plot(Y ~ grid, type="l", ylim=c(0,1), col="red", 
     main="Probability of disease", ylab="Probability", xlab="Standardised PRS")
Ydat <- 1/(1 + exp(-cdat[1] - cdat[2]*X))
points (Ydat ~ grid, col="red", pch=4, cex = 0.5)


GWASwithAPOE<-PRS






#PRS_noAPE
PRS_noAPOE<-PRS
for (i in (1:(Nc+Nn))) PRS_noAPOE[i] <- PRS[i]-lor[37]*gen[i,37]  #
cdat <- glm(status ~ PRS_noAPOE, family=binomial)$coefficients

# calculate means, variances and st.dev. of PRS in cases and in non-cases
m1 <- mean(PRS_noAPOE[1:Nc])
m0 <- mean(PRS_noAPOE[(Nc+1):(Nc+Nn)])
v1 <- var(PRS_noAPOE[1:Nc])
v0 <- var(PRS_noAPOE[(Nc+1):(Nc+Nn)])
s1 <- sqrt(v1)
s0 <- sqrt(v0)

# hence population data
mp <- K*m1 + (1-K)*m0  # mean PRS in population
vp <- K*s1*s1 + (1-K)*s0*s0 + K*(1-K)*(m1 - m0)*(m1 - m0) # PRS variance in population

# theoretical regression coefficients
r0 <- (s1*s1 + (m0 - m1)*(m0 - m1))/s0/s0
r1 <- (s0*s0 + (m0 - m1)*(m0 - m1))/s1/s1
beta <- (m1 - m0)*(K*(1-K)*((r0+r1)/2 - 1) + K*s1*s1/s0/s0 + (1-K)*s0*s0/s1/s1)/vp
alpha <- log(K*s0/(1-K)/s1) + ((r0-1)*K + (1-r1)*(1-K))/2 - mp*beta

c(alpha, beta) # for comparison

# plots:
# black = theoretical
# green = from logistic regression (actual genotype simulation done in beginning)
grid <- (-40:40)/10
X <- sqrt(vp)*grid+mp
Y <- 1/(1 + exp(-alpha - beta*X))

points(Y ~ grid, type="l", col="purple")
Ydat <- 1/(1 + exp(-cdat[1] - cdat[2]*X))
points (Ydat ~ grid, col="purple", pch=4, cex = 0.5)










#APOE only
APOE <- array(dim=Nc+Nn)
for (i in (1:(Nc+Nn))) APOE[i] <- gen[i,37]  # polygenic risk score

# do logistic regression on the generated PRS from simulation
cdat <- glm(status ~ APOE, family=binomial)$coefficients

# calculate means, variances and st.dev. of PRS in cases and in non-cases
m1 <- mean(APOE[1:Nc])
m0 <- mean(APOE[(Nc+1):(Nc+Nn)])
v1 <- var(APOE[1:Nc])
v0 <- var(APOE[(Nc+1):(Nc+Nn)])
s1 <- sqrt(v1)
s0 <- sqrt(v0)

# hence population data
mp <- K*m1 + (1-K)*m0  # mean PRS in population
vp <- K*s1*s1 + (1-K)*s0*s0 + K*(1-K)*(m1 - m0)*(m1 - m0) # PRS variance in population

# theoretical regression coefficients
r0 <- (s1*s1 + (m0 - m1)*(m0 - m1))/s0/s0
r1 <- (s0*s0 + (m0 - m1)*(m0 - m1))/s1/s1
beta <- (m1 - m0)*(K*(1-K)*((r0+r1)/2 - 1) + K*s1*s1/s0/s0 + (1-K)*s0*s0/s1/s1)/vp
alpha <- log(K*s0/(1-K)/s1) + ((r0-1)*K + (1-r1)*(1-K))/2 - mp*beta


# plots:
grid <- (-40:40)/10
X <- sqrt(vp)*grid+mp
Y <- 1/(1 + exp(-alpha - beta*X))

#RARE
RARE<-array(0, dim=Nc+Nn)
lor<-c(log(7.22), log(2.46), log(1.67),  log(1.43), log(0.68), log(0.189))
f1<-c(0.0062,    0.004,     0.014,      0.011,     0.006,      0.0013)
f0<-f1/(exp(lor)*(1-f1)+f1)
#SORL1 (OR=7.22, Bellenguez), TREM2 (OR=2.46 and OR=1.67, Sims), ABI3 (1.43, sims), 
#PLCG2 (OR=0.68, Sims), (APP, OR=0.189, Jonsson, Nature 488:96–99(2012))

M<-length(lor)
gen <- matrix(nrow=Nc+Nn, ncol=M)
for (m in (1:M)) {
  gen[1:Nc, m] <- rbinom(Nc, 1, f1[m])+rbinom(Nc, 1, f1[m]);
  gen[(Nc+1):(Nc+Nn), m] <- rbinom(Nc, 1, f0[m])+rbinom(Nc, 1, f0[m])
}
for (i in (1:(Nc+Nn))) RARE[i] <- sum(lor*gen[i,])  # polygenic risk score

GWASwithAPOEwithRARE<-GWASwithAPOE+RARE

cdat <- glm(status ~ GWASwithAPOEwithRARE, family=binomial)$coefficients

# calculate means, variances and st.dev. of PRS in cases and in non-cases
m1 <- mean(GWASwithAPOEwithRARE[1:Nc])
m0 <- mean(GWASwithAPOEwithRARE[(Nc+1):(Nc+Nn)])
v1 <- var(GWASwithAPOEwithRARE[1:Nc])
v0 <- var(GWASwithAPOEwithRARE[(Nc+1):(Nc+Nn)])
s1 <- sqrt(v1)
s0 <- sqrt(v0)

# hence population data
mp <- K*m1 + (1-K)*m0  # mean PRS in population
vp <- K*s1*s1 + (1-K)*s0*s0 + K*(1-K)*(m1 - m0)*(m1 - m0) # PRS variance in population

# theoretical regression coefficients
r0 <- (s1*s1 + (m0 - m1)*(m0 - m1))/s0/s0
r1 <- (s0*s0 + (m0 - m1)*(m0 - m1))/s1/s1
beta <- (m1 - m0)*(K*(1-K)*((r0+r1)/2 - 1) + K*s1*s1/s0/s0 + (1-K)*s0*s0/s1/s1)/vp
alpha <- log(K*s0/(1-K)/s1) + ((r0-1)*K + (1-r1)*(1-K))/2 - mp*beta

grid <- (-40:40)/10
X <- sqrt(vp)*grid+mp
Y <- 1/(1 + exp(-alpha - beta*X))








#Simulation of additional PRS
Nmar<-  c(5000,  10000, 25000,  20000,  20000, 20000)
Bmean<- c(0.005, 0.003, 0.002, 0.0015, 0.001,  0.0005)
Bsdev <-c(0.02,  0.015, 0.01,   0.01,   0.01,   0.01)
RAN   <-c(0.9,   0.8,   0.57,  0.7,    0.6,   0.5)
mafL <-0.01
mafU <- 0.45
sum(Nmar)

GWASwithAPOEwithRAREwithPRS<-GWASwithAPOE+RARE
for (SC in 1:length(Nmar))
{
  KK<-floor(Nmar[SC]*RAN[SC])
  KKminorAllele<-floor(KK*0.7)
  lor<-rnorm(KKminorAllele, Bmean[SC], Bsdev[SC])
  f1<- runif(KKminorAllele, mafL, mafU)
  f0<-f1/(exp(lor)*(1-f1)+f1)
  gen <- matrix(nrow=Nc+Nn, ncol=KKminorAllele)
  for (m in (1:KKminorAllele)) {
    gen[1:Nc, m] <- rbinom(Nc, 1, f1[m])+rbinom(Nc, 1, f1[m]);
    gen[(Nc+1):(Nc+Nn), m] <- rbinom(Nc, 1, f0[m])+rbinom(Nc, 1, f0[m])
  }
  for (i in (1:(Nc+Nn))) GWASwithAPOEwithRAREwithPRS[i] <- GWASwithAPOEwithRAREwithPRS[i]+sum(lor*gen[i,])
    
  KKmajorAllele<-KK-KKminorAllele
  lor<-rnorm(KKmajorAllele, Bmean[SC], Bsdev[SC])
  f1<- runif(KKmajorAllele, (1-mafU), (1-mafL))
  f0<-f1/(exp(lor)*(1-f1)+f1)
  gen <- matrix(nrow=Nc+Nn, ncol=KKmajorAllele)
  for (m in (1:KKmajorAllele)) {
    gen[1:Nc, m] <- rbinom(Nc, 1, f1[m])+rbinom(Nc, 1, f1[m]);
    gen[(Nc+1):(Nc+Nn), m] <- rbinom(Nc, 1, f0[m])+rbinom(Nc, 1, f0[m])
  }
  for (i in (1:(Nc+Nn))) GWASwithAPOEwithRAREwithPRS[i] <- GWASwithAPOEwithRAREwithPRS[i]+sum(lor*gen[i,])  # polygenic risk score
  
  #RANDOM
    KK<-Nmar[SC]-KK
    lor<-rnorm(KK, Bmean[SC], Bsdev[SC])
    f1<-runif(KK, mafL, mafU)
    f0<-f1
    gen <- matrix(nrow=Nc+Nn, ncol=KK)
    for (m in (1:KK)) {
      gen[1:Nc, m] <- rbinom(Nc, 1, f1[m])+rbinom(Nc, 1, f1[m]);
      gen[(Nc+1):(Nc+Nn), m] <- rbinom(Nc, 1, f0[m])+rbinom(Nc, 1, f0[m])
    }
    for (i in (1:(Nc+Nn))) GWASwithAPOEwithRAREwithPRS[i] <- GWASwithAPOEwithRAREwithPRS[i]+sum(lor*gen[i,])
}#SCENARIO

cdat <- glm(status ~ GWASwithAPOEwithRAREwithPRS, family=binomial)$coefficients
cdat

# calculate means, variances and st.dev. of PRS in cases and in non-cases
m1 <- mean(GWASwithAPOEwithRAREwithPRS[1:Nc])
m0 <- mean(GWASwithAPOEwithRAREwithPRS[(Nc+1):(Nc+Nn)])
v1 <- var(GWASwithAPOEwithRAREwithPRS[1:Nc])
v0 <- var(GWASwithAPOEwithRAREwithPRS[(Nc+1):(Nc+Nn)])
s1 <- sqrt(v1)
s0 <- sqrt(v0)

# hence population data
mp <- K*m1 + (1-K)*m0  # mean PRS in population
vp <- K*s1*s1 + (1-K)*s0*s0 + K*(1-K)*(m1 - m0)*(m1 - m0) # PRS variance in population

# theoretical regression coefficients
r0 <- (s1*s1 + (m0 - m1)*(m0 - m1))/s0/s0
r1 <- (s0*s0 + (m0 - m1)*(m0 - m1))/s1/s1
beta <- (m1 - m0)*(K*(1-K)*((r0+r1)/2 - 1) + K*s1*s1/s0/s0 + (1-K)*s0*s0/s1/s1)/vp
alpha <- log(K*s0/(1-K)/s1) + ((r0-1)*K + (1-r1)*(1-K))/2 - mp*beta
c(alpha, beta) # for comparison

grid <- (-40:40)/10
X <- sqrt(vp)*grid+mp
Y <- 1/(1 + exp(-alpha - beta*X))

points(Y ~ grid, type="l", ylim=c(0,1), col="black")
Ydat <- 1/(1 + exp(-cdat[1] - cdat[2]*X))
points (Ydat ~ grid, col="black", pch=4, cex = 0.5)


legend(x="topleft", 
       legend=c("GWAS SNPs without APOE (formula)", "GWAS SNPs without APOE (regression)",
                "APOE and GWAS SNPs (formula)", "APOE and GWAS SNPs (regression)", 
                "Full PRS (formula)", "Full PRS (regression)"), 
        col=c("purple", "purple", "red", "red", "black", "black"), 
        lty=c(1, 0, 1, 0, 1, 0), 
        pch=c(NA,4, NA,4, NA,4),
cex=0.8)


