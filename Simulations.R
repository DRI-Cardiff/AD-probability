grid <- (-40:40)/10

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
for (i in (1:(Nc+Nn))) PRS[i] <- sum(lor*gen[i,]) # polygenic risk score
sm<-mean(PRS)
ssd<-sd(PRS)
sPRS<-(PRS-sm)/ssd
cdat <- glm(status ~ sPRS, family=binomial)$coefficients

# calculate means, variances and st.dev. of PRS in cases and in non-cases
m1 <- mean(sPRS[1:Nc])
m0 <- mean(sPRS[(Nc+1):(Nc+Nn)])
v1 <- var(sPRS[1:Nc])
v0 <- var(sPRS[(Nc+1):(Nc+Nn)])
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

X <- sqrt(vp)*grid+mp
Y <- 1/(1 + exp(-alpha - beta*X))

plot(Y ~ grid, type="l", ylim=c(0,1), col="red", 
     main="Probability of disease", ylab="Probability", xlab="Standardised PRS")
Ydat <- 1/(1 + exp(-cdat[1] - cdat[2]*X))
points (Ydat ~ grid, col="red", pch=4, cex = 0.5)

GWASwithAPOE<-PRS
PRS_noAPOE<-PRS
for (i in (1:(Nc+Nn))) PRS_noAPOE[i] <- PRS[i]-lor[37]*gen[i,37]  #
sm<-mean(PRS_noAPOE)
ssd<-sd(PRS_noAPOE)
sPRS_noAPOE<-(PRS_noAPOE-sm)/ssd
cdat <- glm(status ~ sPRS_noAPOE, family=binomial)$coefficients

# calculate means, variances and st.dev. of PRS in cases and in non-cases
m1 <- mean(sPRS_noAPOE[1:Nc])
m0 <- mean(sPRS_noAPOE[(Nc+1):(Nc+Nn)])
v1 <- var(sPRS_noAPOE[1:Nc])
v0 <- var(sPRS_noAPOE[(Nc+1):(Nc+Nn)])
s1 <- sqrt(v1)
s0 <- sqrt(v0)

mp <- K*m1 + (1-K)*m0  # mean PRS in population
vp <- K*s1*s1 + (1-K)*s0*s0 + K*(1-K)*(m1 - m0)*(m1 - m0) # PRS variance in population

# theoretical regression coefficients
r0 <- (s1*s1 + (m0 - m1)*(m0 - m1))/s0/s0
r1 <- (s0*s0 + (m0 - m1)*(m0 - m1))/s1/s1
beta <- (m1 - m0)*(K*(1-K)*((r0+r1)/2 - 1) + K*s1*s1/s0/s0 + (1-K)*s0*s0/s1/s1)/vp
alpha <- log(K*s0/(1-K)/s1) + ((r0-1)*K + (1-r1)*(1-K))/2 - mp*beta

X <- sqrt(vp)*grid+mp
Y <- 1/(1 + exp(-alpha - beta*X))

points(Y ~ grid, type="l", col="purple")
Ydat <- 1/(1 + exp(-cdat[1] - cdat[2]*X))
points (Ydat ~ grid, col="purple", pch=4, cex = 0.5)










#Simulation of additional PRS
dat<-read.table(file="Kunkle_published_noATCG_clump_1000GLDref_r0.1_1000kb_exclAPOE_exclGWAS1MB.txt", header=T)
r<-seq(1:nrow(dat))
a<-sample(x=r, size=10000,  replace=FALSE)
dat<-dat[a,]
Nmar<-nrow(dat); Nmar
GWASwithAPOEwithPRS<-GWASwithAPOE

K<-22000/64000 #as in Kunkle et al 2019
BETAs<-dat$BETA
f <-dat$Effect_allele_freq
f1<-dat$Effect_allele_freq
f0<-dat$Effect_allele_freq
for (i in 1:length(BETAs)) {BETAs[i]<-rnorm(1,dat$BETA[i], dat$SE[i])}
OdR<-exp(BETAs)
for (i in 1:length(f))
{
  a<- K*(OdR[i]-1)
  b<- (f[i]+K)*(1-OdR[i])-1
  c<- f[i]*OdR[i]
  D<-b^2-4*a*c
  f1[i]<- min ( (-b+sqrt(D))/2/a, (-b-sqrt(D))/2/a)
  if (is.na(f1[i])) {f1[i]<-f[i]}
  if (f1[i]<0) {f1[i]<- max ( (-b+sqrt(D))/2/a, (-b-sqrt(D))/2/a)}
  f0[i]<- f1[i]/(OdR[i]*(1-f1[i])+f1[i])
}

lor<-dat$BETA
gen <- matrix(nrow=Nc+Nn, ncol=Nmar)
for (m in (1:Nmar)) {
  gen[1:Nc, m] <- rbinom(Nc, 1, f1[m])+rbinom(Nc, 1, f1[m]);
  gen[(Nc+1):(Nc+Nn), m] <- rbinom(Nc, 1, f0[m])+rbinom(Nc, 1, f0[m])
}
for (i in (1:(Nc+Nn))) GWASwithAPOEwithPRS[i] <- GWASwithAPOEwithPRS[i]+sum(lor*gen[i,])
sm<-mean(GWASwithAPOEwithPRS)
ssd<-sd(GWASwithAPOEwithPRS)
sGWASwithAPOEwithPRS<-(GWASwithAPOEwithPRS-sm)/ssd
model<- glm(status ~ sGWASwithAPOEwithPRS, family=binomial)
cdat <- model$coefficients

# calculate means, variances and st.dev. of PRS in cases and in non-cases
m1 <- mean(sGWASwithAPOEwithPRS[1:Nc])
m0 <- mean(sGWASwithAPOEwithPRS[(Nc+1):(Nc+Nn)])
v1 <- var(sGWASwithAPOEwithPRS[1:Nc])
v0 <- var(sGWASwithAPOEwithPRS[(Nc+1):(Nc+Nn)])
s1 <- sqrt(v1)
s0 <- sqrt(v0)

K<-Nc/(Nc+Nn) #as in the simulatedsample
mp <- K*m1 + (1-K)*m0  # mean PRS in population
vp <- K*s1*s1 + (1-K)*s0*s0 + K*(1-K)*(m1 - m0)*(m1 - m0) # PRS variance in population

# theoretical regression coefficients
r0 <- (s1*s1 + (m0 - m1)*(m0 - m1))/s0/s0
r1 <- (s0*s0 + (m0 - m1)*(m0 - m1))/s1/s1
beta <- (m1 - m0)*(K*(1-K)*((r0+r1)/2 - 1) + K*s1*s1/s0/s0 + (1-K)*s0*s0/s1/s1)/vp
alpha <- log(K*s0/(1-K)/s1) + ((r0-1)*K + (1-r1)*(1-K))/2 - mp*beta
c(alpha, beta) # for comparison

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


