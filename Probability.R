#Function
coeff <- function(m0, m1, v0, v1, K){
  r1<-(v0+(m1-m0)^2)/v1
  r0<-(v1+(m1-m0)^2)/v0
  mp<-K*m1+(1-K)*m0
  vp<-K*v1+(1-K)*v0+K*(1-K)*(m1-m0)^2
  beta<-(m1-m0)*(K*(1-K)*((r0+r1)/2-1) + K*v1/v0 + (1-K)*v0/v1)/vp
  alpha<-log(K*sqrt(v0/v1)/(1-K)) + ((r0-1)*K + (1-r1)*(1-K))/2 - mp*beta
  return(c(alpha, beta))}
grid<-(-40:40)/10

Kp<-0.02 # AD prevalence in the population
K_EO<-0.1 #AD prevalence in 65+
K_LO<-0.3 #AD prevalence in85+

#PRS mean and variance as reported in Leonenko et al (2021), Nat Comms
PRS.ADm1<- 0.630
PRS.ADv1<- 1.200^2
PRS.ADm0<- -0.370
PRS.ADv0<- 0.946^2

#RARE###########################################################################
OR_SORL1<-7.22 #SORL1 early onset

#mean and variance in population
mpp<-Kp*PRS.ADm1+(1-Kp)*PRS.ADm0
vpp<-Kp*PRS.ADv1+(1-Kp)*PRS.ADv0+Kp*(1-Kp)*(PRS.ADm1-PRS.ADm0)^2
x<- mpp + sqrt(vpp)*grid

#SORL1
OR_Rare<-OR_SORL1
p_Rarep <-Kp*(OR_Rare-1)/(Kp*(OR_Rare-1)+1)
p_RareEO<-K_EO*(OR_Rare-1)/(K_EO*(OR_Rare-1)+1)

co<-coeff(PRS.ADm0, PRS.ADm1, PRS.ADv0, PRS.ADv1, Kp)
alpha<-co[1]
beta<-co[2]
y<-1/(1+exp(-alpha-beta*x))
plot(y ~ grid, ylim=c(0, 1), type="l", col="black", lwd=2, lty=3, 
     xlab="PRS", ylab="Probability", 
     main="Probability of Alzheimers disese\n based on PRS and PRS+SORL1")
y_Rare<-y+p_Rarep*(1-y)
points(y_Rare ~ grid, type="l", col="black", lwd=3)

#ADJUSTING for prevalence in Early onset sporadic AD
co<-coeff(PRS.ADm0, PRS.ADm1, PRS.ADv0, PRS.ADv1, K_EO)
alpha<-co[1]
beta<-co[2]
y<-1/(1+exp(-alpha-beta*x))
points(y ~ grid, type="l", col="red", lwd=2, lty=3)
y_Rare<-y+p_RareEO*(1-y)
points(y_Rare ~ grid, type="l", col="red", lwd=3)

legend("topleft", 
       legend=c("General population (K=0.02): PRS", 
                "General population (K=0.02): PRS+SORL1", 
                "Population of age 65+ (K=0.1): PRS", 
                "Population of age 65+ (K=0.1): PRS+SORL1"),
       col=c("black", "black",  "red", "red"), lty=c(3,1,3,1), 
       cex=0.8, lwd=c(2,3,2,3))






#APOE##################################################################################
Kp<-0.02 #AD prevalence in the population
K <- 0.1 #AD prevalence in age group: K
f <- 0.2  #allele 1 frequency in age group: f,
X <- 3.064854 #APOE OR


#PRS.noAPOE as iin Leonenko et al (2021), Nat Comm
PRS.noAPOEm1<- 0.320
PRS.noAPOEv1<- 1.089^2
PRS.noAPOEm0<- -0.136
PRS.noAPOEv0<- 1.070^2

mb <- 2*(1+(1-X)*(K-f))
nn0 <- (mb - sqrt(mb*mb-16*(1-X)*(1-f)*K))/(2*(1-X))
nn1 <- 2*K - nn0
mm0 <- 2*(1-f) - nn0
mm1 <- 2*f - nn1
nn1*mm0/nn0/mm1 # to check OR

# table allele 0, allele 1 in cases (row 1) and noncases (row 2)
matrix(c(nn0, mm0, nn1, mm1), nrow=2)

# getting genotype fractions as a function of n1 (heterozygotes in cases)
n0 <- function(n1){return((nn0-n1)/2)}
n2 <- function(n1){return((nn1-n1)/2)}
m0 <- function(n1){return((mm0-2*f*(1-f)+n1)/2)}
m1 <- function(n1){return(2*f*(1-f)-n1)}
m2 <- function(n1){return((mm1-2*f*(1-f)+n1)/2)}
# table of genotype 0, 1, 2 fractions in cases (row 1), noncases (row 2)
tabl <- function(n1){return(matrix(c(n0(n1),m0(n1),n1,m1(n1),n2(n1),m2(n1)), nrow=2))}

# Now chose n1 to fit HWE in noncases (which may make sense if K is small)
tabl(2*f*(1-f) - mm0*mm1/(mm0+mm1))

# So we get the prevalence in 00 homozygotes and in risk allele carriers
n1 <- 2*f*(1-f) - mm0*mm1/(mm0+mm1)
K0 <- n0(n1)/(n0(n1)+m0(n1))
K1 <- n1/(n1+m1(n1))
K2 <- n2(n1)/(n2(n1)+m2(n1))
c(K0, K1, K2) # print prevalences for genotype 0, 1, 2 in age group

#mean and variance in population
mpp<-Kp*PRS.noAPOEm1+(1-Kp)*PRS.noAPOEm0
vpp<-Kp*PRS.noAPOEv1+(1-Kp)*PRS.noAPOEv0+Kp*(1-Kp)*(PRS.noAPOEm1-PRS.noAPOEm0)^2
x<- mpp + sqrt(vpp)*grid

co<-coeff(PRS.noAPOEm0, PRS.noAPOEm1, PRS.noAPOEv0, PRS.noAPOEv1, K)
alpha<-co[1]
beta<-co[2]
y<-1/(1+exp(-alpha-beta*x))

plot(y ~ grid, ylim=c(0, 1), type="l", col="black", xlab="PRS",lwd=2, lty=2, 
     ylab="Probability", main="Probability of Alzheimers disease with/without APOE\nEarly onset (K=0.1)")
points(y ~ grid, ylim=c(0, 1), type="l", col="black", xlab="PRS",lwd=2, lty=2)


co<-coeff(PRS.noAPOEm0, PRS.noAPOEm1, PRS.noAPOEv0, PRS.noAPOEv1, K0)
alpha<-co[1]
beta<-co[2]
y<-1/(1+exp(-alpha-beta*x))
points(y ~ grid, type="l", col="blue", lwd=2, lty=1)

co<-coeff(PRS.noAPOEm0, PRS.noAPOEm1, PRS.noAPOEv0, PRS.noAPOEv1, K1)
alpha<-co[1]
beta<-co[2]
y<-1/(1+exp(-alpha-beta*x))
points(y ~ grid, type="l", col="blue", lwd=3, lty=1)

co<-coeff(PRS.noAPOEm0, PRS.noAPOEm1, PRS.noAPOEv0, PRS.noAPOEv1, K2)
alpha<-co[1]
beta<-co[2]
y<-1/(1+exp(-alpha-beta*x))
points(y ~ grid, type="l", col="blue", lwd=4, lty=1)


legend("topleft", 
       legend=c(
#         "PRS.AD",
         "PRS.no.APOE",
         "e4 non-carriers", 
         "e4 carriers",
         "e4e4 carriers"
       ),
       col=c("black", "blue", "blue", "blue"), 
       lty=c(2,1,1,1), 
       lwd=c(2,2,3,4))




# AD prevalence in age group K, allele 1 frequency in age group f, OR X
K_LO <- 0.3
f <- 0.05
X <- 3.064854

K<-K_LO
mb <- 2*(1+(1-X)*(K-f))
nn0 <- (mb - sqrt(mb*mb-16*(1-X)*(1-f)*K))/(2*(1-X))
nn1 <- 2*K - nn0
mm0 <- 2*(1-f) - nn0
mm1 <- 2*f - nn1
nn1*mm0/nn0/mm1 # to check OR

# table allele 0, allele 1 in cases (row 1) and noncases (row 2)
matrix(c(nn0, mm0, nn1, mm1), nrow=2)

# getting genotype fractions as a function of n1 (heterozygotes in cases)
n0 <- function(n1){return((nn0-n1)/2)}
n2 <- function(n1){return((nn1-n1)/2)}
m0 <- function(n1){return((mm0-2*f*(1-f)+n1)/2)}
m1 <- function(n1){return(2*f*(1-f)-n1)}
m2 <- function(n1){return((mm1-2*f*(1-f)+n1)/2)}
# table of genotype 0, 1, 2 fractions in cases (row 1), noncases (row 2)
tabl <- function(n1){return(matrix(c(n0(n1),m0(n1),n1,m1(n1),n2(n1),m2(n1)), nrow=2))}

# Now chose n1 to fit HWE in noncases (which may make sense if K is small)
tabl(2*f*(1-f) - mm0*mm1/(mm0+mm1))

# So we get the prevalence in 00 homozygotes and in risk allele carriers
n1 <- 2*f*(1-f) - mm0*mm1/(mm0+mm1)
K0 <- n0(n1)/(n0(n1)+m0(n1))
K1 <- n1/(n1+m1(n1))
K2 <- n2(n1)/(n2(n1)+m2(n1))
c(K0, K1, K2) # print prevalences for genotype 0, 1, 2 in age group

m1<- PRS.noAPOEm1
v1<- PRS.noAPOEv1
m0<- PRS.noAPOEm0
v0<- PRS.noAPOEv0

mpp<-Kp*m1+(1-Kp)*m0
vpp<-Kp*v1+(1-Kp)*v0+Kp*(1-Kp)*(m1-m0)^2
x<- mpp + sqrt(vpp)*grid
co<-coeff(m0, m1, v0, v1, K_LO)
alpha<-co[1]
beta<-co[2]
y<-1/(1+exp(-alpha-beta*x))

plot(y ~ grid, ylim=c(0, 1), type="l", col="black", xlab="PRS",lwd=2, lty=2, 
     ylab="Probability", main="Probability of Alzheimers disease with/without APOE\nLate onset (K=0.3)")

co<-coeff(m0, m1, v0, v1, K0)
alpha<-co[1]
beta<-co[2]
y<-1/(1+exp(-alpha-beta*x))
points(y ~ grid, type="l", col="blue", lwd=2, lty=1)

co<-coeff(m0, m1, v0, v1, K1)
alpha<-co[1]
beta<-co[2]
y<-1/(1+exp(-alpha-beta*x))
points(y ~ grid, type="l", col="blue", lwd=3, lty=1)

co<-coeff(m0, m1, v0, v1, K2)
alpha<-co[1]
beta<-co[2]
y<-1/(1+exp(-alpha-beta*x))
points(y ~ grid, type="l", col="blue", lwd=4, lty=1)

legend("topleft", 
       legend=c(
         "PRS.no.APOE",
         "e4 non-carriers", 
         "e4 carriers",
         "e4e4 carriers"
       ),
       col=c("black", "blue", "blue", "blue"), 
       lty=c(2,1,1,1), 
       lwd=c(2,2,3,4))
