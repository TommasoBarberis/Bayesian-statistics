# 1/ installer JAGS à partir du lien: http://sourceforge.net/projects/mcmc-jags/

# 2/ installer le paquet rjags depuis R 
#   (sous Windows, si le test ne fonctionne pas et que vous avez installé la version 4.11 de rjags et 4 ou sup de JAGS, c'est probablément à cause de votre version de R, desintallés R et installés une version inférieure à la v.4)

# 3/ lancer le script suivant pour tester si ça fonctionne

# CHARGEMENT DU PACKAGE
library(rjags)

# MODELE
model <- 
  "model
{
  for (i in 1:Ndose)
  {
  pmal[i] <- 1 - (1 - r)^dose[i]
  Nmal[i] ~ dbin(pmal[i], N[i])
  }
  log10r ~ dunif(-15, -5)
  r <- 10^log10r
}
"

# DONNEES
data4jags <- list(dose = c(8e7, 4e8, 1e9, 2.5e9, 1.2e10, 1.6e11, 3.2e11),
                  N = rep(10,7),
                  Nmal = c(0, 0, 3, 2, 8, 9, 10),
                  Ndose = 7)

# VALEURS INITIALES
ini <- list(list(log10r = -12), 
            list(log10r = -11), 
            list(log10r = -10)) 

# MCMC
m <- jags.model(file=textConnection(model),
                data = data4jags,
                n.chains = 3,
                inits = ini)
update(m, 3000)
mc <- coda.samples(m, c("r"), n.iter = 5000)

plot(mc)
gelman.diag(mc)
autocorr.plot(mc[[1]])
summary(mc)


