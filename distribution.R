#############################################
##### First derivitive of link function #####
#############################################
dggaussian = function(mu) rep(1,length(mu))

dgpoisson = function(mu) 1/mu

dgbinomial = function(mu,n.size) (n.size)/(mu*(n.size-mu))

##############################################
###### Inverse funtion of link function ######
##############################################
MUgaussian = function(theta) theta

MUpoisson = function(theta) exp(theta)

MUbinomial = function(theta,n.size)   n.size*exp(theta)/(1+exp(theta))

###################################
##### Log likelihood function #####
###################################
loglikgaussian = function(x,theta) sum(na.omit(as.vector(-(theta^2-2*theta*x)/2)))

loglikpoisson = function(x,theta) sum(na.omit(as.vector(x*theta-exp(theta))))

loglikbinomial = function(x,theta,N) sum(na.omit(as.vector(x*theta-N*log(1+exp(theta)))))
