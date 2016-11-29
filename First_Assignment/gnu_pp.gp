set terminal postscript
set output "histogram.ps"
set samples 1000

# t (Student's t) CDF
nu = 3.0
mu = 1.0 
sigma = 1.0
t(x) = gamma((nu+1)/2)/(gamma(nu/2)*sqrt(nu*pi)/sigma)*(1+((x-mu)/sigma)**2/nu)**(-(nu+1)/2)

set title "Student t MCMC Samples"

plot [-5:6] t(x)


