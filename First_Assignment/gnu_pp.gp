set terminal postscript color
set output "histogram.ps"
set samples 1000


# t (Student's t) PDF
nu = 3.0
mu = 1.0 
sigma = 1.0
t(x) = gamma((nu+1)/2)/(gamma(nu/2)*sqrt(nu*pi)/sigma)*(1+((x-mu)/sigma)**2/nu)**(-(nu+1)/2)

# plot the target distribution


# plot the histogram of sampled data
bin_width = 0.2
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width*(bin_number(x)+0.5)

stats 'dummy.txt' u 1

rand_walk_sig(x) = 1/sqrt(t(x)*1000)

plot [-5:6] t(x) linetype 4 linecolor "blue" linewidth 2 notitle,\
     'dummy.txt' using (rounded($1)):(1)/(bin_width*STATS_records) \
     	 smooth frequency with boxes linecolor "black" linewidth 2 notitle,\
     '+' using 1:(t(x)*(1+3*rand_walk_sig(x))):(t(x)*(1-3*rand_walk_sig(x))) with filledcurves closed

     	 



