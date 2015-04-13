from pylab       import *
from scipy.stats import scoreatpercentile, distributions

poisson = distributions.poisson.pmf
chi2cdf = distributions.chi2.cdf

def iqr(arr):
  arr           = sort(arr.copy())
  upperQuartile = scoreatpercentile(arr,.75)
  lowerQuartile = scoreatpercentile(arr,.25)
  iqr           = upperQuartile - lowerQuartile
  return iqr

num      = arange(7)
tapeworm = hstack((0 * zeros(235),
                   1 * ones(168),
                   2 * ones(75),
                   3 * ones(32),
                   4 * ones(7), 
                   5 * ones(2),
                   6 * ones(1)))
mu       = mean(tapeworm)                # mean
med      = median(tapeworm)              # median
sigma    = std(tapeworm)                 # standard deviation
fe_iqr   = iqr(tapeworm)                 # IQR
v_m_rat  = sigma**2 / mu                 # variance-to-mean ratio
px       = poisson(range(7), mu)         # compute probabilities
expfreq  = px * len(tapeworm)            # computes expected frequencies

#===============================================================================
# plotting :
fig      = figure()
ax       = fig.add_subplot(111)

ax.hist(tapeworm, num, histtype='stepfilled')                      
ax.plot(num + .5, expfreq, 'gs', label='expected freq')
ax.set_xlabel('# Tapeworms per perch')
ax.set_ylabel('Frequency')
ax.set_title('Histogram of Tapeworm Counts')
ax.legend(loc='center right')
ax.grid()
textstr =   '$\mu = %.2f$\n' \
          + '$\mathrm{median} = %.2f$\n' \
          + '$\sigma = %.2f$\n' \
          + '$\mathrm{IQR} = %.3f$\n' \
          + '$\sigma^2 / \mu = %.2f$'
textstr = textstr % (mu, med, sigma, fe_iqr, v_m_rat)

# these are matplotlib.patch.Patch properties
props   = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

# place a text box in upper left in axes coords
ax.text(0.65, 0.95, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
show()

#===============================================================================
# perform chi^2 test :
mle             = mean(tapeworm)
px              = poisson(range(4), mle)             # for results <= 3
px              = append(px, 1-sum(px))              # add on results > 3
expfreq         = px * len(tapeworm)                 # expected frequency
obsfreq, n, mid = hist(tapeworm, [0,1,2,3,4,6])      # histogram
D               = sum((obsfreq-expfreq)**2/expfreq)  # chi squared
m               = len(n) - 1                         # number of outcomes
pval            = 1 - chi2cdf(D, m-2)                # probability not Poissson
print pval


