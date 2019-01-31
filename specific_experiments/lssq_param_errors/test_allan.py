

import numpy as np 
import allantools
import allan_dataset as ad
import allan_plot as ap


# Compute a deviation using the Dataset class
# a = ad.Dataset(data=np.random.rand(1000))
# a.compute("mdev")

# # Plot it using the Plot class
# b = ap.Plot()
# b.plot(a, errorbars=True, grid=True)
# # You can override defaults before "show" if needed
# b.ax.set_xlabel("Tau (s)")
# b.savefig('testallan.eps');

data=np.random.rand(1000);

taus_used, ad, [ade_l, ade_h], adn = allantools.gradev(data);

print("\n");
print(ad);

    # taus: np.array
    #     list of tau vales in seconds
    # adev: np.array
    #     deviations
    # [err_l, err_h] : list of len()==2, np.array
    #     the upper and lower bounds of the confidence interval taken as
    #     distances from the the estimated two sample variance.
    # ns: np.array
    #     numper of terms n in the adev estimate.