09:53 2015-10-29
*encaounted a strange error with low values of k_off
*appears to be sharp drops at the transition points between runs
*checking into it using an older version of the 1 site which also exhibits the bug

11:25 2015-10-29
*the error is corrected
*the problem was due to cumulative walking times exceeding the injection time points and not being set to the timepoint
**instead holding times were added on, so extremely long holding times would cause cumulative times to have extremely long entries rather than being reset at the injection points
*fitting is still not working
**i think it is due to noise in the monte carlo
**it is taking too long to increase the iterations/sites to reduce noise
**i think instead, i should attempt to adjust the minimization algorithm parameters to take bigger steps and be less sensitive to small changes

12:55 2015-10-29
*manual fitting seems to be working though...