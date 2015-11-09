09:10 2015-10-30
*ran it overnight:
**brute force with 25 points
**followed by fmin
*unfortunately it didnt work :C that is, it never converged to SSE<1 even with long run time smooth curves
*going to try recursive brute force instead now
**brute force search returns a winning point
**narrow the brute force search to a smaller grid centered around the winning point
**repeat
