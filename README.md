# minuit1

- expFit.cpp: example using the more general fitting interface in minuit
- expFit.ipynb: equivalent example using lmfit
- SimultaneousExps(lm).ipynb: generation of histograms with correlated signals for simultaneous fit exercise
- rootExample.cpp: just another example of using ROOT classes in a C++ program
- *.root files: various input histograms for fitting exercises

-----

Team member names and computing IDs: Sophia Spaner (VPQ5KD)

-----
**Disclaimer:** Rather than put a chatGPT credit in each of my files, I am putting one in the README. Like other assignments, I used chatGPT to aid in both the development of this work and in my understanding of this assignment. I am able to explain any part of my code if needed.

Exercise 2 comments:
--
Extracted Mean = 54.742 +/- 3.169
Extracted Sigma = 22.63 +/- 1.818
Chi^2 Value = 169.854
P Value = 1.454 * 10^-6
Reduced Chi^2  = 1.846

Though my P and reduced Chi^2 values appear to represent statistically significant data, my fits do not appear to appropriately highlight the signal. For this reason, though I may be interpreting this data incorrectly, I do not believe that I achieved a good fit for the data.

-----

Exercise 3 comments:
--
Parameters:
Amplitude = 54.1704 +/- 0.519731
Mean x = 3.51695 +/- .00541317
Sigma x = 0.989424 +/- 0.00705671
Mean y = 1.89914 +/- 0.0150996
Sigma y =1.95374 +/- 0.0177029
Background Scale Value = 0.242787 +/- 0.00227717

Signal Events:
Number of Signal Events = 658 +/- 10

Explanation: If we consider each signal event to make up a portion of the signal, then the number of signal events is simply the volume under the gaussian signal. This is defined, per the wikipedia article on Gaussian functions, to be 2 * pi * Amplitude * sigma_x * sigma_y. For the error I simply took the error of each of these variable quantities and used the appropriate error propagation formula to come to a result. 


-----
