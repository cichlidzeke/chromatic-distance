# chromatic-distance
This repository contains R scripts for three vision-related calculations, used specifically for east African cichlid fish species. One, calculating visual pigments absorption using Govardovskii et al. 2000 templates. Two, calculating the stimulation of cone photoreceptors in terms of quantum catch. And three, calculating JND values between two sets of targets
 
This project utilizes R.
 
## Contents:

VisualFxns.R

ScriptEx.R
 
## Suggested Usage:
 
1. Use the Opsin_gov function from VisualFxns.R to calculate visual pigment absorption. 
 
2. Use the Qcatch_calc function from VisualFxns.R to calculate the stimulation of the S, M, and L cones in terms of quantum catch.

3. Calculate your Weber fractions for the relevant cone photoreceptors of your species. 

4. Use the JND_calc function from VisualFxns.R to calculate the JND values between two sets of color targets.

5. Manipulate the resulting table of JND values to remove duplicates, stack the data, and name the data columns.

6. Calculate your summary statistics, visualize the JNDs, or otherwise process your resultant chromatic distances.


**See ScriptEx.R for an example of this process for the species *Metriaclima zebra*.**
