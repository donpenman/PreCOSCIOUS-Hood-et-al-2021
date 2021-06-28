# PreCOSCIOUS-Hood-et-al-2021
MATLAB code for the 3 runs presented in the main text of Hood et al 2021 ("Cryogenian syn-glacial carbonate precipitation and implications for a snowball Earth")

Modified from the PreCOSCIOUS model in Penman & Rooney (2019, Geology)

Each will run a snowball earth simulation, with different parameters for deglacation threshold and silicate weathering rates. 

All runs require the file "SVSTRT.txt" in the run directory.

Higher Deglaciation Threshold: Same as the P&R19 run, but threshold to trigger deglaciation is higher: 0.8bar, rather than the standard 0.12bar
9% weathering feedback: Same as P&R19 run, but during the snowball interval silicate weathering operates as a function of pCO2 but is scaled down to 9% of its default (non-glacial) magnitude. 
80% Fvc weathering: Same as P&R19 run, but during the snowball interval, silicate weathering is set constant to value equal to 80% of the background volanic degassing rate. 
