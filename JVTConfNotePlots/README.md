JVTConfNotePlots
=======================

Instructions on how to reproduce the plots from the JVT Conf Note. 

In [macros](macros) you find the python macros used to make the various plots. The macros are named after the actual plot from the Conf note they produce. I.e. the macro that produces Fig 1 from the note is called Fig01.py

To run this, you need to have the following setup:
 * reasonably recent root version with pyroot support, e.g. 5.34/05 or later
 * access to /atlas/output/pnef on atlint (e.g. atlint03.slac.stanford.edu)

Here's an example:
```
python Fig01a.py
```
produces the figure 
![Fig 1a] (https://raw.githubusercontent.com/pnef/JVTCodeAndDocumentation/master/JVTConfNotePlots/plots/fig_01a.png)
