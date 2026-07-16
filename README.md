# GHATS Timing Software

GHATS General High-Energy Aperiodic Timing Software originally written by Tomaso M. Belloni (INAF-OAB)

More than 20 years ago, Tomaso Belloni (INAF-OAB) started writing a package devoted to RXTE/PCA data. The philosophy behind the software was:

1. The program had to work directly from RXTE raw files, without forcing the user to produce light curves before the analysis. In other words, the program itself would accumulate its own light curves, which would be discarded at the end.
2. The analysis could be divided into two separate parts. (i) The production of Power Density Spectra (PDS) or Fast Fourier Transforms (FFT). (ii) The analysis of those products. 

The first task is essentially an almost non-interactive process: one gives input parameters and the program delivers the final product (there is some interaction as the range of possible parameters must depend on the specific data types and a decision can only be made after reading the input files). The second task is highly interactive: read PDS and FFT files, manipulate them, fit them, etc.

Given these starting points, Tomaso decided to write the software in an integrated (graphical) environment. His choice fell upon IDL because it was a very widely used package and there were existing libraries which made life easier. (At the time, Python, with all the currently existing libraries to do this kind of analysis, was not yet widely available.) The downside was that IDL was not a freely available program (this is still the case), however, users have access to a public-domain clone of IDL, the GNU Data Language (GDL). 

From all this, a package called MU was born (the name comes from “Merate” and “Utrecht” since, besides Tomaso being in Merate, Mariano Méndez, at the time in Utrecht, was tangentially involved in the project). Initially only used within Tomaso’s institute, it was moderately exported so that it is was used by people in Groningen, Rome, Mumbai and Kolkata with whom Tomaso shared the code. 

During several years, Sara Motta and Mariano Méndez contributed to changes in the software; the official version was always under Tomaso’s control. In the last 5 years Federico García got more and more involved, and in Groningen Federico and Mariano developed several routines aimed at analysing the cross spectra. From all this, version 3.3.0 was born.

To recognise Tomaso’s work and honour his legacy, we decided to upload the last version of the code that we received from him, which he called 3.2.0, and immediately after that uploaded the latest version, 3.3.0. In doing so, we also updated the Manual to reflect the latest additions.

Originally, the GHATS package could only be obtained directly from Tomaso, and was not distributed in any way. This was so because he did not have the resources to support the software for a broad community. From 3.3.0 we decided to make the code available on this Github repository to the community. We, however, cannot support all the request for updates, upgrades, etc., and will allow the community to add contributed routines or modify parts of the code (e.g., to correct bugs) on a development branch.

We miss Tomaso, and his enthusiasm in adding stuff to the code whenever we came with ideas. We would have very much liked to share this new step with him. In the meantime, we try to honour his memory by trying to do good work using this software. We are sure that you, the new generation of users, will do the same.
