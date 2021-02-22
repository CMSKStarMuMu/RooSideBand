# RooSideBand

Table of contents:

* testSidebandFit       -> code to test RooBernsteinSideband models (pdf): 
                           
			   fits&plots, re-compute the p-value, pulls studies...
			   
		           All the files with sideband models will be provided.
			   
		           However, if you need to save again sb models in rooworkspace wsb 
			   (in file savesb_*root), compile the code and run: 
		   
		   			./testSidebandFit [Q2Bin] [Era]
					
		           [ testSidebandFit reads: 
			     * file namelist-SB3DB0-[Era]-[Q2Bin]-3.lis 
			     * ntuples on EOS
			     * coefficients saved in ListParValues-[Era]--[Q2Bin]-*-SigmaProb.txt ]	
        
* readSideband           -> simple program to read&plots the bkg pdfs (angular&mass)

* namelist*.lis         -> lists of command parameters for testSidebandFit [degrees of Bernstein polynomial, bins for plots...]

* ListParValues*.txt    -> lists of coefficients for bkg models (read by testSidebandFit)

* Makefile              -> compile 
