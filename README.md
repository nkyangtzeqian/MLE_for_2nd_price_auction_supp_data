# MLE_for_2nd_price_auction_supp_data
Supplement codes and data for paper _A non-parametric approach for estimating consumer valuation distributions using second price auctions_.

## The structure of folders:

**MLE_2ndprice_code**: contain all R files and data.

  &nbsp;&nbsp;&nbsp;&nbsp; **SIMU**: all simulation results in RData.

  &nbsp;&nbsp;&nbsp;&nbsp; **Plots**: all plots, R file and relevant RData.
  
  &nbsp;&nbsp;&nbsp;&nbsp; **Empirical**: all results about XBox dataset on Ebay.
  

## Explanation of R files:

### Under MLE_2ndprice_code folder:
  **DataGeneration.R**: functions to generate either auctions with all unobserved bids or (seal-price) 2nd-price auctions with only observed bidds, in format used for MLE, Init and PT estimate. Also contain functions to slice dataset.\
                    &nbsp;&nbsp;&nbsp;&nbsp;It is also used to generate results in Empirical study to calculate MSE.
                    
  **Main.R**: definition of MLE and initial estimate function, and important function used by them (e.g. step-wise update function). Also include F_fp separately, and the 2-step MLE used in Supplement _Additional simulation experiments: Settings with high expected number of participants per auction_.
  
  **Initialization_lambda_F_Simp.R**: definition of the simplified Init estimate function only used in Supplement _Additional simulation experiments: Settings with high expected number of participants per auction_.
  
  **PT_MCMC.R**: Polya Tree (PT) estimate function (from Prof. Sam K. Hui), as well as an interface function to accommodate the format we used.


  **distributions.R**: Contain some distributions defined by us, include piecewise uniform CDF, and the important log-likelihood_PA function of 2nd-price auction model in Formula (3.3)
  
  **SettingGeneration.R**: auxiliary file, defining all settings in Simulations, e.g. the parameters of F_0 and reserve prices. 
  
  **DisCompare.R**: auxiliary file, include package and R function to compare the KS divergence and TV distance, as well as linear interpolation cdf used in AbscontDistribution class.


  **HulC_main.R**: main functions to calculate the HulC interval for MLE and Init.


  **Simu1.R**: main function to conduct (one of) the simulation in Table 2 (comparison of MLE, Init and PT estimate on KS and TV, with K=100 and 1000; originally Table 1).
  
  **Simu1_Table.R**: automatically conducts all simulations correspond to the previous file with different settings and present Table 2.


  **Simu2_HulC_Coverage.R**:  main function to conduct (one of) the simulation in Table 3 (comparison of HulC band of MLE, Init and credible band of PT on band area and coverage, with K=100 and 1000; originally Table2).

  **Simu2_2_HulC_TruncCoverage.R**:  main function to conduct (one of) the simulation in Table 3 (comparison of HulC band of MLE, Init and credible band of PT on truncated coverage, with K=100 and 1000; originally Table2).
  
  **Simu2_Table.R**: automatically conducts all simulations correspond to the previous two files with different settings and present Table 3. First run coverage, then run truncated coverage.


  **Simu1_2_morebidder.R**:  main function to conduct (one of) the simulation in Table S.1 in Supplement _Additional simulation experiments: Settings with high expected number of participants per auction_ (comparison of MLE, Init and PT estimate on KS and TV, with lambda=10 and 50).
  
  **Simu1_2_Table.R**: automatically conducts all simulations correspond to the previous file with different settings and present Table S.1.

  ### Under Plots folder (and their plots):
  **Instability_Boundary_Correction.R**: 2 functions, one shows the instability on boundary due to the non-identifiability when theta_1 is not fixed; another shows the slow convergence when lambda*tau is large (thus we need the two-step MLE).\
      &nbsp;&nbsp;&nbsp;&nbsp; Results from it is saved in file Example_MLE_2step_for_large_lambda.RData.
  
  **example_plot.R**: Used to generate the Fig 1 and Fig 2 in the paper.\
       &nbsp;&nbsp;&nbsp;&nbsp; Fig 2 is saved as example_plot.pdf.

  **HulC_plot.R**: main function to draw the fitted curves of MLE, Init and PT estimates and their 90% HulC band/credible band in (one of) the settings. 
  
  **HulC_plot_simu_all.R**: automatically conducts all plottings correspond to the previous file with different settings and generate the plots.\
    &nbsp;&nbsp;&nbsp;&nbsp; Results from them is saved as HulC_Simu_unif.pdf and the others in the same folder.

  **median_bias_simu.R**: main function simulate the median bias and difference in median of MLE, Init and PT (which is not meaningful) estimate to the true F0 with K=20 (since we have batches in HulC methods; when K grows, they will be even smaller) in (one of) the settings. Could be seen that the median bias is usually small for MLE excluding the boundary for all settings we used in the plots. Use as initial study, only for reference.

### Under Empirical folder (and their plots):

  **XBOX_Pre.R**: extract organized data from a real (raw) Xbox Auctions data set and covert to our format.
  
  **MSE_selling_price.R**: function to fit models on train set, and generate (simulate) the final selling price on test set and compare them to the true value in MSE. Corresponding to Table 4. Note the PT here uses estimated result from MLE for the unknown total number of bidders, only for reference.
  
  **HulC_XBOX.R**: function to generate fitted CDF curves and bands in the whole dataset.\
    &nbsp;&nbsp;&nbsp;&nbsp; Results saved in HulC_XBOX_all.pdf.

    
  **check_seq_auc.R**: function to check the multiple-bidding issue.
  
  **dis_selling_price.R**: comparison of the final selling price generated by each method, and calculate The Average Overlap Ratio and the Average Range Ratio for each methods. Not used in the paper, just for reference. Further investigate shows heterogeneity in some of the XBOX products which cannot be explained by observed bid only.


## Explanation of R Data file

### Under SIMU folder:
  **Table1_unif_K=100.RData** and so on: results of Simu1.R and Simu1_Table.R.
  
  **HulC_Cov_unif With No. Auction = 100.RData** and so on: results of Simu2_HulC_Coverage.R and Simu2_Table.R.

  **HulC_TrunCov_unif With No. Auction = 100.RData** and so on: results of Simu2_2_HulC_TruncCoverage.R and Simu2_Table.R (truncated coverage).
  
  **Table1_unif_Lbd=10.RData** and so on: results of Simu1_2_morebidder.R and Simu1_2_Table.R.

### Under plots folder:
  **Example_MLE_2step_for_large_lambda.RData**: results of Instability_Boundary_Correction.R.

### Under Empirical folder:
  **Xbox_7day_auctions.csv** Original dataset after the tie-breaking procedure.
  
  **Xbox_7day_auctions.RData**: dataset we used, produced by XBOX_Pre.R

  **Xbox_simu_0.5.RData** and **Xbox_simu_0.67.RData**: results of MSE_selling_price.R.
  
  **Xbox_simu_dist.RData**: result of dis_selling_price.R.
  

