# FinancialEconometrics
Finmetrics Project

This project was initially to examine portfolio optimization. I also wanted to look at GAS modeling in the multivariate space and compare DCC models of different specification to examine the differences and whether any it is possible to outperform the market using portfolio optimization. The GAS modeling section was eventually mostly cut as multivariate GAS takes considerable time to calculate due to parameter explostion in the updating equation. As this is not mentioned elsewhere or at least emphasised as a problem I include the section on how GAS adjust DCC and why it is not possible in the case of portfolio optimization

Regarding the file Write_Up. This is the file that contains the document information, maths graph construction, table construction and the document can be knit to pdf with no problems. It has 4 datasets as inputs namely: the DailyTop40.rds, ModelReturns.RData, ModelReturnsSpread.RDS and LowerleverageModelReturns.RDS.  The DailyTop40.rds is the dataset of all stocks and is used in PortfolioGeneration.R to generate portfolios. PortfoliosGeneration.R is explained at the end

Introduction is left open initially to filled in when models are estimated and basics are written up. Basic explanation of portfolio optimization, where it comes from and what problems exist added.  The reasons for portfolio optimization are explained and why DCC modeling can be used to extend GARCH modeling. Why GARCH is an improvement over traditional variance.  What are benifits of DCC vs orther models and what are some of the problems with DCC. Also touch on GAS modeling

To start with an introduction into portfolio theory and statistics is useful showing how portfolios statistics are constructed helps in the explanation of why Mean Variance optimization works through minimizing the return over variance equation as developed in the Classical Portfolio Statistics Section. 

The next step is obviously to show how the weights can be found to maximize return over variance. This Mean Variance Portfolio is what Markowitz developed. From this some of the requirements for modeling the portfolio is found namely return and variance. Risk free rate was attempted to be incorporated but at certain times the portfolio found in estimation would have global minimum variance below that of the risk free rate which causes errors in estimation. Thus risk free rate is just set to be 0 across the paper. It has minimal effect but might have better sharpe ratios where the sharpe ratio has a given risk free rate.  Also explained what equation needs to be optimized and by what method.

The efficient market hypothesis makes numerous statements about the potential for certain portfolios to outperform. Namely that in the long term no portfolio will outperform the market under perfect information. It is explained why market cap weighting is thus theoretically unbeatable as it is the market portfolio. This is true to some degree however irrationality may lead to differences or market structure. But theoretically unbeatable. This must be tested.

The next section explains the mathematical underpinnings of the various models examined.
The Application of DCC-GARCH Model section simply explains the mathematical underpinnings of DCC and why it can be used to estimate correlations for the variables. Correlatoin is requires covariance. Covariance can be used in Markowitz optimization. Thus optimization by DCC is possible. This method also allows many specifications for the GARCH estimation and is explained.  

The next subsection of Analysis of the potential of GAS-GARCH DCC is a section explaining how one would go about estimating a multivariate GAS and why GAS modeling may be desirable. It explains how it also starts at a DCC specification. However while other DCC specification of GARCH will collapse back into a DCC structure GAS modeling does not. This means that estimates of one covartiance term is dependent on other terms which means the model parameters expand and this makes estimating the models much more data and computationally expensive.

The alternative specifications of GARCH that are used are then explained in the next subsection section and why they are useful in GARCH estimation. 

The Data section aims to explain the limitations of the data. We are trying to beat the JSE based on real returns and correlations. Imputing returns or using returns for the stock when they are not in the index might lead to overperformance or under performance biasing the results, drawing from a distribution can also mean that real returns are not used.  The rebalancing periods are thus explained and set up such that it is in accordance to the JSE rebalancing period. This means fees will be in line with JSE at least for the non short portfolios. Further explaination of this is given

Then a Graph of the various stocks making up the index is shown to show the stocks price over time and how they are removed at various points. Johansen to identify autocrellation not shown.

The results section is next. This explains the various optimization limitations placed on the various portfolios before showing the portfolio results. The portfolio results come from the PortflioGeneration.R script that contains the code that can find the various models monthly returns. It is saved in that script as the calculations have many convergence errors that must be fixed and cannot be knit into the document for that reason.
The results section then show how the various models perform. The portfolio optimizations that allow for large shorts are the most volitile and dominate the returns and thus returns and cumulative returns graphs are split to remove them for better analysis of the remaining stocks. It is also clear that some portfolios have massive negative returns in 2008 as the stocks fall rapidly. This results in negative cumulative returns.  For this reason it is also usefull to split up the models for the entire period and a post 2008 period and examine whether there are large diffences over time. Correlation analysis also tells an interesting story here.
 The final part of the results section continues with this time split but examines downside risks, sharpe ratios and annualized returns to find what are the best performing portfolios and whether they can beat the JSE

Finally the conclustion summarises the findings over the entire paper focussing on the correlation analysis results the problem with GAS multivariate, potential for diversification and an optimal set of conditions as well as robustness of the models. Also admitting that across the entire period portfolio optimization failed but that post crisis it massively outperformed. A method to seperate crisis from non crisis might lead to significant improvement

Regarding the PortfolioGeneration.R script. This code could potentially be significantly improved through the use of functions.  In the code the full explanation of each piece of code is provided in the comments. The problem with this code is that errors of failure to converge happen a few times with every model.  You set the type of model that you need in the univariate specification and then run the entire while loop. This will run through the models and provide monthly return data for optimal weight portfolios in the order of almost no limit shorts, limited shorts and then minimum investment maximim investment limitations.  These results can then be saved in Model1,2 and 3. Change specification and then run again and save into the Model4,5 and 4 etc.. If there were no errors this could be changed easily into functional form which would be much neater and easier. However it will give error convergence. This means that for that period that is printed it could not find a solution to the multivariate GARCH fit. If this is the case running the monthly returns estimation part will provide the returns at the same weights that were in place prior to the error of no convergence. Change whiledate to the next monthly return that is needed minus 3 months. Run the code from where stated across the loop. sGARCH should give 1-2 errors same for iGARCH, EGARCH will give 3 or 4, while GJRGARCH will give a large number. If I knew how to catch functions errors it is possible to write if statements to automatically do the correct piece of code and continue the estimation.  Once all models are saved into models the rest of the code can be run. This will simply combine the various models into a single dataframe and save that dataframe for use in the results sections. Finding the optimizations is computationaly intensive and will take considerable time to achieve.
