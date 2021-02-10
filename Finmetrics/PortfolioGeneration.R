
library(rmsfuns)
library(tidyverse)
library(xtable)
library(IntroCompFinR)
pacman::p_load(tbl2xts)
pacman::p_load("xts", "tidyverse", "tbl2xts", "PerformanceAnalytics",
               "lubridate", "glue","roll")
pacman::p_load("dplyr")



DailyTop40 <- read_rds( "data/DailyTop40.rds")

pacman::p_load("MTS", "robustbase")
pacman::p_load("tidyverse", "devtools", "rugarch", "rmgarch",
               "forecast", "tbl2xts", "lubridate", "PerformanceAnalytics",
               "ggthemes", "parallel" )
#Create a dataset that can be adjusted without potentially change original dataset
WorkData <- DailyTop40

#first calculate monthly returns of JSE Top 40 to provide a control portfolio this is done by multiplying weights by the value of a stock to get index value and then finding the log returns
JSESeries <- WorkData %>% group_by(date) %>%
    summarise(IndexPrice = sum(PX_LAST * J400_W_Adj, na.rm = TRUE))
dlogrtn <- JSESeries  %>%
    mutate(dlogret = log(IndexPrice) - log(lag(IndexPrice))) %>%  filter(date > first(date))
# select only needed columns and remove NA
dlogrtn <- dlogrtn %>% select(date, dlogret)
dlogrtn[is.na(dlogrtn)] <- 0
#create the returns series that will be filled and the start date 2008-03 is when the models start so that is start for JSE Top 40
Returnseries <- as_tibble()
whiledate <- as.Date("2008-03-01")
# While loop to calculate monthly returns over full time period
while (whiledate <"2019-09-01")
{
    # create a start date for the month
    print( whiledate)
    aftertime <- paste(whiledate,"/", sep="")

    # add 1 month to loop variable
    nodate <- whiledate %m+% months(1)
    whiledate <- as.Date(nodate)

    # Here the variable that will be rbinded to returns is created with the date of the month also create a temporary dataset that can be manipulated of returns
    prereturnseries <- whiledate %>% as_tibble()
    dlogperiod<- dlogrtn %>%tbl_xts()
    # create end date for the month
    beforetime<- paste("/",whiledate, sep="")

    # remove dates not in the month
    dlogperiod <- dlogperiod[beforetime]
    dlogperiod <- dlogperiod[aftertime]
    # Can sum over returns of the month thanks to log returns and then bind into series
    monthreturn <- colSums(dlogperiod)
    prereturnseries <- cbind(prereturnseries, monthreturn)
    Returnseries <- rbind(Returnseries,prereturnseries)

}


# Remove nonessensial data and spread to allow for future calculation
Stocksseries <-  WorkData %>%  tbl_xts(cols_to_xts = "PX_LAST", spread_by = "Tickers") %>% xts_tbl()
#find the returns series
Stocksreturns <- Stocksseries %>% gather(Tickers, PX_LAST, -date) %>% group_by(Tickers) %>%
    mutate(dlogret = log(PX_LAST) - log(lag(PX_LAST)))

Spreadreturns <- Stocksreturns %>%  tbl_xts(cols_to_xts = "dlogret", spread_by = "Tickers")
returnloop <- Return.annualized(Spreadreturns, scale = NA, geometric = TRUE)

#remove NA in spreadreturns note stockreturns will maintain NA values
SpreadreturnsNoNA<- Spreadreturns
SpreadreturnsNoNA[is.na(SpreadreturnsNoNA)] <- 0



pacman::p_load(RiskPortfolios)
# dplyr causing problem with finding first for a date for some reason so manually remove the NA of the first day
Spreadreturns <-  Spreadreturns["2005-01-04/"]
SpreadreturnsNoNA <- SpreadreturnsNoNA["2005-01-04/"]



#Creaste the models return series
Returnseries1 <- as_tibble()
Returnseries2 <- as_tibble()
Returnseries3 <- as_tibble()
#creates a loop variable for the start of the period
whiledate <- as.Date("2008-01-01")

#alternative start date for start for continuing after a failure to converge or other model error
whiledate <- as.Date("2019-01-01")

while (whiledate <"2019-07-01")
{
    print( whiledate)
    # move date forward in 3 month intervals
    nodate <- whiledate %m+% months(3)
    whiledate <- as.Date(nodate)
    # create lookback period
    back1year <- whiledate %m+% months(-12)

    timeSigma <- paste(back1year,whiledate, sep="/", collapse=NULL)
    loopStockseries <- Spreadreturns[timeSigma]
    # this will remove all stocks for which no observations have been noted in lookback period
    naremoved <- loopStockseries  %>% xts_tbl() %>%  select_if(~ !any(is.na(.))) %>% tbl_xts()

    # get a list of the stocks names that are observed in period
    loopnames <-colnames(naremoved, do.NULL = TRUE, prefix = "col")



    # Select from returns
    returnsSelect <- SpreadreturnsNoNA %>% xts_tbl() %>% select(date, all_of(loopnames)) %>% tbl_xts

    beforetime <- paste("/",whiledate, sep="")

    returnsSelectbefore <- returnsSelect[beforetime]

    returnloop <- Return.annualized(returnsSelectbefore, scale = NA, geometric = TRUE)

    #COVARMATs <- cov( naremoved, y = NULL,
    # method = "pearson")*252


    # Specify univariate specification
    # change between the differnt types of specifications by changing model variable
    # this entire loop could be made into a function however it is probably impossible due to the errors that are sometimes generated that require certain steps to be done. It could be done by catching errors and running the correct code in event of error but this is above my coding ability at present
    uspec <- ugarchspec(variance.model = list(model = "gjrGARCH",
                                              garchOrder = c(1, 1)), mean.model = list(armaOrder = c(1,
                                                                                                     0), include.mean = TRUE), distribution.model = "sstd")


    multi_univ_garch_spec <- multispec(replicate(ncol(naremoved), uspec))


    spec.dcc = dccspec(multi_univ_garch_spec, dccOrder = c(1, 1),
                       distribution = "mvnorm", lag.criterion = c("AIC", "HQ", "SC",
                                                                  "FPE")[1], model = c("DCC", "aDCC")[1])  #

    cl = makePSOCKcluster(10)

    # Specify the multivariate specification

    multf = multifit(multi_univ_garch_spec, naremoved, cluster = cl)

    # fit the dcc model

    fit.dcc = dccfit(spec.dcc, data = naremoved, solver = "solnp",
                     cluster = cl, fit.control = list(eval.se = FALSE), fit = multf)
    # extract the covariance matrix
    RcovList <- rcov(fit.dcc)
    covmat = matrix(RcovList, nrow(naremoved), ncol(naremoved) * ncol(naremoved),
                    byrow = TRUE)
    # test for MC
    mc1 = MCHdiag(naremoved, covmat)
    # # find the covariance matrix for the future ie at the end of a 3 month period approximatly 60 trading days multiply by time in order to get a yearly covariance matrix
    dcc.time.var.cov <- rcov(fit.dcc)
    COVARMATsdcc <-  dcc.time.var.cov[, , 60]
    COVARMATs <- COVARMATsdcc*252



    # use mean variance optimization to find optimal portfolio under constraints
    Type = "mv"
    Opt_W1 <-
        optimalPortfolio(mu = returnloop, Sigma = COVARMATs,
                         control = list(type = Type, constraint = 'user',
                                        LB = rep(-0.5, ncol(COVARMATs)),
                                        UB = rep(0.5, ncol(COVARMATs))))

    Opt_W2 <-
        optimalPortfolio(mu = returnloop, Sigma = COVARMATs,
                         control = list(type = Type, constraint = 'user',
                                        LB = rep(-0.05, ncol(COVARMATs)),
                                        UB = rep(0.3, ncol(COVARMATs))))

    Opt_W3 <-
        optimalPortfolio(mu = returnloop, Sigma = COVARMATs,
                         control = list(type = Type, constraint = 'user',
                                        LB = rep(0.01, ncol(COVARMATs)),
                                        UB = rep(0.2, ncol(COVARMATs))))



    # attempt to find tangency portfolio. However numerous errors crop up due to risk free rate being higher than some stocks returns or the minimum variance portfolio for certain time periods
    # tan.port <- tangency.portfolio(returnloop, COVARMATs, 0.03)
    # print(tan.port)




    print(Opt_W1)
    # start of the monthly return calculation. Should an error of covergence occur the same weights can be used by running from here until end of loop manually. Ensure that the whiledate variable is correct in this event or there may be duplicate dates or missing dates in returns series nd the model inaccurate
    # create a new loop for 3 monthly loop
    looptime <- whiledate
    #forecast is 3 months or 1 quater of monthly returns
    forecast <- 1:3
    # Just reduce data to 3 month timeframe of interest
    aftertime <- paste(whiledate,"/", sep="")
    returnsSelectFor <- returnsSelect[aftertime]
    Afterdate <- whiledate %m+% months(3)
    forendtime <- paste( "/", Afterdate, sep="")
    returnsSelectFor <- returnsSelectFor[forendtime]

    # Create three returns weights from above optimization. This entire process can be streamlined through using another for loop to reduce code lines
    returnsweights1 <- Return.portfolio(returnsSelectFor,Opt_W1)
    returnsweights2 <- Return.portfolio(returnsSelectFor,Opt_W2)
    returnsweights3 <- Return.portfolio(returnsSelectFor,Opt_W3)


    for(i in forecast)
    {
        #create the variables to be rbind to returns series with the month of the returns
        prereturnseries1 <- looptime %>% as_tibble()
        prereturnseries2 <- looptime %>% as_tibble()
        prereturnseries3 <- looptime %>% as_tibble()
        # set columnames correctly for rbind
        colnames(prereturnseries1) <-( "date")
        colnames(prereturnseries2) <-( "date")
        colnames(prereturnseries3) <-( "date")
        # remove data before month of interest
        aftertime <- paste(looptime,"/", sep="")

        returnsweightsFor1 <- returnsweights1[aftertime]
        returnsweightsFor2 <- returnsweights2[aftertime]
        returnsweightsFor3 <- returnsweights3[aftertime]

        # remove data after month of interest
        looptime <- looptime %m+% months(1)
        forendtime <- paste( "/", looptime, sep="")

        returnsweightsFor1 <- returnsweightsFor1[forendtime]
        returnsweightsFor2 <- returnsweightsFor2[forendtime]
        returnsweightsFor3 <- returnsweightsFor3[forendtime]
        # sum over columns to get monthly returns thanks to logreturns
        monthreturn1 <- colSums(returnsweightsFor1)
        monthreturn2 <- colSums(returnsweightsFor2)
        monthreturn3 <- colSums(returnsweightsFor3)


        # cbind to get full variables to be rbound to returnsseries and rbind
        prereturnseries1 <- cbind(prereturnseries1, monthreturn1)
        Returnseries1 <- rbind(Returnseries1,prereturnseries1)
        prereturnseries2 <- cbind(prereturnseries2, monthreturn2)
        Returnseries2 <- rbind(Returnseries2,prereturnseries2)
        prereturnseries3 <- cbind(prereturnseries3, monthreturn3)
        Returnseries3 <- rbind(Returnseries3,prereturnseries3)

    }






}


#specifications 1-3 is sGARCH 4-6 is iGARCH 7 to 9 is EGARCH and 10-12 is GJRGARCH
# Save returns series for first model
Model1 <-Returnseries1
Model2 <-Returnseries2
Model3 <-Returnseries3
# save returns series of second model
Model4<- Returnseries1
Model5<- Returnseries2
Model6<- Returnseries3


# save returns series of third model
Model7<- Returnseries1
Model8<- Returnseries2
Model9<- Returnseries3

#  GJRGARCH gave considerable problems so sometimes saving the data to later return it to Returnseries just made getting the correct model easier also saves the GJRGARCH
Model10<-Returnseries1
Model11<-Returnseries2
Model12<-Returnseries3

Returnseries1 <-  Model10
Returnseries2 <- Model11
Returnseries3 <- Model12
# Changing the variable name to ensure that there is no overwrite by accident
SGARCHn<- Model1
SGARCHl<-Model2
SGARCHns<-Model3

IGARCHn<- Model4
IGARCHl<-Model5
IGARCHns<-Model6

EGARCHn<- Model7
EGARCHl<-Model8
EGARCHns<-Model9

GJRGARCHn<- Model10
GJRGARCHl<-Model11
GJRGARCHns<-Model12



# this is irrevalant as columnames in most cases are not kept makes neat variables though
SGARCHn <- SGARCHn %>% rename(sGARCHnolim = monthreturn1)
SGARCHl <- SGARCHl %>% rename(sGARCHlim = monthreturn2)
SGARCHns <- SGARCHns %>% rename(sGARCHnoshort = monthreturn3)
IGARCHn <- IGARCHn %>% rename(iGARCHnolim = monthreturn1)
IGARCHl <- IGARCHl %>% rename(iGARCHlim = monthreturn2)
IGARCHns <- IGARCHns %>% rename(iGARCHnoshort = monthreturn3)
EGARCHn <- EGARCHn %>% rename(eGARCHnolim = monthreturn1)
EGARCHl <- EGARCHl %>% rename(eGARCHlim = monthreturn2)
EGARCHns <- EGARCHns %>% rename(eGARCHnoshort = monthreturn3)
GJRGARCHn <- GJRGARCHn %>% rename(gjrGARCHnolim = monthreturn1)
GJRGARCHl <- GJRGARCHl %>% rename(gjrGARCHlim = monthreturn2)
GJRGARCHns <- GJRGARCHns %>% rename(gjrGARCHnoshort = monthreturn3)

Returnseries <- Returnseries %>% rename(JSETop40 = monthreturn)

# create composite of all models
AllModelReturn<- cbind(SGARCHn,sGARCHlim = SGARCHl[,2],sGARCHnoshort = SGARCHns[,2],iGARCHnolim = IGARCHn[,2],iGARCHlim = IGARCHl[,2],iGARCHnoshort = IGARCHns[,2],eGARCHnolim = EGARCHn[,2],eGARCHlim = EGARCHl[,2],eGARCHnoshort = EGARCHns[,2],gjrGARCHnolim = GJRGARCHn[,2], gjrGARCHlim =GJRGARCHl[,2],gjrGARCHnoshort = GJRGARCHns[,2], JSETop40 = Returnseries[,2])
# create composite of all low leverage and no leverage models seperation makes later estimation slightly easier
LowerLevModelReturn<- cbind(SGARCHl ,sGARCHnoshort = SGARCHns[,2],iGARCHlim = IGARCHl[,2],iGARCHnoshort = IGARCHns[,2],eGARCHlim = EGARCHl[,2],eGARCHnoshort = EGARCHns[,2], gjrGARCHlim =GJRGARCHl[,2],gjrGARCHnoshort = GJRGARCHns[,2], JSETop40 = Returnseries[,2])

AllModelReturntbl<- AllModelReturn %>%   gather(Portfolios, returns , -date)

# save the models permently as this modeling process is intricate and cannot be redone automatically in the current coding structure due to errors of convergence in the DCC proces. This saved data can be used for further analysis
saveRDS(AllModelReturn, file = "ModelReturnsSpread.RDS")
saveRDS(LowerLevModelReturn, file = "LowerleverageModelReturns.RDS")

