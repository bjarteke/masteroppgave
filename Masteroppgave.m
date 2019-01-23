%////////////////////////////////////////////////////////////////
%///////////////////////////Simulation data//////////////////////
%////////////////////////////////////////////////////////////////
numSims = 10000;
numIntervals = 4;
storeValues = false;
storeValuesMean = false;
 
initProfitShare = 0.8;
initWeightEquity = 0.2;
initTAmax = [0.05];
initKRFmax = 0.09;
g = 0.035;               %guaranteed rate
 
 
reduceKRFInPayoutYears = true;
 
outputFile = '\\sambaad.stud.ntnu.no\bjarteke\Documents\Fordypningsemne - finans\newTesting.xlsx';
 
inputValues = initTAmax;
testInput = initTAmax;
 
%////////////////////////////////////////////////////////////////
%/////////////////////////////Lifetime///////////////////////////
%////////////////////////////////////////////////////////////////
%The reserves are defined as the future cashflows discounted by the guaranteed rate. 
promisedPayment = 100;   %what yearly cash flow the policyholder is promised after retiering
yearsToRetirement = 30;
pensionAge = 67;
deathAge = 100;          
lastPaymentAge = 90;     %the years from 90-100 are not simulated, but the remaining value is paid out
listSize = lastPaymentAge+1-pensionAge;
realPayments = zeros([1 listSize]);
 
yearsFromRetirementToDead = listSize;
numYears = yearsToRetirement + yearsFromRetirementToDead;
 
%////////////////////////////////////////////////////////////////
%////////////////////////Lists to Save///////////////////////////
%////////////////////////////////////////////////////////////////
insurerProfitStored = zeros([1 numYears]);
equityStored = zeros([1 numYears]);
bondReturnsStored = zeros([1 numYears]);
TAStored = zeros([1 numYears]);
KRFStored = zeros([1 numYears]);
reservesStored = zeros([1 numYears]);
bookReturnStored = zeros([1 numYears]);
unrealizedReturn = zeros([1 numYears]);
totalReturnStored = zeros([1 numYears]);
profitShareReservesStored = zeros([1 numYears]);
weightsEquityStored = zeros([1 numYears]);
equityReturnStored = zeros([1 numYears]);
unRealizedStored = zeros([1 numYears]);
preBookCashStored = zeros([1 numYears]);
divStored = zeros([1 numYears]);
portfolioStored = zeros([1 numYears]);
KRFMaxStored = zeros([1 numYears]);
customerPayments = zeros([1 numYears]);
%bondWeightsStored = zeros([1 numYears]);
newEquityStored = zeros([1 numYears]);
equitySoldStored = zeros([1 numYears]);
equityAddedStored = zeros([1 numYears]);
realizedCashStored = zeros([1 numYears]);
couponStored = zeros([1 numYears]);
inputKRFstored = zeros([1 numYears]);
prePortfolioStored = zeros([1 numYears]);
 
promisedProfitShareStored = zeros([1 numYears]);
PVinsurersProfitStored = zeros([1 numYears]);

rateStored = zeros([numYears*numIntervals+1 numSims]);
 
%////////////////////////////////////////////////////////////////
%//////////////////////////Interest Rate/////////////////////////
%////////////////////////////////////////////////////////////////
zeroRate = [0.01172, 0.01452, 0.01640, 0.01780, 0.01901, 0.01981, ...
    0.02054, 0.02117, 0.02173, 0.02221, 0.02270, 0.02323, 0.02376, ...
    0.02428, 0.02480, 0.02531, 0.02579, 0.02626, 0.02672, 0.02715, ...
    0.02756, 0.02796, 0.02834, 0.02870, 0.02905, 0.02938, 0.02970, ...
    0.03000, 0.03029, 0.03056, 0.03082, 0.03108, 0.03132, 0.03155, ...
    0.03177, 0.03198, 0.03219, 0.03238, 0.03257, 0.03275, 0.03292, ...
    0.03308, 0.03324, 0.03340, 0.03354, 0.03369, 0.03382, 0.03395, ...
    0.03408, 0.03420, 0.03432, 0.03444, 0.03455, 0.03465, 0.03475, ...
    0.03475, 0.03475, 0.03475, 0.03475, 0.03475, 0.03475, 0.03475, ...
    0.03475, 0.03475, 0.03475, 0.03475, 0.03475, 0.03475, 0.03475]';
 
%Converting the zeroRates into compounded rates
for x=1:size(zeroRate)
   compoundRate = log(1+zeroRate(x));
   zeroRate(x) = compoundRate;
end
 
maturities = [1/12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, ...
    17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,... 
    33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, ...
    49, 50, 51, 52, 53, 54, 55]';
sigma = 0.0085;
alpha = 0.0304;
 

corr = 0.249; %Correlation Brownian motion between H&W and BS
 

 
 
%////////////////////////////////////////////////////////////////
%///////////////////////Martingale Testing///////////////////////
%////////////////////////////////////////////////////////////////
martingale = zeros([1 numYears]);
martingaleBonds = zeros([1 numYears]);
martingaleBonds2 = zeros([1 numYears]);
martingaleEquity = zeros([1 numYears]);
 
 
%////////////////////////////////////////////////////////////////
%/////////////////////////Output variables///////////////////////
%////////////////////////////////////////////////////////////////
insurersProfitProfitShare = zeros([1 size(testInput,2)]);
customerProfitProfitShare = zeros([1 size(testInput,2)]);
 
 
 
tic
for testing = 1:size(testInput,2)
    disp("Simulation run: " + testing);
    randn('seed', 0);
    
    testingTAKRF_ProfitShareInsurer = zeros([numSims numYears]);
    testingTAKRF_EquityToCoverG = zeros([numSims numYears]);
    testingTAKRF_EquityToCoverNegativeReturn = zeros([numSims numYears]);
    testingTAKRF_RateOnInsProfit = zeros([numSims numYears]);
    testingTAKRF_TotalProfitShares = zeros([numSims numYears]);
    testingTAKRF_KRFMax = zeros([numSims numYears]);
    testingTAKRF_TAMax = zeros([numSims numYears]);
    testingTAKRF_TA = zeros([numSims numYears]);
    testingTAKRF_KRF = zeros([numSims numYears+1]);
    testingTAKRF_totalReturn = zeros([numSims numYears]);
    
    
    
for sim = 1:numSims
    testingTAKRF_ProfitShares = 0;
    testingTAKRF_ToCoverG = 0;
    testingTAKRF_NegativeReturn = 0;
    testingTAKRF_Rate = 0;
 
    promisedPayment = 100;
    
    [reserves, realPayments] = calcReservesMen(yearsToRetirement,promisedPayment,... 
        g, pensionAge, deathAge, lastPaymentAge, realPayments);
    reservesInit = reserves;
 
    mortalityRates = [];
    numYears = yearsToRetirement + yearsFromRetirementToDead;
    
    
    profitShare = initProfitShare;
    
    %Set initial portfolio
    KRF = 0.03*reserves;
    portfolio = reserves + KRF;
    KRFmaxPercentage = initKRFmax;    %maximum percentage of the reserves and proftitShare reserves
    KRFminPercentage = 0;       %minimum percentage of the reserves and profitShare reserves
    KRFmaxConst = 0;
    
    %Set initial asset mix
    weightEquity = initWeightEquity;
    initialEquityPercentage = weightEquity;
    equityMinLimit = weightEquity - 0.1;
    equityMaxLimit = weightEquity + 0.1;
    totalEquityBought = 0;
    
    weightBonds = 1 - weightEquity;
    
    
    %Bonds
    bondHTMYears = [1,3,5,10];
    initHTMYears = [1,3,5,10];
    bondsHTMValue = weightBonds*portfolio;
    bondHTMWeights = [0.25 0.15 0.10 0.5];
    simBondHTMReturn = [0.0117 0.0163 0.0188 0.0220];
    %simBondHTMReturn = [0.035 0.035 0.035 0.035];
   
    profitShareReserves = 0;  %value added to reserves due to profit sharing
    g_addReserves = 0;       %guaranteed rate on additional reserves   
    
  
 
    TA = 0.05*reservesInit;
    TAlowerLimit = 0;
    TAupperLimit = initTAmax(testing);
 
    insurersProfit = 0;
    customerPayment = 0;
    
    %Black & Scholes
    s0 = portfolio * weightEquity;
    yearlyDiv = 0.02;
    div = yearlyDiv / numIntervals; %Dividend will be handled manually
    vol = .25;
    stockPrices = zeros([numIntervals numYears]);
    stockTesting = zeros([numIntervals numYears]);
    sT = s0;
    sT_testing = s0;
    
    %Create stochastic z values to use in HullWhite and BlackScholes
    z = calcZ(numYears*numIntervals,numIntervals);
    
    %Simulate interest rate 
    hwInit = simulateInterestRate(numYears, numIntervals, z,alpha,sigma,zeroRate);
    rf = hwInit(:,1,:); %Using 1/12 years bonds as a proxy for risk free rate
    rfDisk = hwInit(:,1,:);
    
    %rateStored(:,sim) = rf;

    %For each year
    for year = 1: numYears
        %disp("YEAR: " + year);
        if (year > yearsToRetirement && year <= (numYears))
            promisedPayment = realPayments(year-yearsToRetirement);
            
            reserves = reserves - promisedPayment;
            
            promisedProfit = profitShareReserves/(numYears - year + 1);
            profitShareReserves = profitShareReserves - promisedProfit;
            
            if (storeValues)
                promisedProfitShareStored(year) = promisedProfit;
            end
            
            [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);
            
            %Reducing the stocks and bonds according to their weights. 
            portfolio = portfolio - (promisedProfit + promisedPayment);
            sT = sT - (promisedProfit + promisedPayment)*weightEquity;
            bondsHTMValue = bondsHTMValue - (promisedProfit + promisedPayment)*weightBonds;
            
            if (portfolio - (bondsHTMValue + sT) > abs(0.001))
                disp("Error: portfolio does not equal bond + equity")
            end
     
            if (portfolio < 0)
                disp("ERROR: Portfolio is negative" + portfolio);
            end
            
            
            customerPayment = customerPayment + promisedPayment + promisedProfit;
        end
        
        equityReturnTemp = 0;
        divPayments = 0;
        equityBookReturn = 0;
        
        %store before buffer strategy
        prePortfolioStored(year) = sT + bondsHTMValue;
        
        %For each month/quartal
        for interval = 1:numIntervals
            
           %Simulate interest rate and stock prices
            rf_current = rfDisk(interval + (year-1)*numIntervals);
            z_temp = randn*sqrt(1/numIntervals);
            
            %z value correlated
            z_current = sqrt(1-corr^2)*z_temp + corr*z(interval + (year-1)*numIntervals)*sqrt(1/numIntervals);            
            %Simulate stock prices, and handle dividend payments. 
            stockPrices(interval,year) = sT*exp((rf_current-0.5*vol^2)*1/numIntervals + vol*z_current);
            stockTesting(interval,year) = sT_testing*exp((rf_current-0.5*vol^2)*1/numIntervals + vol*z_current);
            
            divPayments = divPayments + stockPrices(interval,year)*div;
            
            sT = stockPrices(interval,year);
            sT_testing = stockTesting(interval,year);
                        
            %Calculate unrealized return from equity.
            equityReturn = calcEquityReturn(stockPrices, year, interval, numIntervals, s0);
            equityReturnTemp = equityReturnTemp + equityReturn;
            
            insurersProfit = insurersProfit*exp(rf_current/numIntervals);
            
        end

        %Reduce stock price with dividends
        sT = sT - divPayments;
        
        %Portfolio update
        portfolio = bondsHTMValue + sT;
        

        
        %Increasing the insurer's profit by the risk-free rate.
        testingTAKRF_Rate = testingTAKRF_Rate + insurersProfit*rf_current;
        
        
        customerPayment = customerPayment*exp(rf_current);
 
        
        
        
        %___________________________________________________________
        %BOND RETURNS
        %___________________________________________________________
        
        [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);

        [newEquityValue, equitySold, simBondHTMReturn, bondHTMWeights, ... 
            bondsHTMValue, bondHTMYears,  newEquityStored, ... 
            equitySoldStored, equityAddedStored, equityBought, couponPayment] = calcBondHTMReturn2(storeValues, bondHTMYears, hwInit,...
            year, interval, simBondHTMReturn, bondHTMWeights, bondsHTMValue, ...
            initHTMYears, portfolio, weightEquity, ...
            initialEquityPercentage, ... 
            sT, equityMinLimit, equityMaxLimit, newEquityStored, ... 
            equitySoldStored, equityAddedStored, KRF);
        
        %all of the bond return is realized each year since it's coupon
        %bonds
        totalEquityBought = totalEquityBought + equityBought;
        cashFromEquitySold = (equitySold / portfolio) * KRF;
  
        %___________________________________________________________
        %EQUITY RETURNS
        %___________________________________________________________
        equityReturn1Y = equityReturnTemp;
        
        %___________________________________________________________
        %TOTAL RETURNS
        %___________________________________________________________   
        totalReturn = equityReturn1Y * weightEquity;
        
        bookReturn = simBondHTMReturn * transpose(bondHTMWeights);
  
        
        totalUnrealizedReturn = totalReturn - bookReturn;
        
        
        if (storeValuesMean)
           testingTAKRF_totalReturn(sim,year) = totalReturn; 
        end
        
        %Calculate new bond and equity weights
        sT = newEquityValue;
        [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);
       
        %Store values
        if (storeValues)
            weightsEquityStored(year) = weightEquity;
            equityReturnStored(year) = equityReturn1Y;
            totalReturnStored(year) = totalReturn;
            bondReturnsStored(year) = bookReturn; 
        end
        
        %used to calculate the decrease of the TA in the payout years
        TAoldMaxValue = TAupperLimit*(reserves + profitShareReserves);
 
        %Allocate book return (g to reserves, and the rest to TA/profit share).
        realizedFromKRF = 0;
        
        %Set the new values of the reserves and reduce equity with the
        %equity sold, calculated in newEquityValue
        reserves = reserves*exp(g);
        
        stockPrices(interval,year) = sT;
        
        portfolio = sT + bondsHTMValue;
        KRF = portfolio - (reserves*exp(g) + profitShareReserves + TA);
        
        if (storeValues)
            preBookCashStored(year) = bookReturn;
        end
  
        %Run buffer strategy
        [TA,insurersProfit, reserves, profitShareReserves, KRFStored, TAStored, reservesStored, bookReturnStored, ... 
            insurerProfitStored, profitShareReservesStored, KRFMaxStored, KRF, TAtoPayment,testingTAKRF_ToCoverG,testingTAKRF_NegativeReturn, testingTAKRF_ProfitShares,...
    testingTAKRF_TA, testingTAKRF_TAMax,testingTAKRF_KRF,testingTAKRF_KRFMax, KRFmaxConst, portfolio, realizedCashStored, couponStored, ...
    inputKRFstored, sT, bondsHTMValue] = ...
            runBufferStrategy(sim, storeValues, storeValuesMean, reduceKRFInPayoutYears, bookReturn, g, TA, TAupperLimit, TAlowerLimit, insurersProfit, ...
            reserves, profitShareReserves, profitShare, KRF, KRFmaxPercentage, KRFminPercentage, totalUnrealizedReturn, ... 
            KRFStored, TAStored, reservesStored, bookReturnStored, insurerProfitStored, profitShareReservesStored, year, ...
            yearsToRetirement, yearsFromRetirementToDead, divPayments, cashFromEquitySold, realizedFromKRF, portfolio, totalReturn, KRFMaxStored, TAoldMaxValue, ...
            testingTAKRF_ToCoverG,testingTAKRF_NegativeReturn, testingTAKRF_ProfitShares,...
    testingTAKRF_TA, testingTAKRF_TAMax,testingTAKRF_KRF,testingTAKRF_KRFMax, KRFmaxConst, couponPayment, cashFromEquitySold, ...
    realizedCashStored, couponStored, inputKRFstored, sT, bondsHTMValue);
        
        liabilities = reserves + TA + profitShareReserves;
        
        portfolio = sT + bondsHTMValue;
        
        if (year > yearsToRetirement && year < (numYears))
            customerPayment = customerPayment + TAtoPayment; 
            weightEquity = sT / portfolio;
            weightBonds = 1 - weightEquity;
            
            %Reducing the stocks and bonds according to their weights. 
            sT = sT - TAtoPayment*weightEquity;
            bondsHTMValue = bondsHTMValue - TAtoPayment*weightBonds;
        end
        
        %update portfolio;
        [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);
        
        %Adding from the last run of buffer strategy (for instance when
        %money is added to Reserves this year). 
        if (year == numYears)
            customerPayment = customerPayment + profitShareReserves + reserves + TA;
            [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);
            portfolio = portfolio - profitShareReserves;
            sT = sT - profitShareReserves*weightEquity;
            bondsHTMValue = bondsHTMValue - profitShareReserves*weightBonds;
            reserves = 0;
            TA = 0;
            profitShareReserves = 0;
        end
        
        customerPayments(year) = customerPayment;
        
        %store portfolio
        if (storeValues)
            unRealizedStored(year) = totalUnrealizedReturn; 
            portfolioStored(year) = portfolio;
        end
        
        if storeValuesMean
           discountFactor = exp(-calcBankAccountRate(rfDisk(1:(year)*numIntervals),numIntervals));
           testingTAKRF_EquityToCoverG(sim,year) = testingTAKRF_ToCoverG * discountFactor;
           testingTAKRF_ProfitShareInsurer(sim,year) = testingTAKRF_ProfitShares * discountFactor;
           testingTAKRF_EquityToCoverNegativeReturn(sim,year) = testingTAKRF_NegativeReturn * discountFactor;
           testingTAKRF_RateOnInsProfit(sim,year) = testingTAKRF_Rate * discountFactor;
        end
                
        %Martingale testing
        %martingale(year) = martingale(year) + (insurersProfit + portfolio)*exp(-calcBankAccountRate(rfDisk(1:(year)*numIntervals),numIntervals));
        martingaleEquity(year) =  martingaleEquity(year) + sT_testing*exp(-calcBankAccountRate(rfDisk(1:(year)*numIntervals),numIntervals));
        martingaleBonds(year) = martingaleBonds(year) + exp(-calcBankAccountRate(rfDisk(1:(year*numIntervals)),numIntervals));
        %martingaleBonds2(year) = exp(-(hwInit(1,year)*year));    
 
        
        
        if storeValuesMean
            testingTAKRF_TotalProfitShares(sim,year) = insurersProfit * exp(-calcBankAccountRate(rfDisk(1:(year)*numIntervals),numIntervals));
        end
           
    end
    
    %testingTAKRF_KRF(sim,year+1) = KRF*discountFactor;
 
 
    %Print values from one simulation to excel
    if sim == 1 && storeValues
        filename = outputFile;
        sheet = 1;
        xlswrite(filename,stockTesting,sheet,'C2');
        xlswrite(filename,equityReturnStored,sheet,'C7');
        xlswrite(filename,bondReturnsStored,sheet,'C8');
        xlswrite(filename,weightsEquityStored,sheet,'C10');
        xlswrite(filename,totalReturnStored, sheet,'C13');
        %xlswrite(filename,unRealizedStored,sheet,'C14');
        %xlswrite(filename,preBookCashStored,sheet,'C15');
        
        xlswrite(filename,realizedCashStored,sheet,'C17');
        xlswrite(filename,inputKRFstored,sheet,'C18')
        xlswrite(filename,couponStored,sheet,'C19');
        xlswrite(filename,prePortfolioStored,sheet,'C21')
        xlswrite(filename,portfolioStored,sheet,'C22')
        
        xlswrite(filename,reservesStored,sheet,'C25'); 
        
        xlswrite(filename,TAStored,sheet,'C26');
        xlswrite(filename,insurerProfitStored,sheet,'C27');
        xlswrite(filename,profitShareReservesStored,sheet,'C28');
        xlswrite(filename,KRFStored,sheet,'C29'); 
        %xlswrite(filename,promisedProfitShareStored,sheet,'C25');
        xlswrite(filename,customerPayments,sheet,'C32');

        
        %xlswrite(filename,portfolioStored,sheet,'C31');
        xlswrite(filename,KRFMaxStored,sheet,'C36');
        xlswrite(filename,newEquityStored,sheet,'C41');
        xlswrite(filename,equitySoldStored,sheet,'C43');
        xlswrite(filename,equityAddedStored,sheet,'C45');
        disp("Simulering skrevet til Excel");
    end
    
    KRF = updateKRF(portfolio,TA,reserves,profitShareReserves); 
    %Discount with the bank account rate
    PVinsurersProfit = (insurersProfit + KRF*(1-profitShare)) * exp(-calcBankAccountRate(rfDisk(1:(year)*numIntervals),numIntervals));
    PVcustomerProfit = (customerPayment + KRF*(profitShare))* exp(-calcBankAccountRate(rfDisk(1:(year)*numIntervals),numIntervals));
    PVinsurersProfitStored(sim) = PVinsurersProfit;
    
    sT = s0;
    fprintf('Simulation %i of %i finished \n',sim,numSims);
    
    insurersProfitProfitShare(testing) = insurersProfitProfitShare(testing) + PVinsurersProfit;
    customerProfitProfitShare(testing) = customerProfitProfitShare(testing) + PVcustomerProfit;
    
    PVinsurersProfit = 0;
    PVcustomerProfit = 0;
    
    
    
    
end
%filename = '\\sambaad.stud.ntnu.no\bjarteke\Documents\Fordypningsemne - finans\TAKRFSims2.xlsx';
%xlswrite(filename,TAmaxes,1,'A1');
%xlswrite(filename,insurersProfitProfitShare,1,'B1');
%xlswrite(filename,customerProfitProfitShare,1,'C1');
%disp(martingaleBonds / numSims);
%disp(martingaleEquity / numSims);
 
insurersProfitProfitShare(testing) = insurersProfitProfitShare(testing)/numSims;
customerProfitProfitShare(testing) = customerProfitProfitShare(testing)/numSims;
 
disp(insurersProfitProfitShare);
disp(customerProfitProfitShare);
disp(testInput);
disp("KRF: " + initKRFmax)
 
martingaleBonds = martingaleBonds / numSims;
martingaleEquity = martingaleEquity / numSims;

%martingale = martingale/numSims;
 
 
if (storeValuesMean)
    mean_EquityToCoverG = mean(testingTAKRF_EquityToCoverG, 1);
    mean_Rate = mean(testingTAKRF_RateOnInsProfit, 1);
    mean_NegativeReturn = mean(testingTAKRF_EquityToCoverNegativeReturn, 1);
    mean_TotalProfitShares = mean(testingTAKRF_TotalProfitShares, 1);
    mean_ProfitShare = mean(testingTAKRF_ProfitShareInsurer, 1);
 
    mean_TAMax = mean(testingTAKRF_TAMax,1);
    mean_KRFMax = mean(testingTAKRF_KRFMax,1);
    mean_TA = mean(testingTAKRF_TA,1);
    mean_KRF = mean(testingTAKRF_KRF,1);
    mean_Return = mean(testingTAKRF_totalReturn,1);
    
    filename = outputFile;
    sheet = 2;
    xlswrite(filename,mean_EquityToCoverG,sheet,'B2');
    xlswrite(filename,mean_NegativeReturn,sheet,'B3');
    xlswrite(filename,mean_ProfitShare,sheet,'B4');
    xlswrite(filename,mean_Rate,sheet,'B5');
    xlswrite(filename,mean_TotalProfitShares,sheet,'B6');
    xlswrite(filename,mean_TA,sheet,'B7');
    xlswrite(filename,mean_TAMax,sheet,'B8');
    xlswrite(filename,mean_KRF,sheet,'B9');
    xlswrite(filename,mean_KRFMax,sheet,'B10');
    disp("Mean values written til Excel file");
end
 
end
toc
 
 
 
 
%Functions used for updating policy specific variables
function equityReturn = calcEquityReturn(stockPrices, year, interval, numIntervals, s0)   
    if (s0 == 0)
        equityReturn = 0;
    else
        if (year==1 && interval==1)
            equityReturn = log(stockPrices(interval, year)/s0);
 
        elseif (interval==1)
            equityReturn = log(stockPrices(interval,year)/stockPrices(numIntervals,year-1));
 
        else
            equityReturn = log(stockPrices(interval, year)/stockPrices(interval-1, year));
        end
    end          
end
 
function [reserves, realPayments] = calcReservesMen(yearsToRetirement,promisedPayment,... 
        guaranteedRate, pensionAge, deathAge, lastPaymentAge, realPayments)
 
    %base year is the year the formula from Finanstilsynet 
    baseYear = 2013;
    startYear = 2018 + yearsToRetirement;
    discountBaseYear = 2018;
    lifetime = deathAge - pensionAge;
    
    currentAge = pensionAge;
    currentYear = startYear;
    
    temp_reserves = 0;
    temp_lastPayment = 0;
    
    for i = 1:(lifetime+1)
        mu_2013 = (0.189948+0.003564*10^(0.051*currentAge))/1000;
        w = min(2.6714548 - 0.17248*currentAge + 0.001485*currentAge^2,0);
        mu_T = mu_2013*(1+w/100)^(currentYear-baseYear);
        
        if (i==1)
            cumProb = 1-mu_T;
        else
            cumProb = (1-mu_T)*cumProb;
        end
        
        payment = promisedPayment*cumProb;
        
        discountFactor = exp(guaranteedRate*(currentYear-discountBaseYear));
        discountedPayment = payment/discountFactor;
        
        %the realPayments list is used in the years of payout later
        
        %if we are in the years of payout
        if (i < (lastPaymentAge+1 - pensionAge))
            realPayments(i) = payment;
            temp_reserves = temp_reserves + discountedPayment;
        
        %if we are the last living year
        elseif (i == (lifetime+1))
            realPayments(end) = temp_lastPayment+payment;
        
        %if we are in the years from 90-100
        else
            temp_lastPayment = temp_lastPayment + payment;
        end
        
        currentAge = currentAge + 1;
        currentYear = currentYear + 1;
    end
    
    %the years from 90-100 is gathered in one sum, which is discounted
    %using the year the pensjonist is 90 years old
    lastPayment = realPayments(end);
    lastDiscountYear = startYear + lastPaymentAge - pensionAge; 
    discountFactor = exp(guaranteedRate * (lastDiscountYear-discountBaseYear));
    
    lastReserve = lastPayment/discountFactor;
    
    reserves = temp_reserves + lastReserve;
end
 
function [newEquityValue, equitySold, currentBondReturns, updatedWeights, totalBondsHTMValue, ...
    bondHTMyears, newEquityStored, equitySoldStored, equityAddedStored, equityBought, couponPayment] = ... 
    calcBondHTMReturn2(storeValues, bondHTMyears, hw, year, interval, currentBondReturns, bondHTMWeights, ... 
    totalBondsHTMValue, initHTMYears, portfolio, weightEquity, ...
    initialEquityPercentage, currentEquityValue, equityMinLimit, equityMaxLimit, ...
     newEquityStored, equitySoldStored, equityAddedStored, KRF)
    
    %Bond value should be constant for bonds traded at par. Coupon payment
    %is returned
    couponPayment = (exp((bondHTMWeights * transpose(currentBondReturns)))-1) * totalBondsHTMValue ;
    
    if storeValues
        %bondWeightsStored(year) = bondHTMWeights;
    end
    
    %rebalancing
    weightEquity = currentEquityValue/portfolio; 
    equityValueNeeded = portfolio * initialEquityPercentage - currentEquityValue;
    equityValueAdded = 0;
    finished = false;
    equitySold = 0;
    equityBought = 0;
    newEquityValue = currentEquityValue;
    equityAlreadySold = false;
    
    %As the size of the list increases when we add a bond, we have to set
    %the number of iterations initially.
    numIterations = size(bondHTMyears, 2);
    
    %needed for rebalancing later
    bondsHTMValue = bondHTMWeights*totalBondsHTMValue;
    
    for bond=1:numIterations
       %Calculate yearly return from HTM bonds
       rate = currentBondReturns(bond);
       weight = bondHTMWeights(bond);
 
       %decrease time to maturity by 1. 
       bondHTMyears(bond) =  bondHTMyears(bond) - 1;
       
       %buy new bond if bondHTMyear = 0
       if (bondHTMyears(bond)==0)
 
           %If equity value is too small, sell what we have in the bond that expires and buy equity
           if (weightEquity < equityMinLimit && finished == false)
                
                if (bondsHTMValue(bond) >= equityValueNeeded)
                    currentAdding = equityValueNeeded;
                    equityValueNeeded = 0;
                    bondsHTMValue(bond) = bondsHTMValue(bond) - currentAdding;
                    
                else
                    currentAdding = bondsHTMValue(bond);
                    equityValueNeeded = equityValueNeeded - currentAdding;
                    bondsHTMValue(bond) = bondsHTMValue(bond) - currentAdding;
                    
                end
                
                equityValueAdded = equityValueAdded + currentAdding;
                totalBondsHTMValue = totalBondsHTMValue - currentAdding;
                equityBought = equityBought + currentAdding;
                
                if (abs(initialEquityPercentage-weightEquity) < 0.005)
                   finished = true;
                end
                
                newEquityValue = currentEquityValue + equityValueAdded;
                weightEquity = newEquityValue / portfolio;
                
           %Sell equity and buy bonds if equity value exceeds the upper value.    
           elseif (weightEquity > equityMaxLimit && ~equityAlreadySold)
                available = (weightEquity - initialEquityPercentage)*portfolio;
                equitySold = available;
                newEquityValue = initialEquityPercentage*portfolio;
                
                %We can invest in bonds the equity sold minus the part that
                %is realized from the KRF and goes into the book return
                available = available - (equitySold/portfolio)*KRF;
                
                
                totalBondsHTMValue = totalBondsHTMValue + available;
               
                bondHTMWeights = [bondHTMWeights, available/totalBondsHTMValue];
                
                bondsHTMValue = [bondsHTMValue, available];
                bondHTMyears = [bondHTMyears, 1];
                
                rateIndex = year*interval + 1;
                %rateIndex = (year+1)*interval + 1;
                
                currentBondReturns = [currentBondReturns, hw(rateIndex,2)];
                equityAlreadySold = true;
                
                
           end
 
           %Selecting which bond to buy
           if (bond <=4)
               buyYear = initHTMYears(bond);
           else
               rnd = randi(4);
               buyYear = initHTMYears(rnd);
           end
           
           bondHTMyears(bond) = buyYear;
           
           %Set new bond return
           rateIndex = year*interval + 1;
           %rateIndex = (year+buyYear)*interval + 1;
           currentBondReturns(bond) = hw(rateIndex,1);
           
       end
 
    end
    
    totalBondsHTMValue = sum(bondsHTMValue);
    
    updatedWeights = bondsHTMValue / totalBondsHTMValue; 
    
    if storeValues
        equityAddedStored(year) = equityValueAdded;
        newEquityStored(year) = newEquityValue;
        equitySoldStored(year) = equitySold;
    end

end
 
 
%Functions used for simulating HW and BS
function zArray = calcZ(number,numIntervals)
    temp = zeros([number+50 1]);
    for x = 1:number
        temp(x) = randn;
    end
    zArray = temp;
end

function rateArray = simulateInterestRate(years, interval, Z, alpha,sigma,zeroRate)
    rateArray = zeros([years*interval 5]);

    nPeriods = interval*years + 50;
    deltaTime = 1/interval;
    count = 1;
    countIntv = 1;
    yPrev = 0;
    
    for t=1:nPeriods
        B0t = 1 - exp(-alpha * (t - 0));
        B02dt = 1 - exp(-alpha * (2*deltaTime -0));

        if  mod(t-1,interval) == 0
            countIntv = 1;
            forwardRate = zeroRate(count);
            count = count + 1;
        else
            rate1 = zeroRate(count-1);
            rate2 = zeroRate(count);
            delta = (rate2 - rate1) / interval;
            forwardRate = rate1 + countIntv * delta;

            countIntv = countIntv + 1;
        end

        a = forwardRate + (sigma^2/2) * B0t^2;
        y = exp(-alpha*(deltaTime))*yPrev + sqrt(0.5*sigma^2*B02dt)*Z(t);
        r = a + y;
        rateArray(t,1) = r;
        rateArray(t,2) = r;
        rateArray(t,3) = r;
        rateArray(t,4) = r;
        rateArray(t,5) = r;

        yPrev = y;
    end
end
 

 
function bankAccountRate = calcBankAccountRate(rf_array, numIntervals)
    temp = 0;
    for i = 1:size(rf_array)
        temp = temp + rf_array(i)*(1/numIntervals);
    end
    
    bankAccountRate = temp;
end 
 
 
 
%Buffer strategy
function [TA,insurersProfit, reserves, profitShareReserves, KRFStored, TAStored, reservesStored, bookReturnStored, ... 
    insurerProfitStored, profitShareReservesStored, KRFMaxStored, KRF, TAtoPayment,testingTAKRF_ToCoverG,testingTAKRF_NegativeReturn, testingTAKRF_ProfitShares...
    testingTAKRF_TA, testingTAKRF_TAMax,testingTAKRF_KRF,testingTAKRF_KRFMax, KRFmaxConst, portfolio, realizedCashStored, couponStored, ...
    inputKRFstored, sT,bondsHTMValue] = ...
    runBufferStrategy(sim, storeValues, storeValuesMean, reduceKRFInPayoutYears, bookReturn, g, TA, TAupperLimit, TAlowerLimit, insurersProfit, reserves, ...
    profitShareReserves, profitShare, KRF, KRFmaxPercentage, KRFminPercentage, unrealizedReturn, KRFStored, TAStored,...
    reservesStored, bookReturnStored, insurerProfitStored, profitShareReservesStored, year, yearsToRetirement, yearsFromRetirementToDead, divPayment, ...
    bookReturnFromEquitySold, realizedFromKRF, portfolio, totalReturn, KRFMaxStored, TAoldMaxValue,testingTAKRF_ToCoverG,testingTAKRF_NegativeReturn, testingTAKRF_ProfitShares,...
    testingTAKRF_TA, testingTAKRF_TAMax,testingTAKRF_KRF,testingTAKRF_KRFMax, KRFmaxConst, couponPayment, equitySold, ...
    realizedCashStored, couponStored, inputKRFstored,sT,bondsHTMValue)

    TAtoPayment = 0;
    TAupperValue = TAupperLimit*(reserves+profitShareReserves);
    
    sT_init = sT;
    bondsHTMValue_init = bondsHTMValue;
    
    [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);

    %old reserves to calculate KRF
    
    reserves = reserves/exp(g);
    
    %Calculating the guaranteed payment
    guaranteedPayment = reserves*(exp(g)-1);

    %The return of the portfolio could be positive or negative
    %KRF > 0 means return is positive
    %KRF < 0 means return is negative
    KRF = portfolio - (reserves + TA + profitShareReserves);
    
    %disp("KRF: " + KRF);
    %disp(" ");
    
    realizedCash = divPayment + couponPayment + equitySold;
    sT = sT + realizedCash*weightEquity;
    bondsHTMValue = bondsHTMValue + realizedCash*weightBonds;
    [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);

    if storeValues
        inputKRFstored(year) = KRF;
        couponStored(year) = couponPayment;
        realizedCashStored(year) = realizedCash;
        
    end
    
    %available for distribution. Could be negative or positive
    available = KRF + realizedCash;
    reserves = reserves*exp(g);
    
    %Setting the maximum value of the KRF
    KRFmax = KRFmaxPercentage * reserves;
    
    %Years of payout
    if (year < yearsToRetirement)
        KRFmax = KRFmaxPercentage * reserves;
    elseif (year == yearsToRetirement)
         KRFmaxConst = KRFmaxPercentage * reserves;
         KRFmax = KRFmaxConst;
    else
         KRFmax = KRFmaxConst;
    end
    
    if (reduceKRFInPayoutYears)
        if (year >= yearsToRetirement)
            if (KRF < KRFmax)
               newLimitWish = max(KRFmax - (KRFmax/yearsFromRetirementToDead) * (year-yearsToRetirement) ,1);
               if (KRF < newLimitWish)
                   KRFmax = newLimitWish;
               end
            end
        end
    end
    
    if storeValues
        %KRFMaxStored(year) = KRFmax;
    end
    
    
    %____________________________________________________________________
    % If available is negative, then the insurer's has to inject money
    % into the portfolio to cover negative return
    if (year < yearsToRetirement)
        TAoldMaxValue = TAupperValue;
    end
    
    
    if (available < 0)
        insurersProfit = insurersProfit + available;
        [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);
        sT = sT - available*weightEquity;
        bondsHTMValue = bondsHTMValue - available*weightBonds;
        [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);

        %from 0 and up to g
        neededCash = guaranteedPayment;
        
        %We have enough in the TA to cover the guaranteed rate.
        if (neededCash < TA)
            TA = TA - neededCash;
            neededCash = 0;
        
        %We empty the TA, and take from the insurer's profit to cover the
        %guaranteed rate.
        else
            neededCash = neededCash - TA;
            TA = 0;
            insurersProfit = insurersProfit - neededCash;
            [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);
            sT = sT + neededCash*weightEquity;
            bondsHTMValue = bondsHTMValue + neededCash*weightBonds;
            [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);
            testingTAKRF_ToCoverG = testingTAKRF_ToCoverG - neededCash;
        end
    
    %The available amount is less than g but bigger than 0    
    elseif (available >= 0 && available < guaranteedPayment)
        neededCash = guaranteedPayment - available;
        
        
        %We have enough in the TA to cover the guaranteed rate.
        if (neededCash < TA)
            TA = TA - neededCash;
            neededCash = 0;
        
        %We empty the TA, and take from the insurer's profit to cover the
        %guaranteed rate.
        else
            neededCash = neededCash - TA;
            TA = 0;
            insurersProfit = insurersProfit - neededCash;
            [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);
            sT = sT + neededCash*weightEquity;
            bondsHTMValue = bondsHTMValue + neededCash*weightBonds;
            [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);
            testingTAKRF_ToCoverG = testingTAKRF_ToCoverG - neededCash;
        end
        
        
        %We have available amount over g. Things should be distributed
    else
        
        if (realizedCash >= guaranteedPayment)
            realizedCash = realizedCash - guaranteedPayment;
            
            %Put additional profit in TA if new TA is below TAupperValue
            if (TA + realizedCash < TAoldMaxValue)    
                TA = TA + realizedCash;

            %Else: Fill TA up to TAupperValue, and profit share the rest
            else
                realizedCash = realizedCash - (TAoldMaxValue - TA);
                TA = TAoldMaxValue;
                
                %Profit sharing the rest
                profitShareReserves = profitShareReserves + realizedCash*profitShare;
                
                insurersProfit = insurersProfit + realizedCash*(1-profitShare);
                [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);
                sT = sT - realizedCash*(1-profitShare)*weightEquity;
                bondsHTMValue = bondsHTMValue - realizedCash*(1-profitShare)*weightBonds;
                [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);
 
                realizedCash = 0;
            end
        
       %the realized cash is not enough to cover g
        else
            notCovered = guaranteedPayment - realizedCash;
            realizedCash = 0;
            
            %As available is large enough, the rest can be realized from KRF
            KRF = portfolio - (reserves + TA + profitShareReserves);
        end     
        
        if (KRF > KRFmax)
            realizedCash = KRF-KRFmax;

            %Realize this amount from the portfolio 

            if (TA + realizedCash < TAoldMaxValue)    
                TA = TA + realizedCash;

            %Else: Fill TA up to TAupperValue, and profit share the rest
            else
                realizedCash = realizedCash - (TAoldMaxValue - TA);
                TA = TAoldMaxValue;

                %Profit sharing the rest
                profitShareReserves = profitShareReserves + realizedCash*profitShare;
                insurersProfit = insurersProfit + realizedCash*(1-profitShare);
                [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);
                sT = sT - realizedCash*(1-profitShare)*weightEquity;
                bondsHTMValue = bondsHTMValue - realizedCash*(1-profitShare)*weightBonds;
                [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);
                realizedCash = 0;
            end
        end
    end
    
    KRF = portfolio - (reserves + TA + profitShareReserves);
    
    %reduce TA according to max value in years of payout
    TAtoPayment = 0;
    if (TA > TAupperValue)
        TAtoPayment = TAoldMaxValue - TAupperValue;
        TA = TAupperValue;
    end
    
    
    
    %Storing values
     if storeValues
        insurerProfitStored(year) = insurersProfit;
        TAStored(year) = TA;
        reservesStored(year) = reserves;
        profitShareReservesStored(year) = profitShareReserves;
        KRFStored(year) = KRF;
    end
    
    if storeValuesMean
        testingTAKRF_KRF(sim, year) = KRF;
        testingTAKRF_KRFMax(sim,year) = KRFmax;
        testingTAKRF_TA(sim, year) = TA;
        testingTAKRF_TAMax(sim, year) = TAupperValue;
    end
    
end

function [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT)
    portfolio = sT + bondsHTMValue;
    weightEquity = sT/portfolio;
    weightBonds = bondsHTMValue/portfolio;
end

function KRF = updateKRF(portfolio,TA,reserves,profitShareReserves)
    KRF = portfolio - (reserves + TA + profitShareReserves);
end

