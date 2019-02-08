testProfitShare = 0;
testProfitShareDiscounted = 0;


%////////////////////////////////////////////////////////////////
%///////////////////////////Simulation data//////////////////////
%////////////////////////////////////////////////////////////////
numSims = 100;
numIntervals = 4;
storeValues = false;
storeValuesMean = false;
 
initProfitShare = 0.8;
initWeightEquity = 0.2;
initTAmax = [0.10];
initKRFmax = 0.05;
g = 0.035;               %guaranteed rate
 
includeConfidenceIntervals = true; 
reduceKRFInPayoutYears = false;
rebalancing = true;
 
outputFile = '\\sambaad.stud.ntnu.no\erlengs\Documents\Fripoliser\Excel\newTesting.xlsx';
 
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

insurerProfitTestStored = zeros([1 numYears]);
insurerProfitTestStoredFullMatrix = zeros([numSims numYears]);
insurerProfitTestStoredIns = zeros([1 numYears]);
insurerProfitTestStoredPor = zeros([1 numYears]);

customerProfitTestStored = zeros([1 numYears]);
totalProfitTestStored = zeros([1 numYears]);

rateStored = zeros([numYears*numIntervals+1 numSims]);
 
%////////////////////////////////////////////////////////////////
%//////////////////////////Interest Rate/////////////////////////
%////////////////////////////////////////////////////////////////
zeroRate = [0.01247,0.01436,0.01538,0.01627,0.01706,0.01781,0.01852,...
    0.01917,0.01976,0.02026,0.0208,0.02137,0.02196,0.02254,0.02312,0.02368,...
    0.02422,0.02475,0.02525,0.02573,0.0262,0.02664,0.02706,0.02747,0.02785,...
    0.02822,0.02857,0.02891,0.02923,0.02953,0.02983,0.03011,0.03038,0.03063,...
    0.03088,0.03111,0.03134,0.03155,0.03176,0.03196,0.03215,0.03233,0.03251,...
    0.03268,0.03284,0.033,0.03315,0.03329,0.03343,0.03357,0.0337,0.03383,....
    0.03395,0.03406,0.03418,0.03429,0.03439,0.0345,0.0346,0.03469,0.03479,...
    0.03488,0.03496,0.03505,0.03513,0.03521,0.03529,0.03537,0.03544,0.03551,....
    0.03558,0.03565,0.03571,0.03578,0.03584,0.0359,0.03596,0.03602,0.03607,...
    0.03613,0.03618,0.03624,0.03629,0.03634,0.03639,0.03643,0.03648,0.03653,...
    0.03657,0.03661,0.03666,0.0367,0.03674,0.03678,0.03682,0.03686,0.03689,...
    0.03693,0.03697,0.037,0.03704,0.03707,0.0371,0.03714,0.03717,0.0372,...
    0.03723,0.03726,0.03729]';
 
%Converting the zeroRates into compounded rates
for x=1:size(zeroRate)
   compoundRate = log(1+zeroRate(x));
   zeroRate(x) = compoundRate;
end
 
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

%////////////////////////////////////////////////////////////////
%/////////////////////////Control Variate////////////////////////
%////////////////////////////////////////////////////////////////
endpricesStock = zeros([numSims 1]);
endpricesBond = zeros([numSims 1]);
endInsurersProfit = zeros([numSims 1]);
endCustomerProfit = zeros([numSims 1]);

 
tic
for testing = 1:size(testInput,2)
    disp("Simulation run: " + testing);
    randn('seed', 10);
    
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
    
    [reservesInit, realPaymentsInit] = calcReservesMen(yearsToRetirement,promisedPayment,... 
        g, pensionAge, deathAge, lastPaymentAge, realPayments);
    
for sim = 1:numSims
    testingTAKRF_ProfitShares = 0;
    testingTAKRF_ToCoverG = 0;
    testingTAKRF_NegativeReturn = 0;
    testingTAKRF_Rate = 0;
 
    promisedPayment = 100;
    
    reserves = reservesInit;
    realPayments = realPaymentsInit;
 
    mortalityRates = [];
    numYears = yearsToRetirement + yearsFromRetirementToDead;
    
    
    profitShare = initProfitShare;
    
    %Set initial portfolio
    KRF = 0.0*reserves;
    portfolio = reserves + KRF;
    KRFmaxPercentage = initKRFmax;    %maximum percentage of the reserves and proftitShare reserves
    KRFminPercentage = 0;       %minimum percentage of the reserves and profitShare reserves
    KRFmaxConst = 0;
    
    %Set initial asset mix
    weightEquity = initWeightEquity;
    initialEquityPercentage = weightEquity;
    equityMinLimit = 0;
    equityMaxLimit = 1;
    
    if (rebalancing)
        equityMinLimit = weightEquity - 0.10;
        equityMaxLimit = weightEquity + 0.10;  
    end
    
    totalEquityBought = 0;
    
    weightBonds = 1 - weightEquity;
        
    %Create stochastic z values to use in HullWhite and BlackScholes
    z = calcZ(numYears*numIntervals,numIntervals);
    
    %Simulate interest rate 
    hwInit = simulateInterestRate(numYears, numIntervals, z,alpha,sigma,zeroRate);
    
    %for testing with guaranteed rate
    %hwInit = simulateGuaranteed(g, numYears, numIntervals);
    
    bondHTMYears = [1,3,5,10];
    initHTMYears = [1,3,5,10];
    bondsHTMValue = weightBonds*portfolio;
    bondHTMWeights = [0.25 0.15 0.10 0.5];
    
    keySet = {1,3,5,10};
    valueSet = [2 3 4 5];
    bondDict = containers.Map(keySet,valueSet);
    
    simBondHTMReturn = zeros([1 4]);
    
    for i=1:4
        simBondHTMReturn(i) = hwInit(numIntervals*bondHTMYears(i)+1, (i+1));
    end
    
    %simBondHTMReturn = [0.01172 0.01640 0.01901 0.02221];     
    %simBondHTMReturn = [0.035 0.035 0.035 0.035];
   
    profitShareReserves = 0;  %value added to reserves due to profit sharing
    g_addReserves = 0;       %guaranteed rate on additional reserves   
    
  
    TA = 0.0*reservesInit;
    TAlowerLimit = 0;
    TAupperLimit = initTAmax(testing);
 
    insurersProfit = 0;
    customerPayment = 0;
    
    %Black & Scholes
    s0 = portfolio * weightEquity;
    b0 = portfolio * weightBonds;
    bT = b0;
    yearlyDiv = 0.0;
    div = yearlyDiv / numIntervals; %Dividend will be handled manually
    vol = .25;
    stockPrices = zeros([numIntervals numYears]);
    stockTesting = zeros([numIntervals numYears]);
    sT = s0;
    sT_testing = s0;

    rfDisk = hwInit(:,1,:);
    
    %For each year
    for year = 1: numYears
        %disp("YEAR: " + year);
        if (year > yearsToRetirement && year <= (numYears))
            %Finding the payments to the customer
            promisedPayment = realPayments(year-yearsToRetirement);
            
            promisedProfit = profitShareReserves/(numYears - year + 1);
            
            %Reducing the reserves and profitShareReserves by the
            %corresponding payment.
            reserves = reserves - promisedPayment;
            profitShareReserves = profitShareReserves - promisedProfit;
            
            %Testing the difference between the nominal and market value of
            %the payments from profit sharing. 
            testProfitShare = testProfitShare + promisedProfit;
            testProfitShareDiscounted = testProfitShareDiscounted + promisedProfit* ...
                exp(-calcBankAccountRate(rfDisk(1:(year)*numIntervals),numIntervals));
            
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
        if (storeValues)
            prePortfolioStored(year) = sT + bondsHTMValue;
        end
        
        %For each month/quartal
        for interval = 1:numIntervals
            
            %Getting the current risk free rate. 
            rf_current = rfDisk(interval + (year-1)*numIntervals);
            
            %Correlated z value. 
            z_temp = randn*sqrt(1/numIntervals);
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
            
            %Increasing the insurersProfit and customerPayment by the
            %risk free rate (for discounting purposes)
            insurersProfit = insurersProfit*exp(rf_current/numIntervals);
            customerPayment = customerPayment*exp(rf_current/numIntervals);

        end

        %Reduce stock price with dividends
        sT = sT - divPayments;
        
        %Portfolio update
        portfolio = bondsHTMValue + sT;
               
        %Increasing the insurer's profit by the risk-free rate.
        testingTAKRF_Rate = testingTAKRF_Rate + insurersProfit*rf_current;

        %disp("----");
        %disp("Year " + year);      
        
        %___________________________________________________________
        %BOND RETURNS
        %___________________________________________________________
        
        [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);

        [newEquityValue, equitySold, simBondHTMReturn, bondHTMWeights, ... 
            bondsHTMValue, bondHTMYears,  newEquityStored, ... 
            equitySoldStored, equityAddedStored, equityBought, couponPayment] = calcBondHTMReturn2(storeValues, bondHTMYears, hwInit,...
            year, numIntervals, simBondHTMReturn, bondHTMWeights, bondsHTMValue, ...
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
        
        bT = bT*exp(bookReturn);

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
               
        %Running buffer strategy
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
        
        [weightBonds, weightEquity,portfolio] = updateWeights(portfolio,bondsHTMValue,sT);
        
        if (year > yearsToRetirement && year < (numYears))
            customerPayment = customerPayment + TAtoPayment; 
            
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
        
        bankAccount = exp(-calcBankAccountRate(rfDisk(1:(year)*numIntervals),numIntervals));
        
        %Martingale testing
        martingaleEquity(year) =  martingaleEquity(year) + sT_testing*...
            bankAccount;
        
        martingaleBonds(year) =  martingaleBonds(year) + bT*...
            bankAccount;

        if storeValuesMean
            testingTAKRF_TotalProfitShares(sim,year) = insurersProfit * bankAccount;
        end
        
        %for testing the yearly development in total value for customer and
        %insurer
        
        insurerProfitTestStored(year) = insurerProfitTestStored(year) + ...  
            (insurersProfit + portfolio + customerPayment)...
            * bankAccount;
        
        if (includeConfidenceIntervals)
            insurerProfitTestStoredFullMatrix(sim,year) = insurerProfitTestStoredFullMatrix(sim,year) + ...  
                (insurersProfit + portfolio + customerPayment)* bankAccount;
        end   
    end
 
    %Print values from one simulation to excel
    if sim == 1 && storeValues
        filename = outputFile;
        sheet = 1;
        xlswrite(filename,stockPrices,sheet,'C3');
        xlswrite(filename,equityReturnStored,sheet,'C8');
        xlswrite(filename,bondReturnsStored,sheet,'C9');
        xlswrite(filename,weightsEquityStored,sheet,'C11');
        xlswrite(filename,totalReturnStored, sheet,'C14');
        %xlswrite(filename,unRealizedStored,sheet,'C14');
        %xlswrite(filename,preBookCashStored,sheet,'C15');
        
        xlswrite(filename,realizedCashStored,sheet,'C21');
        xlswrite(filename,inputKRFstored,sheet,'C22')
        xlswrite(filename,couponStored,sheet,'C23');
        xlswrite(filename,prePortfolioStored,sheet,'C27')
        xlswrite(filename,portfolioStored,sheet,'C28')
        
        xlswrite(filename,reservesStored,sheet,'C31'); 
        
        xlswrite(filename,TAStored,sheet,'C32');
        xlswrite(filename,insurerProfitStored,sheet,'C33');
        xlswrite(filename,profitShareReservesStored,sheet,'C35');
        xlswrite(filename,KRFStored,sheet,'C36'); 
        %xlswrite(filename,promisedProfitShareStored,sheet,'C25');
        xlswrite(filename,customerPayments,sheet,'C39');

        
        %xlswrite(filename,portfolioStored,sheet,'C31');
        xlswrite(filename,KRFMaxStored,sheet,'C43');
        xlswrite(filename,newEquityStored,sheet,'C45');
        xlswrite(filename,equitySoldStored,sheet,'C46');
        xlswrite(filename,equityAddedStored,sheet,'C47');
        disp("Simulering skrevet til Excel");
    end
    
    KRF = updateKRF(portfolio,TA,reserves,profitShareReserves); 
    
    bankAccount = exp(-calcBankAccountRate(rfDisk(1:(year)*numIntervals),numIntervals));
    
    %Discount with the bank account rate
    PVinsurersProfit = (insurersProfit + KRF*(1-profitShare)) * bankAccount;
    PVcustomerProfit = (customerPayment + KRF*(profitShare))* bankAccount;
    PVinsurersProfitStored(sim) = PVinsurersProfit;
    
    %Needed for control variates.
    endpricesStock(sim) = sT_testing * bankAccount;
    endpricesBond(sim) = bankAccount;
    
    endInsurersProfit(sim) = (insurersProfit + KRF*(1-profitShare)) * bankAccount;
    endCustomerProfit(sim) = (customerPayment + KRF*profitShare) * bankAccount;
    %///////////////////////////
    
    sT = s0;
    
    insurersProfitProfitShare(testing) = insurersProfitProfitShare(testing) + PVinsurersProfit;
    customerProfitProfitShare(testing) = customerProfitProfitShare(testing) + PVcustomerProfit;
    
    PVinsurersProfit = 0;
    PVcustomerProfit = 0;
    
    fprintf('Simulation %i of %i finished \n',sim,numSims); 
end
%////////////////////////////////////////////////////////////////////
%////////////////////////Control Variate/////////////////////////////
%////////////////////////////////////////////////////////////////////

covSI = (1/(numSims-1)) * sum((endpricesStock - mean(endpricesStock)).*(endInsurersProfit - mean(endInsurersProfit)));
covBI = (1/(numSims-1)) * sum((endInsurersProfit - mean(endInsurersProfit)).*(endpricesBond - mean(endpricesBond)));
covSC = (1/(numSims-1)) * sum((endpricesStock - mean(endpricesStock)).*(endCustomerProfit - mean(endCustomerProfit)));
covBC = (1/(numSims-1)) * sum((endCustomerProfit - mean(endCustomerProfit)).*(endpricesBond - mean(endpricesBond)));
expectedBond = exp(-sum(zeroRate(1:54)));

%Insurer
betaCoefficientI = inv(cov(endpricesStock,endpricesBond)) * [covSI; covBI];

insurersProfitControlVariate = mean(endInsurersProfit) - ...
   ([mean(endpricesStock) mean(endpricesBond)] - [s0 expectedBond])* betaCoefficientI;

%Customer
betaCoefficientC = inv(cov(endpricesStock,endpricesBond)) * [covSC; covBC];

customerProfitControlVariate = mean(endCustomerProfit) - ...
   ([mean(endpricesStock) mean(endpricesBond)] - [s0 expectedBond])* betaCoefficientC;

sumVariate = insurersProfitControlVariate + customerProfitControlVariate;

disp("-----");
disp("Values (with control variate)")
disp("Insurer: " + insurersProfitControlVariate);
disp("Policyholder: " + customerProfitControlVariate);
disp("Sum: " + sumVariate);
disp(" ");

R2I = (transpose([covSI; covBI])*inv(cov(endpricesStock,endpricesBond))*[covSI; covBI])/var(endInsurersProfit);
R2C = (transpose([covSC; covBC])*inv(cov(endpricesStock,endpricesBond))*[covSC; covBC])/var(endCustomerProfit);

VarianceInsurer = ((numSims - 2)/(numSims - 4)) * (1-R2I) * var(endInsurersProfit);
VarianceCustomer = ((numSims - 2)/(numSims - 4)) * (1-R2C) * var(endCustomerProfit);

varianceReductionInsurer = VarianceInsurer/var(endInsurersProfit);
varianceReductionCustomer = VarianceCustomer/var(endCustomerProfit);

disp("Variance reduction insurer: " + varianceReductionInsurer); 
disp("Variance reduction customer: " + varianceReductionCustomer); 
disp("-----");

%//////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////

insurersProfitProfitShare(testing) = insurersProfitProfitShare(testing)/numSims;
customerProfitProfitShare(testing) = customerProfitProfitShare(testing)/numSims;

insurerProfitTestStored = insurerProfitTestStored/numSims;


%Calculating confidence intervals
if (includeConfidenceIntervals)
    insurerProfitTestStoredFullMatrixPercentiles = prctile(insurerProfitTestStoredFullMatrix,[1 5 10 90 95 99],1);
    insurerProfitTestStoredFullMatrixMean = mean(insurerProfitTestStoredFullMatrix);
    insurerProfitTestStoredFullMatrixStd = std(insurerProfitTestStoredFullMatrix);
    zConf = norminv(1 - 0.05/2);
    confidenceIntervalsLower = insurerProfitTestStoredFullMatrixMean - zConf*insurerProfitTestStoredFullMatrixStd/sqrt(numSims);
    confidenceIntervalsUpper = insurerProfitTestStoredFullMatrixMean + zConf*insurerProfitTestStoredFullMatrixStd/sqrt(numSims);
end
%--------------------------------

disp("Values (without control variate)");
disp("Insurer: " + insurersProfitProfitShare);
disp("Policyholder: " + customerProfitProfitShare);
total = insurersProfitProfitShare + customerProfitProfitShare;
disp("Sum: " + total);
disp("-----");

disp("Test Values: ")
disp(testInput);
disp("KRF: " + initKRFmax)
 
martingaleBonds = martingaleBonds / numSims;
martingaleEquity = martingaleEquity / numSims;

testProfitShareDiscounted = testProfitShareDiscounted / numSims;
testProfitShare = testProfitShare / numSims;

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
           %rateIndex = (year)*interval + 1;
           %rateIndex = (year-1)*interval + 1 + buyYear*interval;
           rateIndex = (year-1)*interval + 1 + min(buyYear,4)*interval;
                
           if (buyYear == 1)
                interestColumn = 2;
           elseif (buyYear == 3)
                interestColumn = 3;
           elseif (buyYear == 5)
                interestColumn = 4;
           elseif (buyYear == 10)
                interestColumn = 5;
           end 
           
           currentBondReturns(bond) = hw(rateIndex,interestColumn);
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

%Function used to find the correct column in the HW-matrix
function columnIndex = findColumn(buyYear)
    if (buyYear == 1)
        columnIndex = 2;
    elseif (buyYear == 3)
        columnIndex = 3;
    elseif (buyYear == 5)
        columnIndex = 4;
    elseif (buyYear == 10)
        columnIndex = 5;
    end 
end 
 
%Functions used for simulating HW and BS
function zArray = calcZ(number,numIntervals)
    temp = zeros([number+10*numIntervals 1]);
    for x = 1:number
        temp(x) = randn;
    end
    zArray = temp;
end

function rateArray = simulateInterestRate(years, interval, Z, alpha,sigma,zeroRate)
    rateArray = zeros([years*interval 5]);

    nPeriods = interval*years + 10*interval;
    deltaTime = 1/interval;
    count = 1;
    countIntv = 1;
    yPrev = 0;
    
    for t=1:nPeriods
        B0t = 1 - exp(-alpha * (t - 0));
        B02dt = 1 - exp(-alpha * (2*deltaTime -0));
        B0dt = 1 - exp(-alpha * (deltaTime -0));


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
        %y = exp(-alpha*(deltaTime))*yPrev + sigma*B0dt*Z(t);
        r = a + y;
        rateArray(t,1) = r;
        rateArray(t,2) = r;
        rateArray(t,3) = r;
        rateArray(t,4) = r;
        rateArray(t,5) = r;

        yPrev = y;
    end
end

%for testing with guaranteed interest rate
function rateArrayGuaranteed = simulateGuaranteed(g, years, interval)
        rateArrayGuaranteed = zeros([years*interval 5]);
        
         nPeriods = interval*years + 50;
        
        for t=1:nPeriods
            rateArrayGuaranteed(t,1) = g;
            rateArrayGuaranteed(t,2) = g;
            rateArrayGuaranteed(t,3) = g;
            rateArrayGuaranteed(t,4) = g;
            rateArrayGuaranteed(t,5) = g;
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
    sT = sT + realizedCash;
    %bondsHTMValue = bondsHTMValue + realizedCash*weightBonds;
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
        
        %If the return is negative, we need to inject and inject everything
        %into stocks 
        sT = sT - available;
        %bondsHTMValue = bondsHTMValue - available*weightBonds;
        
        
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
            sT = sT + neededCash;
            
            %bondsHTMValue = bondsHTMValue + neededCash*weightBonds;
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
            
            %If the return is negative, we need to inject and inject everything
            %into stocks 
            sT = sT + neededCash;
            %bondsHTMValue = bondsHTMValue + neededCash*weightBonds;
            
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
                
                %We realize from stocks in case of profit share
                sT = sT - realizedCash*(1-profitShare);
                
                %OLD:
                %bondsHTMValue = bondsHTMValue - realizedCash*(1-profitShare)*weightBonds;
                
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
                sT = sT - realizedCash*(1-profitShare);
                
                %bondsHTMValue = bondsHTMValue - realizedCash*(1-profitShare)*weightBonds;
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

%test