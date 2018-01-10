function prmOpt = roe_prm_opt(annAvgEraPre, annAvgEraVapSat, eraGpGradMag, wgtLagLr, nPts)


alpha1 = cell(4,1);
[ alpha1{[1,2]} ] = deal(10.^(linspace( -5, -3, nPts)));
[ alpha1{3} ]     = deal(10.^(linspace( -3, -1, nPts)));
[ alpha1{4} ]     = deal(10.^(linspace( -2,  2, nPts)));

[prmOut1, maeRes1, biasRes1] = roe_pre_empirical_v2(annAvgEraVapSat, eraGpGradMag, wgtLagLr, 'cal', alpha1, annAvgEraPre);

[minMae1, indPrm1] = min(maeRes1);

%Display results
%MAE with dem slope (no wind direction): -0.8074
disp(['The first round of Roe orographic precipitation optimziation ' ...
    'has occured. The MAE is ' num2str(minMae1) ' and bias is ' num2str(biasRes1(indPrm1)) '.' char(10) ...
    'This occurred for the paramter set = ' num2str(prmOut1(indPrm1,1)) ', '...
    num2str(prmOut1(indPrm1,2)) ', ' num2str(prmOut1(indPrm1,3)) ', ' ...
    num2str(prmOut1(indPrm1,4)) '.']);

prmOpt1 = prmOut1(indPrm1,:);

%Second rond of optimization (select parameter space based on
%optimum from first round)
order1 = order(prmOpt1);
alpha2 = cell(4,1);
alpha2{1} = linspace( 10^(order1(1)-1), 10^(order1(1)+1), nPts);
alpha2{2} = linspace( 10^(order1(2)-1), 10^(order1(2)+1), nPts);
alpha2{3} = linspace( 10^(order1(3)-1), 10^(order1(3)+1), nPts);
alpha2{4} = linspace( 10^(order1(4)-1), 10^(order1(4)+1), nPts);


[prmOut2, maeRes2, biasRes2] = roe_pre_empirical_v2(annAvgEraVapSat, eraGpGradMag, wgtLagLr, 'cal', alpha2, annAvgEraPre);

[minMae2, indPrm2] = min(maeRes2);

prmOpt2 = prmOut2(indPrm2,:);

%Display results
%MAE with dem slope (no wind direction): -0.8074
disp(['The second round of Roe orographic precipitation optimziation ' ...
    'has occured. The MAE is ' num2str(minMae2) ' and bias is ' num2str(biasRes2(indPrm2)) '.' char(10) ...
    'This occurred for the paramter set = ' num2str(prmOut2(indPrm2,1)) ', '...
    num2str(prmOut2(indPrm2,2)) ', ' num2str(prmOut2(indPrm2,3)) ', ' ...
    num2str(prmOut2(indPrm2,4)) '.']);


%Third round of optimization:        
indPrm2(1) = find(alpha2{1} == prmOpt2(1));
indPrm2(2) = find(alpha2{2} == prmOpt2(2));
indPrm2(3) = find(alpha2{3} == prmOpt2(3));
indPrm2(4) = find(alpha2{4} == prmOpt2(4));

indPrm2Bnds = [indPrm2(:)-1, indPrm2(:)+1];
if indPrm2Bnds(1,1) == 0
    indPrm2Bnds(1,1) = 1;
    warning('EraHimalayaDs:param11', 'Optimized parameter 1 is at the lower bound. This should be avoided.');
elseif indPrm2Bnds(1,2) == nPts + 1
    indPrm2Bnds(1,2) = nPts;
    warning('EraHimalayaDs:param1End', 'Optimized parameter 1 is at the upper bound. This should be avoided.');
end
if indPrm2Bnds(2,1) == 0
    indPrm2Bnds(2,1) = 1;
    warning('EraHimalayaDs:param11', 'Optimized parameter 2 is at the lower bound. This should be avoided.');
elseif indPrm2Bnds(2,2) == nPts + 1
    indPrm2Bnds(2,2) = nPts;
    warning('EraHimalayaDs:param1End', 'Optimized parameter 2 is at the upper bound. This should be avoided.');
end
if indPrm2Bnds(3,1) == 0
    indPrm2Bnds(3,1) = 1;
    warning('EraHimalayaDs:param11', 'Optimized parameter 3 is at the lower bound. This should be avoided.');
elseif indPrm2Bnds(3,2) == nPts + 1
    indPrm2Bnds(3,2) = nPts;
    warning('EraHimalayaDs:param1End', 'Optimized parameter 3 is at the upper bound. This should be avoided.');
end
if indPrm2Bnds(4,1) == 0
    indPrm2Bnds(4,1) = 1;
    warning('EraHimalayaDs:param11', 'Optimized parameter 4 is at the lower bound. This should be avoided.');
elseif indPrm2Bnds(4,2) == nPts + 1
    indPrm2Bnds(4,2) = nPts;
    warning('EraHimalayaDs:param1End', 'Optimized parameter 4 is at the upper bound. This should be avoided.');
end

alpha3 = cell(4,1);
alpha3{1} = linspace( alpha2{1}(indPrm2Bnds(1,1)), alpha2{1}(indPrm2Bnds(1,2)), nPts);
alpha3{2} = linspace( alpha2{2}(indPrm2Bnds(2,1)), alpha2{2}(indPrm2Bnds(2,2)), nPts);
alpha3{3} = linspace( alpha2{3}(indPrm2Bnds(3,1)), alpha2{3}(indPrm2Bnds(3,2)), nPts);
alpha3{4} = linspace( alpha2{4}(indPrm2Bnds(4,1)), alpha2{4}(indPrm2Bnds(4,2)), nPts);

[prmOut3, maeRes3, biasRes3] = roe_pre_empirical_v2(annAvgEraVapSat, eraGpGradMag, wgtLagLr, 'cal', alpha3, annAvgEraPre);

[minMae3, indPrm3] = min(maeRes3);

% prmOpt3 = prmOut3(indPrm3,:);

%Validation run of Roe precipitation model
prmOpt = [{prmOut3(indPrm3,1)}, {prmOut3(indPrm3,2)}, {prmOut3(indPrm3,3)}, {prmOut3(indPrm3,4)}];

%Display results
%MAE with dem slope (no wind direction): -0.8074
disp(['The second round of Roe orographic precipitation optimziation ' ...
    'has occured. The MAE is ' num2str(minMae3) ' and bias is ' num2str(biasRes3(indPrm3)) '.' char(10) ...
    'This occurred for the paramter set = ' num2str(prmOut3(indPrm3,1)) ', '...
    num2str(prmOut3(indPrm3,2)) ', ' num2str(prmOut3(indPrm3,3)) ', ' ...
    num2str(prmOut3(indPrm3,4)) '.']);