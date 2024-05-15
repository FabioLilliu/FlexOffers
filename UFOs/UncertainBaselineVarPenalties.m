format shortG
yalmip('solver','glpk');
sett = sdpsettings('verbose',0);

%Program definition
cloud = 1; %set to 1 for using it in cloud
par = 8; %number of iterations going in parallel

if cloud == 0
    load('C:\Users\IQ24VR\OneDrive - Aalborg Universitet\Skrivebord\Postdoc\MATLAB files\M.mat');
    load('C:\Users\IQ24VR\OneDrive - Aalborg Universitet\Skrivebord\Postdoc\MATLAB files\N.mat');
elseif cloud == 1
    load('/home/ubuntu/MATLABfiles/M.mat');
    load('/home/ubuntu/MATLABfiles/N.mat');
end

for pt1 = 100:-10:60
parfor parallel = 1:par %change to "for" or "parfor" 

T = 6; %generation/planning horizon 
    
%Parameters of the algorithm
NumSamples = 3 ; %if needed, type of FO/DFO used
NumBatteries = 1;

%Parameters of the battery
InitialSoC = 0; %initial state of charge
L = sqrt(0.9); %loss of the battery (square root of roundtrip efficiency)
Qmin = 0; %minimum possible state of charge
Qmax = 14;%maximum possible state of charge
Pmin = -5; %maximum discharging power
Pmax = 5; %maximum charging power
decay = 1; %decay at each time unit
percentage = 100; %percentage of flexibility exploited by DFO
pt = pt1/100; %probability threshold
gr = 100; %granularity

days = ceil(365/par);%number of days for iteration
Length = floor(days*24/T); %number of cycles needed

RMatrix = repelem(0,19,Length); %initializing result matrix
SoC0 = repelem(0,NumBatteries,Length+1); %initializing "continuity array", 
          %it describes the state of charge at the beginning of next cycle
SoC0(:,1) = InitialSoC'; %initial state of charge
daysahead = 0; %how many days ahead we want to start; mainly for debugging
Unit = eye(T);

for counter = 1:Length
    Q0 = SoC0(:,counter); %initial state of charge for the cycle
    spos = T*(counter - 1) + T*Length*(parallel-1);        
    spotprices = M(1+spos+daysahead:T+spos+daysahead,1:1)'; %spot prices
    imbalanceprices = N(1+spos+daysahead:T+spos+daysahead,1:1)'; %imbalance prices

    %%%%% Instantiate LTI models, one for each battery
    
    %First kind
    k1 = 1;
    A = decay;
    B3 = [L,1/L]; %FOR DFOS ONLY
    f = [0];
    C = [0];
    D1 = [1,1]; %FOR DFOS ONLY
    PminV = [0;L*Pmin]; %FOR DFOS ONLY
    PmaxV = [Pmax;0]; %FOR DFOS ONLY
    lti = LTISystem('A', A, 'B', B3, 'C', C, 'D', D1); %FOR DFOS ONLY
    lti.u.min = PminV;
    lti.u.max = PmaxV;
    lti.x.min = Qmin;
    lti.x.max = Qmax;
    lti.initialize(Q0);
    data1 = struct();
    data1.InitialSoC = Q0;
    data1.L = L;
    data1.Pmin = Pmin;
    data1.Pmax = Pmax;
    data1.Qmin = Qmin;
    data1.Qmax = Qmax;
    data1.percentage = 100;
    data1.pt = pt;

    P = Polyhedron('lb', [Q0], 'ub', [Q0]);
    lti.x.with('initialSet');
    lti.x.initialSet = P;

    lti.instantiate(T);
    F1 = FlexSystem(lti);
    
    %Calculate profit - LTI model
    
    constraints = F1.getAllConstraints;
    costbattery = spotprices(1:T) * F1.getEnergyVars()';
    optimize(constraints,costbattery,sett);
    aggrcost = value(costbattery);
    EnergyValuesUnion = value(F1.getEnergyVars); %Save the energy values for calculating Q
    
    %Calculate profit - baseline
    
    %create constraints
    evr = sdpvar(2,T);
    Qa = sdpvar(T+1);
    Qa(1) = Q0;
    ctrs = [repelem(L*Pmin,T) <= evr(1,:) <= repelem(0,T); repelem(0,T) <= evr(2,:) <= repelem(Pmax,T)];
    pdist = (1-L)*Pmax;
    for t = 1:T
        Qa(t+1) = Qa(t) + L*evr(2,t) + evr(1,t)/L;
        QmaxU = Qmax + (1-pt)*pdist*floor((t-1)/2);
        ctrs = [ctrs; Qmin <= Qa(t+1) <= QmaxU];
    end
    basecost = spotprices(1:T) * transpose(evr(1,:)+evr(2,:));
    
    %optimize
    baseSolution = optimize(ctrs,basecost,sett);
    basecostValue = value(basecost);
    baseEnergyValue = value(evr(1,:)+evr(2,:));

    %imbalance calculation and disaggregation
    constraints = F1.model.constraints;
    Imbalance = imbalance_calculation_n(F1.getEnergyVars(),baseEnergyValue,imbalanceprices(1:T),constraints);
    baseEnergyImbalance = Imbalance(1);
    baseCostImbalance = Imbalance(2);
    baseTotalCost = basecostValue + baseCostImbalance;
    
    %Calculate profit - DFO
    
    %create constraints
    DFO = DFOSystem(F1,0.3,data1);
    DFOcost = spotprices(1:T) * transpose(DFO.getEnergyVars());
    
    %optimize
    DFOSolution = optimize(DFO.getAllConstraints(),DFOcost,sett);
    DFOcostValue = value(DFOcost);
    DFOEnergyValue = value(DFO.getEnergyVars);

    %imbalance calculation and disaggregation
    constraints = F1.model.constraints;
    Imbalance = imbalance_calculation_n(F1.getEnergyVars(),DFOEnergyValue,imbalanceprices(1:T),constraints);
    DFOEnergyImbalance = Imbalance(1);
    DFOCostImbalance = Imbalance(2);
    DFOTotalCost = DFOcostValue + DFOCostImbalance;
    
    %Calculate profit - ProbDFO
    
    %constraints
    PFO = DFOSystem(data1,T,pt,"prob");
    PFOcost = spotprices(1:T) * transpose(PFO.getEnergyVars());
    
    %optimization
    PFOSolution = optimize(PFO.getAllConstraints,PFOcost,sett);
    PFOcostValue = value(PFOcost);
    PFOEnergyValue = value(PFO.getEnergyVars);
    PFOEnergyValue1 = repelem(0,T);
    PFOEnergyValue1(1) = DFOEnergyValue(1);
    cont = [1];
    for t = 2:T
        kh = sum(PFOEnergyValue1(1:t-1));
        s = DFOEnergyValue(t);
        prob = prob_sum(data1,t,kh,s);
        if prob >= pt || imbalanceprices(t) < 0
            PFOEnergyValue1(t) = DFOEnergyValue(t);
        else
            [t1,t2] = find_thresholds(data1,t,kh,pt);
            if s < 0
                PFOEnergyValue1(t) = t1;
            else
                PFOEnergyValue1(t) = t2;
            end
        end
        if PFOEnergyValue1(t)*PFOEnergyValue1(t-1) >= 0.01 && t <= T-1
            cont = [cont,t];
        else
            if PFOEnergyValue1(t)*PFOEnergyValue1(t-1) >= 0.01 && t == T
                cont = [cont,t];
            end
            sc = size(cont);
            if sc(2) > 1
                v = sdpvar(1,sc(2));
                co = [sum(v) == sum(PFOEnergyValue1(cont))];
                if PFOEnergyValue1(t-1) >= 0.01
                    co = [co, 0 <= v <= Pmax];
                else
                    co = [co, L*Pmin <= v <= 0];
                end
                optimize(co,spotprices(cont)*v',sett);
                PFOEnergyValue1(cont) = value(v);
            end
            cont = [t];
        end
    end
    PFOcostValue1 = value(spotprices(1:T) * PFOEnergyValue1');
    
    %imbalance calculation and disaggregation
    
    %Only first step
    constraints = F1.model.constraints;
    ImbalanceP = imbalance_calculation_n(F1.getEnergyVars(),PFOEnergyValue,imbalanceprices(1:T),constraints);
    PFOEnergyImbalance = ImbalanceP(1);
    PFOCostImbalance = ImbalanceP(2);
    PFOTotalCost = PFOcostValue + PFOCostImbalance;
    
    %Second step
    ImbalanceP1 = imbalance_calculation_n(F1.getEnergyVars(),PFOEnergyValue1,imbalanceprices(1:T),constraints);
    PFOEnergyImbalance1 = ImbalanceP1(1);
    PFOCostImbalance1 = ImbalanceP1(2);
    PFOTotalCost1 = PFOcostValue1 + PFOCostImbalance1;
    
    %Calculate profit - DFOadjusted
    
    %create constraints
    DFOU = DFOSystem(F1,0.4,data1);
    DFOUcost = spotprices(1:T) * transpose(DFOU.getEnergyVars());
    
    %optimize
    DFOUSolution = optimize(DFOU.getAllConstraints(),DFOUcost,sett);
    DFOUcostValue = value(DFOUcost);
    DFOUEnergyValue = value(DFOU.getEnergyVars);

    %imbalance calculation and disaggregation
    constraints = F1.model.constraints;
    ImbalanceU = imbalance_calculation_n(F1.getEnergyVars(),DFOUEnergyValue,imbalanceprices(1:T),constraints);
    DFOUEnergyImbalance = ImbalanceU(1);
    DFOUCostImbalance = ImbalanceU(2);
    DFOUTotalCost = DFOUcostValue + DFOUCostImbalance;
    
    % DETERMINE SoC
    
    Q = repelem(0,NumBatteries,T+1);
    Q(:,1) = Q0;
    for t = 1:T
        for k = 1:NumBatteries
            Energy = EnergyValuesUnion;
            Q(k,t+1) = decay(k)*Q(k,t) + L(k)*max(Energy(t),0) + (1/L(k))*min(Energy(t),0);
        end
    end
    
    for k = 1:NumBatteries
        SoC0(k,counter+1) = max(0,min(Qmax(k),Q(k,T+1)));
    end
    
    %%%%% ADD RESULTS
    
    RMatrix(1,counter) = aggrcost;
    RMatrix(2,counter) = basecostValue;
    RMatrix(3,counter) = baseCostImbalance;
    RMatrix(4,counter) = baseTotalCost;     
    RMatrix(5,counter) = DFOUcostValue;
    RMatrix(6,counter) = DFOUCostImbalance;
    RMatrix(7,counter) = DFOUTotalCost;
    RMatrix(8,counter) = PFOcostValue1;
    RMatrix(9,counter) = PFOCostImbalance1;
    RMatrix(10,counter) = PFOTotalCost1;  
    
    disp('Counter is');
    disp(value(counter));
    disp('Parallel is');
    disp(value(parallel));
end
    %filename = sprintf('/home/ubuntu/imbalance_Q0is0_Tis6_Typeis0_03_daysare1to372part%d.xlsx',parallel);
    if cloud == 1
        filename = sprintf("/home/ubuntu/UncBaselineSwitching%d_%d_%d_%d.xlsx",T,NumBatteries,parallel,pt1);
    else
        filename = "C:\Users\Fab\Documents\Postdoc\Disaggregate"+T+"_"+NumBatteries+".xlsx";
    end
    writematrix(RMatrix,filename)

end
end