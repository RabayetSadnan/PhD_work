% Version: April 30, 2024
% Author: Rabayet Sadnan (c)
% -------------------------------------------------------------------------
% This is an LP OPF Solver for 123 Bus test system
% -------------------------------------------------------------------------
% Variables
% V_sub            : Substation voltage magnitude square
% CVR              : CVR factor for voltage dependent load
% Area             : 13 or 123 bus system
% v_max            : Upper limit of the nodal votlages
% v_min            : lower limit of the nodal votlages
% kV               : Base kV
% load_mult        : load multiplier
% gen_mult         : Dg generation multiplier/ PV irradiance
% S_mult           : The % of nominal kVA rating-- 1.2 means the rating of
%                    DG is 120% of the nominal real power rating
%                    e.g., 1kW DG has a rating of 1.2 kVA.
%
% Called Functions
% objfun           : The objective/cost function of the local OPF
% -------------------------------------------------------------------------
% This local solvers use the 3-ph OPF formulation that has been developed
% in <a href="https://ieeexplore.ieee.org/document/8598813">this Paper</a>.
% -------------------------------------------------------------------------

%%

clc;
clear all;

root_directory = "C:\Users\sadn725\OneDrive - PNNL\PNNL\AppDeconfliction\lindistflow_matlab_code\Test_system";
save_variable_directory = strcat(root_directory,"\DSS_validation\");

%% Variables:
V_sub = 1.05^2;
CVR = [0 0];
Area = 123;
v_max = 1.05^2;
v_min = 0.95^2;
kV = 4.16;
load_mult = 1;
gen_mult = 1;
S_mult = 1.2;


%% Input Formatting:

VsA = V_sub ;
VsB = V_sub ;
VsC = V_sub ;
CVR_Pe = CVR(1);                %%% CVR factor for P = 0.6
CVR_Qe = CVR(2);                  %%% CVR factor for Q = 3

s1 = strcat(root_directory,'\Area_',num2str(Area),'\branchdata.txt');
load (s1);

branch=sortrows(branchdata,2);

lineA = branch((find(branch(:,3)~=0)),1:2); % 3rd column for rAA
lineB = branch((find(branch(:,6)~=0)),1:2); % 6th column for rBB
lineC = branch((find(branch(:,8)~=0)),1:2); % 8th column for rCC

s1 = strcat(root_directory,'\Area_',num2str(Area),'\powerdata.txt');
load (s1);
powerdata = sortrows(powerdata,1);

global TableA TableB TableC;
global nb;
global TA TB TC;
global RAA RBB RCC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Calculation Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fb = branch(:,1);
tb = branch(:,2);
fbA = lineA(:,1);
tbA = lineA(:,2);
G = graph(fb,tb);
nbusA =  length(lineA) +1;
GA = graph(fbA,tbA);
VA_G = dfsearch(GA,1);


fbB = lineB(:,1);
tbB = lineB(:,2);

nbusB =  length(lineB) +1;
GB = graph(fbB,tbB);
VB_G = dfsearch(GB,1);

fbC = lineC(:,1);
tbC = lineC(:,2);
nbusC =  length(lineC) +1;
GC = graph(fbC,tbC);
VC_G = dfsearch(GC,1);

tnbA = length(fbA);
tnbB = length(fbB);
tnbC = length(fbC);

nb = size(powerdata,1);
bKVA = 1000;
bKV = kV/sqrt(3);
bZ = ((bKV)^2)*1000/bKVA;
rAA = ((branch(:,3))/bZ);
rAB = ((branch(:,4))/bZ);
rAC = ((branch(:,5))/bZ);
rBB = ((branch(:,6))/bZ);
rBC = ((branch(:,7))/bZ);
rCC = ((branch(:,8))/bZ);
xAA = ((branch(:,9))/bZ);
xAB = ((branch(:,10))/bZ);
xAC = ((branch(:,11))/bZ);
xBB = ((branch(:,12))/bZ);
xBC = ((branch(:,13))/bZ);
xCC = ((branch(:,14))/bZ);
PLA = (powerdata(:,2).*load_mult)/bKVA;           %%% Pload of phase A
PLB = (powerdata(:,4).*load_mult)/bKVA;           %%% Pload of phase B
PLC = (powerdata(:,6).*load_mult)/bKVA;           %%% Pload of phase C
QLA = (powerdata(:,3).*load_mult)/bKVA;           %%% Qload of phase A
QLB = (powerdata(:,5).*load_mult)/bKVA;           %%% Qload of phase B
QLC = (powerdata(:,7).*load_mult)/bKVA;           %%% Qload of phase C
QCA = powerdata(:,8)/bKVA;
QCB = powerdata(:,9)/bKVA;
QCC = powerdata(:,10)/bKVA;
PDERA = powerdata(:,11)/bKVA;         %%% PRated PU of phase A
PDERB = powerdata(:,12)/bKVA;         %%% PRated of phase B
PDERC = powerdata(:,13)/bKVA;          %%% PRated of phase C

PGA = (PDERA.*gen_mult);         %%% Pgen of phase A
PGB = (PDERB.*gen_mult);         %%% Pgen of phase B
PGC = (PDERC.*gen_mult);          %%% Pgen of phase C



%% DER Configuration:
SderA = S_mult*PDERA;
% find(Sder(:)~=0)--->   %DER connected BUS
S_DERA = SderA(find(SderA(:)~=0));   %in PU
P_DERA = PGA(find(SderA(:)~=0));   %in PU
l_bQA(:,1) = -sqrt((S_DERA.^2)-(P_DERA.^2));
u_bQA(:,1) = sqrt((S_DERA.^2)-(P_DERA.^2));

SderB = S_mult*PDERB;
% find(Sder(:)~=0)--->   %DER connected BUS
S_DERB = SderB(find(SderB(:)~=0));   %in PU
P_DERB = PGB(find(SderB(:)~=0));   %in PU
l_bQB(:,1) = -sqrt((S_DERB.^2)-(P_DERB.^2));
u_bQB(:,1) = sqrt((S_DERB.^2)-(P_DERB.^2));

SderC = S_mult*PDERC;
% find(Sder(:)~=0)--->   %DER connected BUS
S_DERC = SderC(find(SderC(:)~=0));   %in PU
P_DERC = PGC(find(SderC(:)~=0));   %in PU
l_bQC(:,1) = -sqrt((S_DERC.^2)-(P_DERC.^2));
u_bQC(:,1) = sqrt((S_DERC.^2)-(P_DERC.^2));

l_bQ = [l_bQA;l_bQB;l_bQC];
u_bQ = [u_bQA;u_bQB;u_bQC];


%%
T=dfsearch(G,1,'edgetonew');
TA=dfsearch(GA,1,'edgetonew');
TB=dfsearch(GB,1,'edgetonew');
TC=dfsearch(GC,1,'edgetonew');

RAA = zeros(nb);
XAA  = zeros(nb);
RAB = zeros(nb);
XAB  = zeros(nb);
RAC = zeros(nb);
XAC  = zeros(nb);
RBB = zeros(nb);
XBB  = zeros(nb);
RBC = zeros(nb);
XBC  = zeros(nb);
RCC = zeros(nb);
XCC  = zeros(nb);

for i = 1:(nb-1)
    RAA(fb(i), tb(i)) = rAA(i);
    RAA(tb(i) ,fb(i)) = RAA(fb(i), tb(i)) ;
    XAA(fb(i), tb(i))= xAA(i);
    XAA(tb(i) ,fb(i)) = XAA(fb(i), tb(i)) ;
    RAB(fb(i), tb(i)) = rAB(i);
    RAB(tb(i) ,fb(i)) = RAB(fb(i), tb(i)) ;
    XAB(fb(i), tb(i))=  xAB(i);
    XAB(tb(i) ,fb(i)) = XAB(fb(i), tb(i)) ;
    RAC(fb(i), tb(i)) = rAC(i);
    RAC(tb(i) ,fb(i)) = RAC(fb(i), tb(i)) ;
    XAC(fb(i), tb(i))=  xAC(i);
    XAC(tb(i) ,fb(i)) = XAC(fb(i), tb(i)) ;
    RBB(fb(i), tb(i)) = rBB(i);
    RBB(tb(i) ,fb(i)) = RBB(fb(i), tb(i)) ;
    XBB(fb(i), tb(i))= xBB(i);
    XBB(tb(i) ,fb(i)) = XBB(fb(i), tb(i)) ;
    RBC(fb(i), tb(i)) = rBC(i);
    RBC(tb(i) ,fb(i)) = RBC(fb(i), tb(i)) ;
    XBC(fb(i), tb(i))= xBC(i);
    XBC(tb(i) ,fb(i)) = XBC(fb(i), tb(i)) ;
    RCC(fb(i), tb(i)) = rCC(i);
    RCC(tb(i) ,fb(i)) = RCC(fb(i), tb(i)) ;
    XCC(fb(i), tb(i))= xCC(i);
    XCC(tb(i) ,fb(i)) = XCC(fb(i), tb(i)) ;
end

Ap = 1:nbusA-1;                          % defining the unknowns for phaseA
Aq = nbusA:2*(nbusA-1);
Av = 2*(nbusA):3*(nbusA-1)+1;

Bp = (3*(nbusA-1)+1)+1:(3*(nbusA-1)+1)+(nbusB-1);                    % defining the unknowns for phaseB
Bq = (3*(nbusA-1)+1)+nbusB :(3*(nbusA-1)+1) + 2*(nbusB-1);
Bv = (3*(nbusA-1)+1)+2*(nbusB):(3*(nbusA-1)+1) + 3*(nbusB-1) +1 ;


Cp = ((3*(nbusA-1)+1) + 3*(nbusB-1) +1)+1:((3*(nbusA-1)+1) + 3*(nbusB-1) +1) + (nbusC-1);             % defining the unknowns for phaseC
Cq = ((3*(nbusA-1)+1) + 3*(nbusB-1) +1) + nbusC :((3*(nbusA-1)+1) + 3*(nbusB-1) +1) + 2*(nbusC-1);
Cv = ((3*(nbusA-1)+1) + 3*(nbusB-1) +1)+2*(nbusC):((3*(nbusA-1)+1) + 3*(nbusB-1) +1) + 3*(nbusC-1) +1 ;


TableA = [TA(:,1) TA(:,2) Ap'  Aq'   Av'];
TableB = [TB(:,1) TB(:,2) Bp'  Bq'   Bv'];
TableC = [TC(:,1) TC(:,2) Cp'  Cq'   Cv'];

Da = 2*(nbusA)-1:3*(nbusA-1)+1;
Db = (3*(nbusA-1)+1)+2*(nbusB)-1:(3*(nbusA-1)+1) + 3*(nbusB-1) +1;
Dc = ((3*(nbusA-1)+1) + 3*(nbusB-1) +1)+2*(nbusC)-1:((3*(nbusA-1)+1) + 3*(nbusB-1) +1) + 3*(nbusC-1) +1 ;

VolttableA = Da';
VolttableB = Db';
VolttableC = Dc';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Aeq and Beq Formation%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
Aeq =[];
for i =2:nb

    CVR_P = CVR_Pe;
    CVR_Q = CVR_Qe;

    row = find(i == TA(:,1));
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
    if isempty(ParentA)
        Aeq = Aeq;
    else
        Poc = find(TA(ParentA,1) == TA(:,1));
        Aeq(ParentA,TableA(ParentA,3))= 1;
        if ~isempty(row)
            Aeq(ParentA,TableA(row(1),3))= -1;
            for j = 1:length(row)-1
                Aeq(ParentA, TableA(row(j+1),3)) =   - 1;
            end
        end
        Aeq(ParentA,TableA(ParentA,5))= -(CVR_P/2)*PLA(TableA(ParentA,2));
        Aeq(ParentA+(nbusA-1),TableA(ParentA,4))= 1;
        if ~isempty(row)
            Aeq(ParentA+(nbusA-1),(nbusA-1)+ TableA(row(1),3))= -1;
            for j = 1:length(row)-1
                Aeq(ParentA+(nbusA-1),(nbusA-1)+TableA(row(j+1),3)) =  -  1;
            end
        end
        Aeq(ParentA+(nbusA-1),(TableA(ParentA,5)))= -(CVR_Q/2)*QLA(TableA(ParentA,2));
        Aeq(ParentA+2*(nbusA-1),TableA(ParentA,5))= 1;
        Aeq(ParentA+2*(nbusA-1),VolttableA(Poc(1)))= -1;
        Aeq(ParentA+2*(nbusA-1),TableA(ParentA,3))= 2*(RAA(TA((ParentA),1),TA((ParentA),2)));
        Aeq(ParentA+2*(nbusA-1),TableA(ParentA,4))= 2*(XAA(TA((ParentA),1),TA((ParentA),2)));
        Aeq(ParentA+2*(nbusA-1),TableB(ParentB,3))= -RAB(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(XAB(TA((ParentA),1),TA((ParentA),2)));
        Aeq(ParentA+2*(nbusA-1),TableB(ParentB,4))= -XAB(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(RAB(TA((ParentA),1),TA((ParentA),2)));
        Aeq(ParentA+2*(nbusA-1),TableC(ParentC,3))= -RAC(TA((ParentA),1),TA((ParentA),2))- sqrt(3)*(XAC(TA((ParentA),1),TA((ParentA),2)));
        Aeq(ParentA+2*(nbusA-1),TableC(ParentC,4))= -XAC(TA((ParentA),1),TA((ParentA),2))+ sqrt(3)*(RAC(TA((ParentA),1),TA((ParentA),2)));
        beq(ParentA)= (1-(CVR_P/2))*PLA(TableA(ParentA,2))-PGA(TableA(ParentA,2));
        beq(ParentA+(nbusA-1)) =  (1-(CVR_Q/2))*QLA(TableA(ParentA,2))-QCA(TableA(ParentA,2));

    end
    Aeq(3*(nbusA-1)+1,VolttableA(1)) = 1;
    beq(3*(nbusA-1)+1) = VsA;

end


%%
for i =2:nb

    CVR_P = CVR_Pe;
    CVR_Q = CVR_Qe;

    row = find(i == TB(:,1));
    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
    if isempty(ParentB)
        Aeq = Aeq;
    else
        Poc = find(TB(ParentB,1) == TB(:,1));
        Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,3))= 1;

        if ~isempty(row)
            Aeq(3*(nbusA-1)+1+ParentB,TableB(row(1),3))= -1;
            for j = 1:length(row)-1
                Aeq(3*(nbusA-1)+1+ParentB, TableB(row(j+1),3)) =   - 1;
            end
        end

        Aeq(3*(nbusA-1)+1+ParentB,TableB(ParentB,5))= -(CVR_P/2)*PLB(TableB(ParentB,2));
        Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),TableB(ParentB,4))= 1;

        if ~isempty(row)
            Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+ TableB(row(1),3))= -1;
            for j = 1:length(row)-1
                Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(nbusB-1)+TableB(row(j+1),3)) =  -1;
            end
        end

        Aeq(3*(nbusA-1)+1+ParentB+(nbusB-1),(TableB(ParentB,5)))= -(CVR_Q/2)*QLB(TableB(ParentB,2));
        Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,5))= 1;
        Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),VolttableB(Poc(1)))= -1;
        Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,3))= 2*(RBB(TB((ParentB),1),TB((ParentB),2)));
        Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableB(ParentB,4))= 2*(XBB(TB((ParentB),1),TB((ParentB),2)));
        Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,3))= -RAB(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(XAB(TB((ParentB),1),TB((ParentB),2)));
        Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableA(ParentA,4))= -XAB(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(RAB(TB((ParentB),1),TB((ParentB),2)));
        Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,3))= -RBC(TB((ParentB),1),TB((ParentB),2))+ sqrt(3)*(XBC(TB((ParentB),1),TB((ParentB),2)));
        Aeq(3*(nbusA-1)+1+ParentB+2*(nbusB-1),TableC(ParentC,4))= -XBC(TB((ParentB),1),TB((ParentB),2))- sqrt(3)*(RBC(TB((ParentB),1),TB((ParentB),2)));
        beq(3*(nbusA-1)+1+ParentB)= (1-(CVR_P/2))*PLB(TableB(ParentB,2))-PGB(TableB(ParentB,2));
        beq(3*(nbusA-1)+1+ParentB+(nbusB-1)) = (1-(CVR_Q/2))*QLB(TableB(ParentB,2))-QCB(TableB(ParentB,2));
    end

Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1,VolttableB(1)) = 1;
beq(3*(nbusA-1)+1+3*(nbusB-1)+1) = VsB;
end


%%
for i =2:nb
    CVR_P = CVR_Pe;
    CVR_Q = CVR_Qe;

    row = find(i == TC(:,1));

    ParentA = find(i == TA(:,2));
    ParentB = find(i == TB(:,2));
    ParentC = find(i == TC(:,2));
    if isempty(ParentC)
        Aeq = Aeq;
    else
        Poc = find(TC(ParentC,1) == TC(:,1));
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,3))= 1;

        if ~isempty(row)
            Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(row(1),3))= -1;
            for j = 1:length(row)-1
                Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC, TableC(row(j+1),3)) =   -1;
            end
        end

        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC,TableC(ParentC,5))= -(CVR_P/2)*PLC(TableC(ParentC,2));
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),TableC(ParentC,4))= 1;

        if ~isempty(row)
            Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+ TableC(row(1),3))= -1;
            for j = 1:length(row)-1
                Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(nbusC-1)+TableC(row(j+1),3)) =  -1;
            end
        end

        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1),(TableC(ParentC,5)))= -(CVR_Q/2)*QLC(TableC(ParentC,2));
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,5))= 1;
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),VolttableC(Poc(1)))= -1;
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,3))= 2*(RCC(TC((ParentC),1),TC((ParentC),2)));
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableC(ParentC,4))= 2*(XCC(TC((ParentC),1),TC((ParentC),2)));
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,3))= -RAC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(XAC(TC((ParentC),1),TC((ParentC),2)));
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableA(ParentA,4))= -XAC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(RAC(TC((ParentC),1),TC((ParentC),2)));
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,3))= -RBC(TC((ParentC),1),TC((ParentC),2))- sqrt(3)*(XBC(TC((ParentC),1),TC((ParentC),2)));
        Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+2*(nbusC-1),TableB(ParentB,4))= -XBC(TC((ParentC),1),TC((ParentC),2))+ sqrt(3)*(RBC(TC((ParentC),1),TC((ParentC),2)));
        beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC)= (1-(CVR_P/2))*PLC(TableC(ParentC,2))-PGC(TableC(ParentC,2));
        beq(3*(nbusA-1)+1+3*(nbusB-1)+1+ParentC+(nbusC-1)) = (1-(CVR_Q/2))*QLC(TableC(ParentC,2))-QCC(TableC(ParentC,2));
    end


    Aeq(3*(nbusA-1)+1+3*(nbusB-1)+1+3*(nbusC-1)+1,VolttableC(1)) = 1;
    beq(3*(nbusA-1)+1+3*(nbusB-1)+1+3*(nbusC-1)+1) = VsC;
end

%% DER equation addition
DER_BusA = find(PDERA(:)~=0);
for k22 = 1:size(DER_BusA,1)
    %     (TableA(find(DER_BusA(k22)==TableA(:,2)),4)); % 3 for P control variables
    Aeq((TableA(find(DER_BusA(k22)==TableA(:,2)),4)),end+1) = 1;
end

DER_BusB = find(PDERB(:)~=0);
for k22 = 1:size(DER_BusB,1)
    %     (TableB(find(DER_BusB(k22)==TableB(:,2)),4));
    Aeq((TableB(find(DER_BusB(k22)==TableB(:,2)),4)),end+1) = 1;
end

DER_BusC = find(PDERC(:)~=0);
for k22 = 1:size(DER_BusC,1)
    %     (TableC(find(DER_BusC(k22)==TableC(:,2)),4));
    Aeq((TableC(find(DER_BusC(k22)==TableC(:,2)),4)),end+1) = 1;
end

%% formation of objective function
Tnvar = size(Aeq,2);         %% total number of variables
f = zeros(Tnvar,1);
% f(TableA(1,3)) = 1;
% f(TableB(1,3)) = 1;
% f(TableC(1,3)) = 1;


%% Limit set:
lb1(TableA(1,3):TableA(end,4),1)= (-100*ones((2*size(TableA,1)),1));
lb1(TableA(1,5)-1:TableA(end,5),1)= ((v_min)*ones((size(TableA,1)+1),1));
lb1(TableB(1,3):TableB(end,4),1)= (-100*ones((2*size(TableB,1)),1));
lb1(TableB(1,5)-1:TableB(end,5),1)= ((v_min)*ones((size(TableB,1)+1),1));
lb1(TableC(1,3):TableC(end,4),1)= (-100*ones((2*size(TableC,1)),1));
lb1(TableC(1,5)-1:TableC(end,5),1)= ((v_min)*ones((size(TableC,1)+1),1));

lb =  [lb1;l_bQ];

ub1(TableA(1,3):TableA(end,4),1)= (100*ones((2*size(TableA,1)),1));
ub1(TableA(1,5)-1:TableA(end,5),1)= ((v_max)*ones((size(TableA,1)+1),1));
ub1(TableB(1,3):TableB(end,4),1)= (100*ones((2*size(TableB,1)),1));
ub1(TableB(1,5)-1:TableB(end,5),1)= ((v_max)*ones((size(TableB,1)+1),1));
ub1(TableC(1,3):TableC(end,4),1)= (100*ones((2*size(TableC,1)),1));
ub1(TableC(1,5)-1:TableC(end,5),1)= ((v_max)*ones((size(TableC,1)+1),1));

ub =  [ub1;u_bQ];

%% Nonlinear Solver:
% x0 = 0*ub+1;
% algo1 = 'sqp';
% algo2 = 'active-set';
% tic
% options = optimoptions('fmincon','Display','iter','MaxIterations',2500,'MaxFunctionEvaluations',1000000,'Algorithm',algo2);
% [x,fvalnonlin,exitflag,output] = fmincon(@objfun,x0,[],[],Aeq,beq,lb,ub, [], options);
% time_dist = toc;

%% Linear Solver:

tic;
optimoptions(@intlinprog,'Display','iter');
[x,fval,exitflag,output] = intlinprog(f,[],[],[],Aeq,beq,lb,ub);

mic_iter = fval;
time_dist = toc;

%% Solution of Voltage and S:
for j = 1:size(TableA,1)
    VA(TableA(j,2),1)= x(TableA(j,5));
    S_allA((TableA(j,2))-1,1) = complex(x(TableA(j,3)),x(TableA(j,4)));
end
VA(1) = (VsA);
S1A = S_allA(1);

for j = 1:size(TableB,1)
    VB(TableB(j,2),1)= x(TableB(j,5));
    S_allB((TableB(j,2))-1,1) = complex(x(TableB(j,3)),x(TableB(j,4)));
end
VB(1) = (VsB);
S1B = S_allB(1);

for j = 1:size(TableC,1)
    VC(TableC(j,2),1)= x(TableC(j,5));
    S_allC((TableC(j,2))-1,1) = complex(x(TableC(j,3)),x(TableC(j,4)));
end
VC(1) = (VsC);
S1C = S_allC(1);
V_Linear = sqrt([VA VB VC]);

%% Decission Variable result:
DER_BusABC = [DER_BusA;DER_BusB;DER_BusC];
DER_BusBC = [DER_BusB;DER_BusC];
Dec_Var_Q_onlyA = (x(end-(size(DER_BusABC,1))+1:end-((size(DER_BusBC,1)))));
Dec_Var_Q_onlyB = (x(end-(size(DER_BusBC,1))+1:end-((size(DER_BusC,1)))));
Dec_Var_Q_onlyC = (x(end-(size(DER_BusC,1))+1:end));

Dec_Var = zeros(nb,3);
for j=1:length(DER_BusA)
    Dec_Var(DER_BusA(j),1) = Dec_Var_Q_onlyA(j);
end

for j=1:length(DER_BusB)
    Dec_Var(DER_BusB(j),2) = Dec_Var_Q_onlyB(j);
end
for j=1:length(DER_BusC)
    Dec_Var(DER_BusC(j),3) = Dec_Var_Q_onlyC(j);
end
Dec_VarA=Dec_Var(:,1);
Dec_VarB=Dec_Var(:,2);
Dec_VarC=Dec_Var(:,3);

Dec_var_file = strcat(save_variable_directory,"\Dec_Var");
V_Linear_file = strcat(save_variable_directory,"\V_Linear");


save(Dec_var_file, "Dec_Var")
save(V_Linear_file, "V_Linear")

disp('C-OPF substation Power flow (kW)-')
real(S_allA(1)+S_allB(1)+S_allC(1))*1000









