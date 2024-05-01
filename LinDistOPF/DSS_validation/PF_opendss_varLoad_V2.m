clear all
% clc


%%
%Input:
load Dec_Var.mat
load V_Linear.mat

load_mult = 1;
gen_mult = 1;

% import Dec_Var.
Area = 123;
Q_ctrl  = 1;  % 1=Q ontrol, 0: P control
root_directory = "C:\Users\sadn725\OneDrive - PNNL\PNNL\AppDeconfliction\lindistflow_matlab_code\Test_system";
data_location = strcat(root_directory,"\Area_",num2str(Area));


%% Create Line Code and Line dss files:

linedss(data_location)

%%
% CVR = [0.6 3.0];
% CVR = [1.6 6];
CVR = [0 0];

DSSObj=actxserver('OpenDSSEngine.DSS');
if ~DSSObj.Start(0)
    disp('Unable to start openDSS Engine');
    return
end


DSSText=DSSObj.Text;
DSSCircuit=DSSObj.ActiveCircuit;
DSSBus = DSSCircuit.ActiveBus;

masterfile = strcat(root_directory,"\DSS_validation\itest_system_varL.dss");
masterfile_compile_command = strcat(" Compile (",masterfile,")");
% DSSText.Command='Compile (C:\Users\sadn725\OneDrive - PNNL\PNNL\AppDeconfliction\lindistflow_matlab_code\Test_system\DSS_validation\itest_system_varL.dss)';  % Create a master file for this
DSSText.Command=masterfile_compile_command;  % Create a master file for this

text_powerdata_r = strcat(data_location,'\powerdata.txt');


Power_data = dlmread(text_powerdata_r);
Power_data(:,2:7) = load_mult.*Power_data(:,2:7);
Power_data(:,11:13) = gen_mult.*Power_data(:,11:13);

%% Load dss generation
for i = 1:size(Power_data,1)   
    if (Power_data(i,2)&&Power_data(i,3)~=0)
        % Constant Power load; k = 1-CVR/2:
        sta1 = strcat('New load.S', num2str(Power_data(i,1)),'aP    Bus=',num2str(Power_data(i,1)),...
            '.1   Phases=1   Conn=Wye  Model=1   kv=2.402   kw=',num2str(Power_data(i,2)*(1-(CVR(1))/2)),...
            '  kvar=',num2str(Power_data(i,3)*(1-(CVR(2))/2)),'       Vminpu=0.2  Vmaxpu=1.1'  );
        
        % Constant Z load; k = CVR/2:
        sta2 = strcat('New load.S', num2str(Power_data(i,1)),'aZ    Bus=',num2str(Power_data(i,1)),...
            '.1   Phases=1   Conn=Wye  Model=2   kv=2.402   kw=',num2str(Power_data(i,2)*((CVR(1))/2)),...
            '  kvar=',num2str(Power_data(i,3)*((CVR(2))/2)),'       Vminpu=0.2  Vmaxpu=1.1'  );
        
       DSSText.Command =sta1;
       DSSText.Command =sta2;
    end

    if (Power_data(i,4)&&Power_data(i,5)~=0)
        % Constant Power load; k = 1-CVR/2:
        stb1 = strcat('New load.S', num2str(Power_data(i,1)),'bP    Bus=',num2str(Power_data(i,1)),...
            '.2   Phases=1   Conn=Wye  Model=1   kv=2.402   kw=',num2str(Power_data(i,4)*(1-(CVR(1))/2)),...
            '  kvar=',num2str(Power_data(i,5)*(1-(CVR(2))/2)),'       Vminpu=0.2  Vmaxpu=1.1'  );
        
        % Constant Z load; k = CVR/2:
        stb2 = strcat('New load.S', num2str(Power_data(i,1)),'bZ    Bus=',num2str(Power_data(i,1)),...
            '.2   Phases=1   Conn=Wye  Model=2   kv=2.402   kw=',num2str(Power_data(i,4)*((CVR(1))/2)),...
            '  kvar=',num2str(Power_data(i,5)*((CVR(2))/2)),'       Vminpu=0.2  Vmaxpu=1.1'  );
        
        DSSText.Command =stb1;
        DSSText.Command =stb2;
    end

    if (Power_data(i,6)&&Power_data(i,7)~=0)
        % Constant Power load; k = 1-CVR/2:
        stc1 = strcat('New load.S', num2str(Power_data(i,1)),'cP    Bus=',num2str(Power_data(i,1)),...
            '.3   Phases=1   Conn=Wye  Model=1   kv=2.402   kw=',num2str(Power_data(i,6)*(1-(CVR(1))/2)),...
            '  kvar=',num2str(Power_data(i,7)*(1-(CVR(2))/2)),'       Vminpu=0.2  Vmaxpu=1.1'  );
        
        % Constant Z load; k = CVR/2:
        stc2 = strcat('New load.S', num2str(Power_data(i,1)),'cZ    Bus=',num2str(Power_data(i,1)),...
            '.3   Phases=1   Conn=Wye  Model=2   kv=2.402   kw=',num2str(Power_data(i,6)*((CVR(1))/2)),...
            '  kvar=',num2str(Power_data(i,7)*((CVR(2))/2)),'       Vminpu=0.2  Vmaxpu=1.1'  );
        
        DSSText.Command =stc1;
        DSSText.Command =stc2;
    end
end
%% Generation Data
% text_powerdata_r = 'powerdata.txt';
% Power_data = dlmread(text_powerdata_r);
nBus= size(Power_data,1);

if (Q_ctrl==1)
    % For Q control
    Qgen= zeros(nBus,3);   % All the Q values in all the buses
    Qgen= Dec_Var*1000;   % DOPF
    % Qgen= Dec_var_NLP;   % COPF
    % conn=wye
    % Qgen= Qgen*1000;   % DOPF
    Pgen = Power_data(:,11:13);
else
    % For P control
    Pgen = zeros(nBus,3);
    Pgen= Dec_Var*1000;
    Qgen= 0*Pgen;
    
end
%%
for i = 1:nBus
    if (Power_data(i,11)~=0)
        stGa = strcat('New generator.',num2str(Power_data(i,1)),'a  Phases=1   Bus1=',num2str(Power_data(i,1)),...
            '.1   kV = 2.4 conn=wye Model=1  kW=',num2str(Pgen(i,1)),'   kVAr=',num2str(Qgen(i,1)));
        DSSText.Command =stGa;
    end
    
    if (Power_data(i,12)~=0)
        stGb = strcat('New generator.',num2str(Power_data(i,1)),'b  Phases=1   Bus1=',num2str(Power_data(i,1)),...
            '.2   kV = 2.4 conn=wye  Model=1 kW=',num2str(Pgen(i,2)),'   kVAr=',num2str(Qgen(i,2)));
        DSSText.Command =stGb;
        
    end
    
    if (Power_data(i,13)~=0)
        stGc = strcat('New generator.',num2str(Power_data(i,1)),'c  Phases=1   Bus1=',num2str(Power_data(i,1)),...
            '.3   kV = 2.4 conn=wye Model=1  kW=',num2str(Pgen(i,3)),'   kVAr=',num2str(Qgen(i,3)));
        DSSText.Command =stGc;
    end
end

%% Capacitor -

for i = 1:nBus
        if (Power_data(i,8)~=0)
            stGa = strcat('New generator.C',num2str(Power_data(i,1)),'a  Phases=1   Bus1=',num2str(Power_data(i,1)),...
                '.1   kV = 2.402 Model=1   kW=',num2str(0),'   kVAr=',num2str(Power_data(i,8)));
            DSSText.Command =stGa;
            
        end
        
        if (Power_data(i,9)~=0)
            stGb = strcat('New generator.C',num2str(Power_data(i,1)),'b  Phases=1   Bus1=',num2str(Power_data(i,1)),...
                '.2   kV = 2.402 Model=1  kW=',num2str(0),'   kVAr=',num2str(Power_data(i,9)));
            DSSText.Command =stGb;
           
        end
        
        if (Power_data(i,10)~=0)
            stGc = strcat('New generator.C',num2str(Power_data(i,1)),'c  Phases=1   Bus1=',num2str(Power_data(i,1)),...
                '.3   kV = 2.402 Model=1   kW=',num2str(0),'   kVAr=',num2str(Power_data(i,10)));
            DSSText.Command =stGc;
     
        end
end


%%

DSSText.Command = 'Set VoltageBases = [4.16, ]' ; %  ! ARRAY OF VOLTAGES IN KV
DSSText.Command = 'CalcVoltageBases' ; %! PERFORMS ZERO LOAD POWER FLOW TO ESTIMATE VOLTAGE BASES


DSSText.Command = 'Set maxcontroliter=100';
DSSText.Command = 'solve mode=snap';

DSSSolution     = DSSCircuit.Solution;

DSSObj.AllowForms = false;

DSSSolution.Solve;


V_PU = [];
V_PU = DSSCircuit.AllBusVmagPU;
MytotalCircuitLosses= (DSSCircuit.Losses)/1000;

DSSCircuit.SetActiveElement(['Line.L1']);
MyPowerArray = DSSCircuit.ActiveCktElement.Powers;
P_Sub = sum(MyPowerArray(1:2:6))

V_opds = zeros(nBus,3);
BUS_phase = DSSCircuit.AllNodeNames;
for i = 1:length(BUS_phase)
    X = BUS_phase{i};
    X_len = length(X);
    col = str2num(X(end));
    row = str2num(X(1:(X_len-2)));
    V_opds(row,col) = V_PU(i);
end
%%

% Validation plot:

V_opds_A = V_opds(find(V_opds(:,1)~=0),1);
V_opds_B = V_opds(find(V_opds(:,2)~=0),2);
V_opds_C = V_opds(find(V_opds(:,3)~=0),3);

V_Lin_A = V_Linear(find(V_Linear(:,1)~=0),1);
V_Lin_B = V_Linear(find(V_Linear(:,2)~=0),2);
V_Lin_C = V_Linear(find(V_Linear(:,3)~=0),3);

v_diffA = V_opds_A-V_Lin_A;
v_diffB = V_opds_B-V_Lin_B;
v_diffC = V_opds_C-V_Lin_C;

plot(v_diffA);
hold on
plot(v_diffB);
plot(v_diffC);

