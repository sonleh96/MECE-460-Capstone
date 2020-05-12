% This script calculates the energy requirements for the
% overwintering building in different locations with hourly data

function PowerAnalysis_SP(filename)

%Input Variable definitions
Q_total = 10.694; %Total airflow [m3/s]
angle = 21; %angle of tilt of solar heater, [deg]
p_supply = 0.7; %Proportion of supply air that is fresh
Tr_glass = 0.83; %Transmittance of glass
L_wall = 42.7; %Length of wall (m)
H_wall = 4.57; %Height of wall (m)
C_elec = 6.32; %Cost of electricity, cents per kWh
C_NG = 1.91; %Cost of natural gas, dollars per GJ
T_desired = 4; %Desired entering temperature
C_min = 909; %From HEx calcs [W/K]
T_w_in = 65; %Water temperature into HEx [C]
Eff_boiler = 0.9; %Efficiency of boiler
CO2_limit = 3500; %Limit at which system goes into CO2 Exhaust mode
CO2_ext = 408; %Outside CO2 concentration [ppm]
P_f_ret = 0.7755; %Operating power of return fan [kW]
P_f_sup = 11.558; %Operating power of supply fan [kW]
P_f_exh = 0.2237; %Operating power of a single exhaust fan [kW]
P_f_pump = 0.1243; %Operating power of the pump [kW]
CO2_bees = 30000; %Additional CO2 added to system from bees [g/hr]
Vol_bldg = 6905; %Volume of building [m3]
ACH = 8; %Number of air changes per hour
q_bees = 58000; %Energy given off by bees [W]
H2O_bees = 30000; %g/hr of water vapour 
HumLim = 85; %Humidity limit

%Air properties (taken at average winter temperature of -10C)
dvis_air = 16.65e-6; %Dynamic viscosity [Ns/m2]
dens_air = 1.33; %Density of air [kg/m3]
cp_air = 1004; %heat capacity air [J/kgK]
k_air = 0.02375; %Thermal conductivity of air [W/mK]
Pr_air = (dvis_air * cp_air) / k_air; %Prandtl number

%Calculated variables (constants not requiring iteration)
CO2_bees_ppm = ((CO2_bees / ACH / 44.01) / ((dens_air / 0.02897)*...
            Vol_bldg)) * 1000000; %ppm of CO2 generated per ACH

%Finding glass surface area (Solar Heater)
L_glass = H_wall / cosd(angle); %Length of glass (hypotenuse)
A_glass = L_glass * L_wall; %Area of heat transfer
W_SAH = H_wall * tand(angle); %Width of base of triangle

%Transferring data from excel into matrices
Temp = xlsread(filename,'Temp');
Hum = xlsread(filename,'Hum');
WSpd = xlsread(filename,'WSpd');
SolB = xlsread(filename,'SolB');
SolD = xlsread(filename,'SolD');
n = length(Temp);

%Checking to ensure length of matrices is the same
if isequal(length(Temp),length(Hum),length(WSpd),length(SolB)) == 0
    fprintf('Datasets not equal. Check excel lengths.');
    return
end

%Solar Air Heater - Internal Convection coefficient (doesn't change with
%hourly climate data)
Perim = H_wall + W_SAH + L_glass; %Perimeter of triangular x-section [m]
A_SAH = H_wall * W_SAH * 0.5; %X-sectional area of SAH [m2]
Q_SAH = Q_total * p_supply; %Airflow thru SAH [m3/s]
v_SAH = Q_SAH / A_SAH; %velocity of air thru SAH [m/s]
D_h = (4 * A_SAH) / Perim; %Hydraulic diameter of SAH [m]
Re_SAH = (dens_air * v_SAH * D_h) / dvis_air; %Reynolds # of SAH airflow
Nu_SAH = 0.023 * (Re_SAH^0.8) * (Pr_air^0.4); %Nusselt # of SAH airflow
h_int = (Nu_SAH * k_air) / D_h; %Internal SAH convection coefficient
                                %[W/m2K]

%Initial weight of water in building
AHum = zeros(1,n-1);
AHum(1) = (6.112*exp((17.67*Temp(1))/(Temp(1)+243.5))*Hum(1)*2.1674) / ...
    (Temp(1)+273.15);
AHum_old = AHum(1);
                                
%Main FOR loop, iterating per hour of winter (Oct - Mar)
%Initial values for tracked info:
CO2_old = CO2_ext;
Ins = zeros(1,n-1);
P_elec = zeros(1,n-1);
P_gas = zeros(1,n-1);
CO2 = zeros(1,n-1);
T = zeros(1,n-1);
T_ent = zeros(1,n-1);
CostE = zeros(1,n-1);
CostNG = zeros(1,n-1);
HumInt = zeros(1,n-1);
T_ent2 = zeros(1,n-1);
P_gas2 = zeros(1,n-1);
CostNG2 = zeros(1,n-1);
for i = 1:(n-1)
    
    %Gathering solar data via NREL's program
    i_angle = spa_SP(i);
    
    %Picking data out of initial matrices
    v_ext = WSpd(i) / 3.6; %Wind speed, converting to [m/s]
    T_ext = Temp(i);
    Hum_ext = Hum(i);
    Ins(i) = (SolB(i)*sind(i_angle)) + (SolD(i)); %Solar insolation
    %is the sum of direct and diffuse irradiation [kJ/m2/hr]
    
    %Finding convection coefficient of outside of SAH (glass)
        Re_ext = (dens_air * v_ext * L_wall) / dvis_air;
        Nu_ext = 0.037 * (Re_ext^0.8) * (Pr_air^(1/3));
        h_ext = (Nu_ext * k_air) / L_wall;

        %Total thermal resistance
        R = (1/(h_int*A_glass)) + (1/(h_ext*A_glass)); %[K/W]

        %Finding Solar energy transmitted into SAH
        q_solar = (Ins(i)*(1000/3600)) * A_glass * Tr_glass; %[W]

        %Finding heat loss out of SAH due to convection on glass
        q_SAH_tot = q_solar;
        %Iterative process, from initial calculations it takes 5
        %iterations to converge
        for j = 1:5
            T_SAH_f = (q_SAH_tot / (Q_total * dens_air * cp_air)) + ...
                T_ext; %Calculating final temperature of SAH
            T_avg = (T_SAH_f + T_ext) / 2; %[C]
            q_SAH_tot = q_solar - abs((T_avg - T_ext)/R); %Updating total q
        end
        
        T_ent(i) = (T_SAH_f*p_supply) + (T_desired*(1-p_supply));
        T_ent2(i) = (T_ext*p_supply) + (T_desired*(1-p_supply));
        
        %Entering absolute humidity                                 
        AHum_ent = (6.112*exp((17.67*T_ext)/(T_ext+243.5))*Hum_ext*...
            2.1674) / (T_ext+273.15); 
   
        
    %If statement depending on temperature
    if (T_ent(i) < T_desired)
            
        
        %Power required from HEx
        Eff_HEx = (T_desired - T_ent(i)) / (T_w_in - T_ent(i)); %HEx Eff.
        q_HEx = Eff_HEx * C_min * (T_w_in - T_ent(i)); %Heat required for H.
                                                    %coil [W]
        P_boiler = q_HEx / Eff_boiler; %Fuel usage in boiler [W]
        
        %Power if SAH not used
        Eff_HEx2 = (T_desired - T_ent2(i)) / (T_w_in - T_ent2(i)); %HEx Eff.
        q_HEx2 = Eff_HEx2 * C_min * (T_w_in - T_ent2(i)); %Heat required for H.
                                                    %coil [W]
        P_boiler2 = q_HEx2 / Eff_boiler; %Fuel usage in boiler [W]
        
        %Fan Power
        P_fans = P_f_ret + P_f_sup + (4*P_f_exh);
        
        %Calculating the increase in the CO2 concentration
        CO2(i) = CO2_old;
        for k = 1:ACH %air changes per hour
        CO2(i) = CO2(i) + (p_supply*CO2_ext) - (p_supply*CO2(i)) + ...
            CO2_bees_ppm;
        end
        
        CO2_old = CO2(i);
        T(i) = T_desired;
        
        P_elec(i) = P_fans + P_f_pump; %Fans + pump
        P_gas(i) = P_boiler * 0.0000036; %Boiler is the source of natural
                                         %gas usage, convert to GJ
        P_gas2(i) = P_boiler2 * 0.0000036;
        
                                        
        %Humidity Changes:
        for m = 1:ACH
            AHum(i) = (AHum_ent*p_supply) + AHum_old + (H2O_bees/(ACH*Vol_bldg)) - ...
                (p_supply*AHum_old);
            AHum_old = AHum(i);
        end
        
        %Relative humidity
        HumInt(i) = abs((AHum(i)*(273.15*T(i)))/(6.112*exp((17.67*T(i))/...
            (T(i)+243.4))*10.1674));
                                         
    else
        T_ent(i) = T_ext; %If entering temperature too high, open up solar
                          %heater bypass.
    
        P_gas(i) = 0;
        
        CO2(i) = CO2_old;
        T(i) = T_desired;
        P_fans = 0;
        for k = 1:ACH %air changes per hour
            if or((CO2(i) >= CO2_limit),(AHum_old > HumLim))
                CO2(i) = CO2(i) - (p_supply*CO2(i)) + CO2_bees_ppm;
                T(i) = (T(i)*(1-p_supply)) + (T_ent(i)*p_supply);
                P_fans = P_fans + ((4*P_f_exh)/ACH); %No pump
                
                AHum(i) = AHum_old - (p_supply*AHum_old) + ...
                    (H2O_bees/(ACH*Vol_bldg));
                AHum_old = AHum(i);

            else
                CO2(i) = CO2(i) + CO2_bees_ppm;
                T(i) = T(i) + (q_bees/(Q_total*(1-p_supply)*dens_air*...
                    cp_air*(3600/ACH)));
                P_fans = P_fans + (P_f_ret/ACH); %No pump
               
                %Relative humidity
                AHum(i) = AHum_old + (H2O_bees/(ACH*Vol_bldg));
                AHum_old = AHum(i);
                
            end
        end
        
        P_elec(i) = P_fans;
        CO2_old = CO2(i);
        
        HumInt(i) = abs((AHum(i)*(273.15*T(i)))/(6.112*...
                    exp((17.67*T(i))/(T(i)+243.4))*21.674));
        
    end
    
   %Cost calcs
   CostE(i) = (C_elec * P_elec(i)) / 100;
   CostNG(i) = C_NG * P_gas(i);
   CostNG2(i) = C_NG * P_gas2(i);
   
   fprintf('\n %i %.2f %.2f %.2f %.2f %.2f',i,T_ext,T_ent(i),CO2(i),sum(CostE),sum(CostNG)); 
    
end

fprintf('\n');

%Combining data into 1 large array
Data = vertcat(P_gas,CostNG,P_gas2,CostNG2);
Data = transpose(Data);

%Writing data to an excel file
excelfile = 'SP_SAHtest.xlsx';
xlswrite(excelfile,Data);

end
            
        
    




    


