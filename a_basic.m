clear all
clc
close all

% tic
% Definition

% Material
Density_A = 2700;
Modulus_A = 70e9;
Density_PZT = 7500;
% Modulus_PZT = inv(1.65e-011)*(1+0.001*1i); %c11E
Modulus_PZT = inv(1.65e-011); %c11E without damping
e31=-2.74e-010/1.65e-011;
epsilon33=3400*8.8541878176*10^(-12)-(2.74e-010)^2/1.65e-011;
Resistance = logspace(0,6,1000)'; % 저항 Controllable
% Resistance = 1e6; % 저항 Controllable

% Frequency
Freq_t = 10e3;

% Wave
Velocity_A = sqrt(Modulus_A/Density_A);
Wavenumber_t_A = Freq_t*2*pi/Velocity_A;
Wavelength_t_A = 2*pi/Wavenumber_t_A;

% Geometry: PZT/substrate thickness = 0.3까지는 잘 맞음.
Length_PZT = 50e-3;
% Length_PZT = linspace(0.01*Wavelength_t_A,5*Wavelength_t_A,10000);
Height_A = 5e-3;
% Height_PZT = 0.5e-3;
Height_PZT = 0.5e-3;
Width_A = 5e-3;
Area_A = Height_A*Width_A;
Impedance_A = Density_A*Velocity_A*Area_A;

% 
Mass_effective =  Width_A*(Density_A*Height_A+Density_PZT*Height_PZT*2); % per length
Stiffness_effective = Width_A*(Modulus_A*Height_A+Modulus_PZT*Height_PZT*2); % per length
Velocity_effective = sqrt(Stiffness_effective/Mass_effective);
Impedance_effective = Mass_effective*Velocity_effective;
Capacitance_effective = Length_PZT*Width_A*epsilon33/Height_PZT/2; % Cp/2

Wavenumber_effective = Freq_t*2*pi/Velocity_effective;
Wavelength_effective = 2*pi/Wavenumber_effective;

for count = 1:length(Height_PZT)
    for id_res = 1:length(Resistance)
        Kcoeff = (Width_A*e31)^2/sqrt(Mass_effective(count)*Stiffness_effective(count));
        kcoupling(count,id_res) = (1-cos(Wavenumber_effective(count)*Length_PZT))*Kcoeff/(1/Resistance(id_res)+1i*2*pi*Freq_t*Capacitance_effective(count)-(exp(-1i*Length_PZT*Wavenumber_effective(count))-1)*Kcoeff);
        Transfer_effective = 1/(1-kcoupling(count,id_res))*([cos(Wavenumber_effective(count)*Length_PZT) 1i/Impedance_effective(count)*sin(Wavenumber_effective(count)*Length_PZT);
            1i*Impedance_effective(count)*sin(Wavenumber_effective(count)*Length_PZT) cos(Wavenumber_effective(count)*Length_PZT)]+...
            [-kcoupling(count,id_res)*exp(-1i*Wavenumber_effective(count)*Length_PZT) kcoupling(count,id_res)*(1+exp(-1i*Wavenumber_effective(count)*Length_PZT))/Impedance_effective(count)
            Impedance_effective(count)*kcoupling(count,id_res)*(exp(-1i*Wavenumber_effective(count)*Length_PZT)-1) -kcoupling(count,id_res)*exp(-1i*Wavenumber_effective(count)*Length_PZT)]);
%         Teff_11(count,id_res) = Transfer_effective(1,1);
%         k_effective(count,id_res) = 1./Length_PZT.*acos(Teff_11(count,id_res));
%         lambda_effective(count,id_res) = 2*pi./k_effective(count,id_res);
        
        Impedance_Matrix_A = 1i*2*pi*Freq_t*[1 1;-Impedance_A  Impedance_A];
        Impedance_Matrix_effective = 1i*2*pi*Freq_t*[1+kcoupling(count,id_res)*exp(-1i*Wavenumber_effective(count)*Length_PZT) 1-kcoupling(count,id_res);
            Impedance_effective(count)*(-1+kcoupling(count,id_res)*exp(-1i*Wavenumber_effective(count)*Length_PZT))  Impedance_effective(count)*(1-kcoupling(count,id_res))];
        
        S_paramter_PEH = inv(Impedance_Matrix_A)*Transfer_effective*Impedance_Matrix_A;
        Reflection(count,id_res) = -(S_paramter_PEH(2,1)/S_paramter_PEH(2,2));
        Transmission(count,id_res) = (1/S_paramter_PEH(2,2));
%         Localization(:,count,id_res) = inv(Impedance_Matrix_effective)*Impedance_Matrix_A*[1; Reflection(count,id_res)]*10*10^(-9);
        Aeff = 1/2/Impedance_effective(count)*(Impedance_effective(count)+Impedance_A+Reflection(count,id_res)*(Impedance_effective(count)-Impedance_A))*10*10^(-9);
        Beff = 1/2/Impedance_effective(count)/(1-kcoupling(count,id_res))*(Impedance_effective(count)-Impedance_A+Reflection(count,id_res)*((Impedance_effective(count)+Impedance_A)-kcoupling(count,id_res)*(Impedance_effective(count)+Impedance_A+Reflection(count,id_res)*(Impedance_effective(count)-Impedance_A))*exp(-i*Wavenumber_effective(count)*Length_PZT)))*10*10^(-9);
        Voltage(count,id_res) = (1i*Width_A*e31*2*pi*Freq_t)/(1/Resistance(id_res)+1i*2*pi*Freq_t*Capacitance_effective(count)-(Kcoeff*(exp(-1i*Length_PZT*Wavenumber_effective(count))-1)))...
            *(Aeff*(exp(-1i*Wavenumber_effective(count)*Length_PZT)-1)+Beff*(exp(1i*Wavenumber_effective(count)*Length_PZT)-1));
        
        disp_PZT(count,id_res) = Aeff + Beff + (Width_A*e31)^2/sqrt(Mass_effective(count)*Stiffness_effective(count))/2/(1/Resistance(id_res)+1i*2*pi*Freq_t*Capacitance_effective(count)-(Kcoeff*(exp(-1i*Length_PZT*Wavenumber_effective(count))-1)))...
            *(1-exp(-1i*Wavenumber_effective(count)*Length_PZT))*(Aeff*(exp(-1i*Wavenumber_effective(count)*Length_PZT)-1)+Beff*(exp(1i*Wavenumber_effective(count)*Length_PZT)-1));
        %         Voltage(count,id_res) = (1i*Width_A*e31*2*pi*Freq_t)/(1/Resistance(id_res)+1i*2*pi*Freq_t*Capacitance_effective(count)-(Kcoeff*(exp(-1i*Length_PZT(count)*Wavenumber_effective(count))-1)))...
        %             *(Localization(1,count,id_res)*(exp(-1i*Wavenumber_effective(count)*Length_PZT(count))-1)+Localization(2,count,id_res)*(exp(1i*Wavenumber_effective(count)*Length_PZT(count))-1));
        Power(count,id_res) = Voltage(count,id_res)^2/Resistance(id_res);
        Current(count,id_res) = Voltage(count,id_res)/Resistance(id_res);
%         energy(count,id_res) = abs(Reflection(count,id_res))^2 + abs(Transmission(count,id_res))^2;
    end
end
%%
width = 8/2.54*300;
height = 6/2.54*300;
fontsize = 20;
linewidth = 2;
marker_size = 10;
%% Voltage response across resistances
figure;
plot(Resistance, abs(Voltage),'-','linewidth',2); grid on; hold on;
set(gca, 'xscale','log')
set(gca, 'fontsize',fontsize)
ylabel('Voltage')

data0 = load('voltage_resistance.txt');
res_FEM = data0(:,1);
vol_FEM = data0(:,2);
plot(res_FEM, vol_FEM,'o','MarkerSize',marker_size);

legend('Analytical model','COMSOL')


