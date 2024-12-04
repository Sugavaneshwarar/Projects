clear all
clc

% Create Figure

DHTE = figure('Position', [100, 100, 1400, 600], 'Name', 'Design of Shell and Tube Heat Exchanger');
fontsize = 10;
headingsize = 11;

% Shell Fluid Temperatures

uicontrol('Style', 'text', 'Position', [50, 550, 200, 40], 'String', 'Shell Fluid Temperatures', 'FontSize', headingsize);
uicontrol('Style', 'text', 'Position', [0, 525, 200, 40], 'String', 'Inlet Temperature(K)', 'FontSize', fontsize);
uicontrol('Style', 'text', 'Position', [0, 500, 205, 40], 'String', 'Outlet Temperature(K)', 'FontSize', fontsize);
shell_inlet_temp = uicontrol('Style', 'edit', 'Position', [190, 545, 80, 18], 'String', '0', 'FontSize', fontsize, 'HorizontalAlignment', 'center');
shell_outlet_temp = uicontrol('Style', 'edit', 'Position', [190, 520, 80, 18], 'String', '0', 'FontSize', fontsize, 'HorizontalAlignment', 'center');

% Tube Fluid Temperatures

uicontrol('Style', 'text', 'Position', [290, 550, 200, 40], 'String', 'Tube Fluid Temperatures', 'FontSize', headingsize);
uicontrol('Style', 'text', 'Position', [280, 525, 150, 40], 'String', 'Inlet Temperature(K)', 'FontSize', fontsize);
uicontrol('Style', 'text', 'Position', [280, 495, 150, 40], 'String', 'Outlet Temperature(K)', 'FontSize', fontsize);
tube_inlet_temp = uicontrol('Style', 'edit', 'Position', [427, 545, 80, 18], 'String', '0', 'FontSize', fontsize, 'HorizontalAlignment', 'center');
tube_outlet_temp = uicontrol('Style', 'edit', 'Position', [427, 518, 80, 18], 'String', '0', 'FontSize', fontsize, 'HorizontalAlignment', 'center');

% Type of Flow 

uicontrol('Style', 'text', 'Position', [510, 550, 100, 40], 'String', 'Flow Type', 'FontSize', headingsize);
flow_type = uicontrol('Style', 'popupmenu', 'Position', [530, 515, 70, 40], 'String', {'Co-flow', 'Counter-Flow'}, 'FontSize', fontsize);

% Mass Flow Rates

uicontrol('Style', 'text', 'Position', [640, 550, 120, 40], 'String', 'Mass Flow Rates', 'FontSize', headingsize);
uicontrol('Style', 'text', 'Position', [605, 520, 120, 40], 'String', 'Shell Side (kg/s)', 'FontSize', fontsize);
uicontrol('Style', 'text', 'Position', [605, 495, 120, 40], 'String', 'Tube Side (kg/s)', 'FontSize', fontsize);
shell_rate = uicontrol('Style', 'edit', 'Position', [720, 540, 80, 18], 'String', '0', 'FontSize', fontsize, 'HorizontalAlignment', 'center');
tube_rate = uicontrol('Style', 'edit', 'Position', [720, 515, 80, 18], 'String', '0', 'FontSize', fontsize, 'HorizontalAlignment', 'center');

% Exchanger Dimensions 

uicontrol('Style','text','Position', [30,440,200,40],'String','Heat Exchanger Specifications','FontSize',headingsize);
uicontrol('Style', 'text', 'Position', [0, 400, 200, 40], 'String', 'Length(m)', 'FontSize', fontsize);
heat_exchanger_length = uicontrol('Style', 'edit', 'Position', [190, 420, 80, 18], 'String', '0', 'FontSize', fontsize);
uicontrol('Style', 'text', 'Position', [0, 360, 200, 40], 'String', 'Tube Inner Diameter(m)', 'FontSize', fontsize);
uicontrol('Style', 'text', 'Position', [10, 320, 200, 40], 'String', 'Tube Thickness(m)', 'FontSize', fontsize);
tube_inner_diameter = uicontrol('Style', 'edit', 'Position', [190, 380, 80, 18], 'String', '0', 'FontSize', fontsize, 'HorizontalAlignment', 'center');
tube_thickness = uicontrol('Style', 'edit', 'Position', [190, 340, 80, 18], 'String', '0', 'FontSize', fontsize, 'HorizontalAlignment', 'center');
uicontrol('Style', 'text', 'Position', [0, 280, 200, 40], 'String', 'Shell Inner Diameter(m)', 'FontSize', fontsize);
shell_inner_diameter = uicontrol('Style', 'edit', 'Position', [190, 300, 80, 18], 'String', '0', 'FontSize', fontsize, 'HorizontalAlignment', 'center');
uicontrol('Style', 'text', 'Position', [10, 240, 200, 40], 'String', 'Tube Thermal Conductivity(W/m.K)', 'FontSize', fontsize);
tube_material_k = uicontrol('Style', 'edit', 'Position', [190, 255, 80, 18], 'String', '0', 'FontSize', fontsize, 'HorizontalAlignment', 'center');

% Tube Fluid Properties

uicontrol('Style', 'text', 'Position', [320, 460, 140, 40], 'String', 'Tube Fluid Properties', 'FontSize', headingsize);
uicontrol('Style', 'text', 'Position', [290, 420, 100, 40], 'String', 'Density(kg/m^3)', 'FontSize', fontsize);
tube_fluid_density = uicontrol('Style', 'edit', 'Position', [405, 440, 80, 18], 'String', '0', 'FontSize', fontsize, 'HorizontalAlignment', 'center');
uicontrol('Style', 'text', 'Position', [290, 375, 100, 40], 'String', 'Viscosity(Pa.s)', 'FontSize', fontsize);
tube_fluid_viscosity = uicontrol('Style', 'edit', 'Position', [405, 395, 80, 18], 'String', '0', 'FontSize', fontsize, 'HorizontalAlignment', 'center');
uicontrol('Style', 'text', 'Position', [278, 335, 120, 40], 'String', 'Specific Heat Capacity(J/kg.K)', 'FontSize', fontsize);
tube_fluid_Cp = uicontrol('Style', 'edit', 'Position', [405, 350, 80, 18], 'String', '0', 'FontSize', fontsize, 'HorizontalAlignment', 'center');
uicontrol('Style', 'text', 'Position', [270, 290, 140, 40], 'String', 'Thermal Conductivity(W/m.K)', 'FontSize', fontsize);
tube_fluid_k = uicontrol('Style', 'edit', 'Position', [405, 305, 80, 18], 'String', '0', 'FontSize', fontsize, 'HorizontalAlignment', 'center');

% Shell Fluid Properties

uicontrol('Style', 'text', 'Position', [320, 235, 140, 40], 'String', 'Shell Fluid Properties', 'FontSize', headingsize);
uicontrol('Style', 'text', 'Position', [290, 200, 100, 40], 'String', 'Density(kg/m^3)', 'FontSize', fontsize);
shell_fluid_density = uicontrol('Style', 'edit', 'Position', [405, 220, 80, 18], 'String', '0', 'FontSize', fontsize, 'HorizontalAlignment', 'center');
uicontrol('Style', 'text', 'Position', [290, 160, 100, 40], 'String', 'Viscosity(Pa.s)', 'FontSize', fontsize);
shell_fluid_viscosity = uicontrol('Style', 'edit', 'Position', [405, 180, 80, 18], 'String', '0', 'FontSize', fontsize, 'HorizontalAlignment', 'center');
uicontrol('Style', 'text', 'Position', [278, 120, 120, 40], 'String', 'Specific Heat Capacity(J/kg.K)', 'FontSize', fontsize);
shell_fluid_Cp = uicontrol('Style', 'edit', 'Position', [405, 135, 80, 18], 'String', '0', 'FontSize', fontsize, 'HorizontalAlignment', 'center');
uicontrol('Style', 'text', 'Position', [263, 70, 140, 40], 'String', 'Thermal Conductivity(W/m.K)', 'FontSize', fontsize);
shell_fluid_k = uicontrol('Style', 'edit', 'Position', [405, 85, 80, 18], 'String', '0', 'FontSize', fontsize, 'HorizontalAlignment', 'center');

% Fluids

uicontrol('Style', 'text', 'Position', [860, 550, 100, 40], 'String', 'Fluids', 'FontSize', headingsize);
uicontrol('Style', 'text', 'Position', [810, 520, 80, 40], 'String', 'Tube Fluid', 'FontSize', fontsize);
tube_fluid = uicontrol('Style', 'popupmenu', 'Position', [900, 522, 70, 40], 'String', {'Other', 'Water', 'Ethylene Glycol'}, 'FontSize', fontsize, 'Callback', @(src, event) TubeValue(src, tube_fluid_density, tube_fluid_viscosity, tube_fluid_Cp, tube_fluid_k));
uicontrol('Style', 'text', 'Position', [810, 485, 80, 40], 'String', 'Shell Fluid', 'FontSize', fontsize);
shell_fluid = uicontrol('Style', 'popupmenu', 'Position', [900, 490, 70, 40], 'String', {'Other', 'Water', 'Ethylene Glycol'}, 'FontSize', fontsize,'Callback', @(src, event) ShellValue(src, shell_fluid_density, shell_fluid_viscosity, shell_fluid_Cp, shell_fluid_k));

% Fouling Factors 

uicontrol('Style','text','Position',[500,460,200,40],'String','Fouling Factors','FontSize',headingsize);
uicontrol('Style','text','Position',[500,415,120,40],'String','Internal Fouling Factor(m^2.K/W)','FontSize',fontsize);
fouling_factor_inside = uicontrol('Style','edit','Position',[620,430,80,18],'String','0','FontSize',fontsize,'HorizontalAlignment','center');
uicontrol('Style','text','Position',[500,360,120,40],'String','External Fouling Factor(m^2.K/W)','FontSize',fontsize);
fouling_factor_outside = uicontrol('Style','edit','Position',[620,370,80,18],'String','0','FontSize',fontsize,'HorizontalAlignment','center');

% Alternate Option : Directly Input Value of U

uicontrol('Style','text','Position',[100,200,30,40],'String','(OR)','FontSize',fontsize);
uicontrol('Style','text','Position',[50,170,130,50],'String','Overall Heat Transfer Coefficient (W/m^2.K)','FontSize',fontsize,'HorizontalAlignment','center');
overall_heat_transfer_coefficient = uicontrol('Style','edit','Position',[190,200,80,18],'String','0','FontSize',fontsize,'HorizontalAlignment','center');

% Assumptions

uicontrol('Style', 'text', 'Position', [0, 27, 350, 40], 'String', 'Assumptions: 1) Fluid Properties do not change with temperature.', 'FontSize', 8);
uicontrol('Style', 'text', 'Position', [85, 13, 240, 40], 'String', '2) Both fluids enter the exchanger at exactly t = 0', 'FontSize', 8);
uicontrol('Style', 'text', 'Position', [85, 0, 290, 40], 'String', '3) Flow regimes other than laminar and turbulent are ignored.');

% Formula Display 

uicontrol('Style','text','Position',[530,300,120,40],'String','Formulas Used','FontSize',headingsize);
uicontrol('Style','text','Position',[530,260,120,40],'String','Q = U*A*LMTD','FontSize',8);
uicontrol('Style','text','Position',[470,240,260,30],'String','1/U = 1/hi + Rfi + ln(do/di)/2*pi*k*L + Rfo + 1/ho','FontSize',8);


% Axes for temperature profile
axes_position = [0.55, 0.1, 0.4, 0.68];
axes_handle = axes('Parent', DHTE, 'Position', axes_position);
xlabel('Length (m)');
ylabel('Temperature (K)');
title('Temperature Profile');
hold(axes_handle, 'on');
uicontrol('Style','text','Position',[1130,469,200,40],'String','The tip of the triangles point to the direction of flow','FontSize',9);

% Calculation Button

uicontrol('Style', 'pushbutton', 'Position', [50,90,180,40], 'String', 'CALCULATE', 'FontSize',16, 'Callback', {@CalculateLMTD, shell_inlet_temp, tube_inlet_temp, shell_outlet_temp, tube_outlet_temp, ...
    tube_inner_diameter, tube_thickness, shell_inner_diameter, heat_exchanger_length, tube_rate, shell_rate, tube_fluid_density, shell_fluid_density, tube_fluid_viscosity, ...
    shell_fluid_viscosity, tube_fluid_Cp, shell_fluid_Cp, tube_fluid_k, shell_fluid_k, tube_material_k, flow_type, overall_heat_transfer_coefficient, axes_handle,fouling_factor_inside,fouling_factor_outside});

% Function to automatically input shell fluid properties 

function ShellValue(src, shell_fluid_density, shell_fluid_viscosity, shell_fluid_Cp, shell_fluid_k) 
Val = src.Value;
if Val == 2
    shell_fluid_density.String = '1000';
        shell_fluid_viscosity.String = '0.001';
        shell_fluid_Cp.String = '4187';
        shell_fluid_k.String = '0.598';
elseif Val == 3
    shell_fluid_density.String = '1113';
        shell_fluid_viscosity.String = '0.018376';
        shell_fluid_Cp.String = '1000.4';
        shell_fluid_k.String = '0.254';
else
    shell_fluid_density.String = '0';
        shell_fluid_viscosity.String = '0';
        shell_fluid_Cp.String = '0';
        shell_fluid_k.String = '0';
end
end

% Function to automatically input fluid properties of tube fluid properties of tube fluid

function TubeValue(src, tube_fluid_density, tube_fluid_viscosity, tube_fluid_Cp, tube_fluid_k)
    val = src.Value;
    if val == 2
        tube_fluid_density.String = '1000';
        tube_fluid_viscosity.String = '0.001';
        tube_fluid_Cp.String = '4187';
        tube_fluid_k.String = '0.598';
    elseif val == 3
        tube_fluid_density.String = '1113';
        tube_fluid_viscosity.String = '0.018376';
        tube_fluid_Cp.String = '1000.4';
        tube_fluid_k.String = '0.254';
    else
        tube_fluid_density.String = '0';
        tube_fluid_viscosity.String = '0';
        tube_fluid_Cp.String = '0';
        tube_fluid_k.String = '0';
    end
end

% Function to Calculate Area of Heat Exchanger through LMTD 

function CalculateLMTD(~, ~, shell_inlet_temp, tube_inlet_temp, shell_outlet_temp, tube_outlet_temp, tube_inner_diameter, tube_thickness, shell_inner_diameter, ...
    heat_exchanger_length, tube_rate, shell_rate, tube_fluid_density, shell_fluid_density, tube_fluid_viscosity, shell_fluid_viscosity, tube_fluid_Cp, shell_fluid_Cp, tube_fluid_k, ...
    shell_fluid_k, tube_material_k, flow_type, overall_heat_transfer_coefficient, axes_handle,fouling_factor_inside,fouling_factor_outside)

% Calculation of Log Mean Temperature Difference

    T1 = str2double(shell_inlet_temp.String);
    t1 = str2double(tube_inlet_temp.String);
    T2 = str2double(shell_outlet_temp.String);
    t2 = str2double(tube_outlet_temp.String);
    di = str2double(tube_inner_diameter.String);
    t = str2double(tube_thickness.String);
    Di = str2double(shell_inner_diameter.String);
    L = str2double(heat_exchanger_length.String);
    m_tube = str2double(tube_rate.String);
    m_shell = str2double(shell_rate.String);
    rho_tube_fluid = str2double(tube_fluid_density.String);
    rho_shell_fluid = str2double(shell_fluid_density.String);
    tube_Cp = str2double(tube_fluid_Cp.String);
    tube_k = str2double(tube_fluid_k.String);
    shell_Cp = str2double(shell_fluid_Cp.String);
    shell_k = str2double(shell_fluid_k.String);
    k = str2double(tube_material_k.String);
    u = str2double(overall_heat_transfer_coefficient.String);
    mu_tube_fluid = str2double(tube_fluid_viscosity.String);
    mu_shell_fluid = str2double(shell_fluid_viscosity.String);
    Rfi = str2double(fouling_factor_inside);
    Rfo = str2double(fouling_factor_outside);

    if T1 > t1
        Th_in = T1;
        Th_out = T2;
        Tc_in = t1;
        Tc_out = t2;
    else
        Th_in = t1;
        Th_out = t2;
        Tc_in = T1;
        Tc_out = T2;
    end

     ftype = flow_type.Value;

    switch ftype
        case 1
            dT1 = Th_in - Tc_in;
            dT2 = Th_out - Tc_out;
        case 2
            dT1 = Th_in - Tc_out;
            dT2 = Th_out - Tc_in;
    end

    LMTD = (dT1 - dT2) / log(dT1 / dT2);

    if(u==0)

    % Calculation of U (Overall Heat Transfer Coefficient)

    A_tube = (pi / 4) * (di)^2;
    do = di + 2 * t;
    d_annulus = Di - do;
    A_annulus = (pi / 4) * (Di^2 - do^2);

    v_tube_fluid = m_tube / (rho_tube_fluid * A_tube);
    v_shell_fluid = m_shell / (rho_shell_fluid * A_annulus);

    NRe_tube_fluid = (rho_tube_fluid * v_tube_fluid * di) / mu_tube_fluid;
    NRe_shell_fluid = (rho_shell_fluid * v_shell_fluid * d_annulus) / mu_shell_fluid;

    Pr_tube_fluid = (mu_tube_fluid * tube_Cp) / tube_k;
    Pr_shell_fluid = (rho_shell_fluid * shell_Cp) / shell_k;

    if NRe_tube_fluid < 1000
        Nu_tube_fluid = 0.023 * NRe_tube_fluid^0.8 * Pr_tube_fluid^0.4;
    else
        Nu_tube_fluid = 0.642 * NRe_tube_fluid^0.5 * Pr_tube_fluid^(1 / 3); 
    end

    h_tube_side = (Nu_tube_fluid * tube_k) / di;

    if NRe_shell_fluid < 1000
        Nu_shell_fluid = 0.023 * NRe_shell_fluid^0.8 * Pr_shell_fluid^0.4;
    else
        Nu_shell_fluid = 0.642 * NRe_shell_fluid^0.5 * Pr_shell_fluid^(1 / 3);
    end

    h_shell_side = (Nu_shell_fluid * shell_k) / d_annulus;
   
    U = 1 / ((1 / h_tube_side) + Rfi + (1 / h_shell_side) + Rfo + (log(do / di) / (2 * pi * k * L)));

    else

         % Directly Taking Value of U from Alternative Input

        U = u;
    end
    disp(U);
    clc
    % Total Heat Transfer (Q)

    if(T1>t1)
        Q = m_shell*shell_Cp*(Th_in-Th_out);
    else
        Q = m_tube*tube_Cp*(Th_in-Th_out);
    end

    A_heat_exchanger = Q/(U*LMTD);

    x = linspace(0,L,100);

        if ftype == 1
            Th = linspace(Th_in,Th_out,100);
            Tc = linspace(Tc_in,Tc_out,100);
        else
            Th = linspace(Th_in,Th_out,100);
            Tc = linspace(Tc_out,Tc_in,100);
        end

         % Clear previous plot
        cla(axes_handle);
        
        % Plot temperature profiles

        if(ftype==1)
            plot(axes_handle, x, Th, 'r>', 'LineWidth', 0.5);
            plot(axes_handle, x, Tc, 'b>', 'LineWidth', 0.5);
        else
            plot(axes_handle, x, Th, 'r>', 'LineWidth', 0.5);
            plot(axes_handle, x, Tc, 'b<', 'LineWidth', 0.5);
        end

        % Open a new figure to display results

        results_fig = figure('Position', [200, 100, 600, 400], 'Name', 'Calculation Results');
        uicontrol('Style', 'text', 'Position', [100, 350, 400, 40], 'String', 'Heat Exchanger Calculation Results', 'FontSize', 14);
        uicontrol('Style', 'text', 'Position', [200, 250, 200, 40], 'String', sprintf('Total Heat Transfer: %.4f m^2', A_heat_exchanger), 'FontSize', 10);
        uicontrol('Style', 'text', 'Position', [200, 200, 200, 40], 'String', sprintf('Heat Duty: %.2f W', Q), 'FontSize', 10);
        uicontrol('Style', 'text', 'Position', [200, 150, 200, 40], 'String', sprintf('Overall Heat Transfer Coefficient : %.2f W/m^2.K', U), 'FontSize', 10);

end