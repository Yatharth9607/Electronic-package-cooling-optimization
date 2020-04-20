% For gradient based optimization
%Resistance Network Function

   function [total_cost,dT_chip,shear] = res_network_gradient_scaled(Nb,p,db,hb)
    
    %Input Variables
    Nb = Nb*1500;
    p = p*0.0002769;
    db = db*0.0005073;
    hb = hb*5.75044e-05;
    %Solder balls
    kb = 50;
%    hb = 0.0004;

    %Fixed Parameters

    %PCB
    lp = 0.1;
    wp = 0.15;
    tp = 0.0016;
    kp = 0.343;

    %Substrate
    ts = 0.0005;
    ks = 10;

    %Silicon Chip
    lc = 0.015;
    wc = lc;
    tc = 0.0007;
    kc = 148;

    %Dependent Variables

    %Substrate
    ls = (sqrt(Nb)*db) + ((sqrt(Nb)+1)*(p - db));
    if (ls<lc)
        ls = lc + 0.0014;
    end
    ws = ls;

    %Mold
    lm = ls - 0.0014;
    wm = lm;
    tm = tc + 0.0004;
    km = 1.3;

    %Adhesive
    lad = lc;
    wad = wc;
    tad = 0.0001;
    kad = 2.45;

    %Air
    T_amb = 25;
    hcp = 25;
    hcm = 50;
    ka = 0.027;

    %R1 - PCB Resistance
    R1 = tp/(kp*lp*wp);

    %R2 - Radiation from PCB to ambient (Base is not considered) - loop to
    %calculate "hrp"
    hrp = 10;
    R2 = 1/(hrp*((tp*2*(lp + wp)) + ((lp*wp) - (ls*ws))));

    %R3 - Convective from PCB to ambient (Base is not considered)
    R3 = 1/(hcp*((tp*2*(lp + wp)) + ((lp*wp) - (ls*ws))));

    %R4 - Spreading resistance from Ball layer to PCB (Construct a loop to calculate the R4 from the variable values)
    kper_4 = (hb + tp)/((hb/ka) + (tp/kp));
    kpar_4 = ((ka*hb) + (kp*tp))/(hb + tp);
    alpha_4 = sqrt(kper_4/kpar_4);
    a_4 = (p - db)/2;
    b_4 = p/2;
    e_4 = (a_4/b_4);
    lamc_4 = pi + 1/(sqrt(pi)*e_4);
    tau_4 = tp/b_4;
    Bi_4 = 1/(pi*kb*b_4*(1/(hcp*4*b_4*b_4)));

    R4 = 0.5*alpha_4*(1 - e_4)^(1.5)*((tanh(lamc_4*tau_4/alpha_4) + lamc_4/Bi_4)/(1 + (lamc_4/Bi_4)*(tanh(lamc_4*tau_4/alpha_4))))/(sqrt(pi)*kper_4*a_4);

    %R5 - Spreading resistance from single ball to PCB (Construct a loop to calculate the R5 from the variable values)
    kper_5 = (hb + tp)/((hb/ka) + (tp/kp));
    kpar_5 = ((ka*hb) + (kp*tp))/(hb + tp);
    alpha_5 = sqrt(kper_5/kpar_5);
    a_5 = db/2;
    b_5 = p/2;
    e_5 = (a_5/b_5);

    R5 = alpha_5*(1 - 1.40978*e_5 + 0.34406*(e_5)^3 + 0.0435*(e_5)^5 + 0.02271*(e_5)^7)/(4*kper_5*a_5*Nb);

    %R6 - Ball resistance
    R6 = hb/(kb*(pi*((db/2)^2 - (hb - (db/2))^2)));

    %R7 - Resistance of air between balls (Area of air needs to be defined as a function of number of balls and pitch)
    Aa = ls*ws - (Nb*pi*(db^2)/4);
    R7 = hb/(ka*Aa); %ka = thermal conductivity of air and Aa = Area between balls

    %R8 - Resistance of Bottom mask neglected
    R8 = 0;

    %R9 - Constriction resistance from substrate to Ball layer (Construct a loop to
    %calculate the R9 from the variable values)
    kper_9 = (hb + ts)/((hb/ka) + (ts/ks));
    kpar_9 = ((ka*hb) + (ks*ts))/(hb + ts);
    alpha_9 = sqrt(kper_9/kpar_9);
    a_9 = (p - db)/2;
    b_9 = p/2;
    e_9 = (a_9/b_9);
    lamc_9 = pi + 1/(sqrt(pi)*e_9);
    tau_9 = hb/b_9;
    Bi_9 = 1/(pi*ka*b_9*(1/(hcp*4*b_9*b_9)));

    R9 = 0.5*alpha_9*(1 - e_9)^(1.5)*((tanh(lamc_9*tau_9/alpha_9) + lamc_9/Bi_9)/(1 + (lamc_9/Bi_9)*(tanh(lamc_9*tau_9/alpha_9))))/(sqrt(pi)*kper_9*a_9);

    %R10 - Constriction resistance from substrate to single ball
    kper_10 = (hb + ts)/((hb/ka) + (ts/ks));
    kpar_10 = ((ka*hb) + (ks*ts))/(hb + ts);
    alpha_10 = sqrt(kper_10/kpar_10);
    a_10 = p/2;
    b_10 = db/2;
    e_10 = (a_10/b_10);
    R10 = alpha_10*(1 - 1.40978*e_10 + 0.34406*(e_10)^3 + 0.0435*(e_10)^5 + 0.02271*(e_10)^7)/(4*kper_10*a_10*Nb);

    %R11 - Substrate Resistance
    R11 = ts/(ks*ls*ws);

    %R12 - Resistance of Top mask neglected
    R12 = 0;

    %R13 - Resistance of adhesive
    R13 = tad/(kad*lad*wad);

    %R14 - Resistance of silicon chip
    R14 = tc/(kc*lc*wc);

    %R15 - Spreading resistance from silicon chip to substrate
    kper_15 = (tc + ts)/((tc/kc) + (ts/ks));
    kpar_15 = ((kc*tc) + (ks*ts))/(tc + ts);
    alpha_15 = sqrt(kper_15/kpar_15);
    a_15 = lc/2;
    b_15 = ls/2;
    e_15 = (a_15/b_15);
    lamc_15 = pi + 1/(sqrt(pi)*e_15);
    tau_15 = ts/b_15;
    Bi_15 = 1/(pi*ks*b_15*(1/(hcp*4*b_15*b_15)));

    R15 = 0.5*alpha_15*(1 - e_15)^(1.5)*((tanh(lamc_15*tau_15/alpha_15) + lamc_15/Bi_15)/(1 + (lamc_15/Bi_15)*(tanh(lamc_15*tau_15/alpha_15))))/(sqrt(pi)*kper_15*a_15);

    %R16 - Mold resistance
    R16 = tm/(km*lc*wc); %Area of chip is to be considered for thermal resistance R16

    %R17 - Convective resistance from mold to ambient
    R17 = 1/(hcm*lm*wm);    %Increased hcm for considering the effect of heat spreader

    %R18 - Radiative resistance from mold to ambient
    hrm = 10;
    R18 = 1/(hrm*lm*wm);

    %R19 - Mold resistance around the die
    R19 = (tm + tc + tad)/(km*((lm*wm) - (lc*wc)));

    %R20 - Convective resistance around the die surface to ambient
    R20 = 1/hcm*(tm + tc + tad)*2*(lm + wm);

    %R21 - Radiative resistance around the die surface to ambient
    R21 = 1/hrm*(tm + tc + tad)*2*(lm + wm);


    %----------------------------------------------------------%

    %Resistance Network formation

    %Equivalent Resistances
    Req_1 = ((1/R1) + (1/R2) + (1/R3))^(-1);
    Req_2 = ((1/R4) + (1/R5))^(-1);
    Req_3 = ((1/R6) + (1/(R7 + R8)))^(-1);
    Req_4 = ((1/R9) + (1/R10))^(-1) + R11 + R12;
    Req_5 = Req_1 + Req_2 + Req_3 + Req_4;
    Req_6 = ((1/R17) + (1/R18))^(-1);
    Req_7 = R13 + R14 + R15;

    %Heat rate
    q_1 = 2;
    syms q_2 q_3 q_4 q_5 q_6 T_chip T_sb
    eqn1 = q_2 + q_3 == q_1;
    eqn2 = q_4 + q_5 == q_2;
    eqn3 = q_3 + q_4 == q_6;
    eqn4 = -q_2*R16 - q_4*R19 + q_3*Req_7 == 0;
    eqn5 = -q_5*Req_6 + q_4*R19 + q_6*Req_5 == 0;
    eqn6 = T_amb + q_5*Req_6 + q_2*R16 == T_chip;
    eqn7 = T_amb + q_6*(Req_1 + Req_2) == T_sb;

    [A,B] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7], [q_2, q_3, q_4, q_5, q_6, T_chip, T_sb]);
    q = linsolve(A,B);
    dq = double(q);
    %disp(dq)

    %---------------------------------------------------------------------------------------------------------------------------------%

    %Thermal Stress analysis

    %Input Variables
    dT = dq(7,1) - T_amb;
    dT_chip = dq(6,1) - T_amb;
    %disp(dT_chip);
    %Substrate
    v_1 = 0.28; %  Poisson's ratio
    e_1 = 2.3e9; % Young's modulus
    a_1 = 15e-6; % coefficient of thermal expansion

    %Silicon Chip
    v_0 = 0.24; %  Poisson's ratio
    e_0 = 13e9; % Young's modulus
    a_0 = 2.4e-6; % coefficient of thermal expansion

    %Mold
    v_2 = 0.35; %  Poisson's ratio
    e_2 = 0.8e9; % Young's modulus
    a_2 = 18e-6; % coefficient of thermal expansion

    %Design with Organic Lid
    %Mid portion stresses - Axial compliances

    lamda_1 = (1 - v_1)/(e_1*ts); %Substrate
    lamda_2 = (1 - v_2)/(e_2*(tm-tc)); %Lid(Mold)
    lamda_0 = (1 - v_0)/(e_0*tc); %Chip

    D = lamda_0*lamda_1 + lamda_2*lamda_1 + lamda_0*lamda_2; %Determinant

    %Forces acting in mid cross-section of a long assembly:

    Tm_0 = ((a_0 - a_1)*lamda_2 + (a_0 - a_2)*lamda_1)*(dT/D); %In the chip
    Tm_1 = ((a_1 - a_0)*lamda_2 + (a_1 - a_2)*lamda_0)*(dT/D); %In the substrate
    Tm_2 = ((a_2 - a_1)*lamda_0 + (a_2 - a_0)*lamda_1)*(dT/D); %In the lid

    %Normal Stresses in the assembly components

    sm_0 = Tm_0/tc; %Max normal stress in chip, kg/m2
    sm_1 = Tm_1/ts; %Max normal stress in substrate
    sm_2 = Tm_2/(tm-tc); %Max normal stress in lid

    %Peripheral portion stresses - Axial compliances
    v_0 = v_2;
    e_0 = e_2;
    a_0 = a_2;
    lamda_1 = (1 - v_1)/(e_1*ts); %Substrate
    lamda_2 = (1 - v_2)/(e_2*(tm-tc)); %Lid(Mold)'s upper part
    lamda_0 = (1 - v_0)/(e_0*tc); %Lid(Mold)'s lower part

    D = lamda_0*lamda_1 + lamda_2*lamda_1 + lamda_0*lamda_2; %Determinant

    %Forces acting in mid cross-section of a long assembly:

    Tp_0 = ((a_0 - a_1)*lamda_2 + (a_0 - a_2)*lamda_1)*(dT/D); %In the chip
    Tp_1 = ((a_1 - a_0)*lamda_2 + (a_1 - a_2)*lamda_0)*(dT/D); %In the substrate
    Tp_2 = ((a_2 - a_1)*lamda_0 + (a_2 - a_0)*lamda_1)*(dT/D); %In the lid

    %Normal Stresses in the assembly components

    sp_0 = Tp_0/tc; %Max normal stress in lid's lower part
    sp_1 = Tp_1/ts; %Max normal stress in substrate
    sp_2 = Tp_2/(tm-tc); %Max normal stress in lid's upper part

    %Solder joints - Stiffness analysis
    stiff = 0.37e6; % stiffness coefficient of solder ball
    F = (stiff) * (a_1 * dT * ls);
    shear = F /((pi*((db/2)^2 - (hb - (db/2))^2)));
    %disp(shear)

    %----------------------------------------------------------------------------%

    %Objective Function - COST

    %{
    Mass calculation:
    M_sub = rhos*ls*ws*ts;      %Substrate
    M_mold = rhom*lm*wm*tm;     %Mold
    M_sb = rhob*Nb*(0.44e-9);   %Solder
    %}

    %Material Cost

    % substrate
    cost_sub = (15.31e-3)*ls*ws*1e6 ;

    % mold 
    Vol_mold = (lm*wm*tm) - (lc*wc*(tc + tad)); 
    cost_mold = (Vol_mold)*0.00003963011*1e9; 

    %solder
    cost_solder1 = 2e8*Nb*pi*((db/2)^2)*hb; cost_solder2 = 1e7*(p^2); cost_solder3 = 0.004/(hb); cost_solder4 = 0.09*Nb;
    cost_solder = cost_solder1 + cost_solder2 + cost_solder3 + cost_solder4;
    % total cost
    total_cost = cost_sub + cost_mold + cost_solder;
   end