clear;
clc;

%Declaring Macros
global r2p;
r2p = @(x) [abs(x) rad2deg(angle(x))]; % Rectangular -> Phasor


%Declaring Variables
global f; global v; global l; global p; global rdc; global rac;
global d; global t; global k; global ph; global gmd; global gmr;
global tw; global lx; global cx;%#ok<*GVMIS> 
f=0; v=0; l=0; p=0; rdc=0; rac=0; d=0; t=0; k=0;
ph=0; gmd=0; gmr=0; tw=0; lx=0; cx=0;

global pf; global pa; global n; global Ir; global Ic;
pf=0; pa=0; n=0; Ir=0; Ic=0;

global Z; global Y;
Z=0; Y=0;


global A; global B; global C; global D;
A=0; B=0; C=0; D=0;


% Arrangment
global arr;
arr = [20,1];

        %Freq   Vol     Power   Length  D       Phases  C Bun#  GMD     GMR    
        %TempC  Resis   Induc   Cap     Imp     Adm     Param   PF      PA
        %RedC   Current
arr =   [1,     2,      3,      4,      5,      6,      7,      8,      9,...
        10,     11,     12,     13,     14,     15,     16,     17,     18,...
        19,     20];
     


% Getting TL Type & Model
disp("Power Elements Systems Coding Project");
disp("=============================================")
disp("=========== First Choose TL Model ===========")
disp("=============================================")
disp("1- Short TL")
disp("2- Medium TL PI-Model")
disp("3- Medium TL T-Model")

global c;

c = input(">> ");

while (~(c<=3 && c>=1))
    c = input("(1-3) >> ");
end

if c==1
    for i=1:length(arr)
        dispShort();
        calculation(i);
    end
elseif c==2
    for i=1:length(arr)
        dispMedium();
        calculation(i);
    end
elseif c==3
    for i=1:length(arr)
        dispMedium();
        calculation(i);
    end
end

function calculation(i)
    global c; global arr;
    if i==arr(1)
        disp("================= Frequency =================");
        disp("=============================================");
        getFreq();
    elseif i==arr(2)
        disp("================= Line Voltage ==============");
        disp("=============================================");
        getVol();
    elseif i==arr(3)
        disp("================= Load Power ================");
        disp("=============================================");
        getPow();
    elseif i==arr(4)
        disp("================= TL Length =================");
        disp("=============================================");
        if c==1
            getLen("S");
        else
            getLen("M");
        end
    elseif i==arr(5)
        disp("================= Diameter ==================");
        disp("=============================================");
        getDiam();
    elseif i==arr(6)
        disp("================= Phases # ==================");
        disp("=============================================");
        getPhase();
    elseif i==arr(7)
        disp("================= Conducter Bun. # ==========");
        disp("=============================================");
        getK();
    elseif i==arr(8)
        disp("================= GMD =======================");
        disp("=============================================");
        calculateGMD();
    elseif i==arr(9)
        disp("================= GMR =======================");
        disp("=============================================");
        calculateGMR();
    elseif i==arr(10)
        disp("================= Temp. Const. ==============");
        disp("=============================================");
        getTempC();
    elseif i==arr(11)
        disp("================= Resistance ================");
        disp("=============================================");
        getR();
    elseif i==arr(12)
        disp("================= Inductance ================");
        disp("=============================================");
        calculateInd();
    elseif i==arr(13)
        if ~(c==1)
            disp("================= Capacitance ===============");
            disp("=============================================");
            calculateCap();
        end
    elseif i==arr(14)
        disp("================= Impedence =================");
        disp("=============================================");
        calculateImp();
    elseif i==arr(15)
        if ~(c==1)
            disp("================= Admittance ================");
            disp("=============================================");
            calculateAdm();
        end
    elseif i==arr(16)
        disp("================= Parameters ================");
        disp("=============================================");
        calculateParm();
    elseif i==arr(17)
        disp("================= Power Factor ==============");
        disp("=============================================");
        getPF();
    elseif i==arr(18)
        disp("================= Power Angle ===============");
        disp("=============================================");
        getPA();
    elseif i==arr(19)
        disp("================= Redundant Circuit =========");
        disp("=============================================");
        calculateN();
    elseif i==arr(20)
        disp("================= Current  ==================a");
        disp("=============================================");
        calculateI();
    end
end

%Declaring functions;

function calculateI()
    global Ir; global Ic; global p; global v; global pf; global n; global k;
    Ir=p*1e6/(v*1e3*sqrt(3)*pf);
    Ic=Ir/((n-1)*k);

    fprintf("\nCurrent Ref. calculated = %f A\n",Ir);
    fprintf("\nCurrent Cap. calculated = %f A\n",Ic);
    disp("Press enter to continue")
    pause;
end

function calculateN()
    global n; global lx; global p; global v; global pa;
    n = p*1e6*lx/sin(pa*pi/180)./(v*1e3)./(v*1e3);
    n = ceil(n);
    a = input("Please input number of redundant circuits >> ");
    while (a<0)
        a = input("Please input number of redundant circuits (<0) >> ");
    end
    n = n+a;

    fprintf("\nRedundant circuit calculated = %d\n",n);
    disp("Press enter to continue")
    pause;
end

function calculateParm()
    global A; global B; global C; global D; global Z; global Y;
    global c;

    global r2p;

    if c==1
        A=1; D=1; C=0;
        B=Z;
    elseif c==2
        A = 1 + (Z*Y)/2;
        B = Z;
        C = Y*(1 + (Z*Y)/4);
        D = 1 + (Z*Y)/2;
    else
        A = 1 + (Z*Y)/2;
        B = Z*(1 + (Z*Y)/4);
        C = Y;
        D = 1 + (Z*Y)/2;
    end

    fprintf("\nParameters calculated");
    fprintf("\nA = %f%+fj = %.2f<%.2f", real(A), imag(A), r2p(A));
    fprintf("\nB = %f%+fj = %.2f<%.2f Ohm", real(B), imag(B), r2p(B));
    fprintf("\nC = %f%+fj = %.2f<%.2f Siemens", real(C), imag(C), r2p(C));
    fprintf("\nD = %f%+fj = %.2f<%.2f\n", real(D), imag(D), r2p(D));

    disp("Press enter to continue")
    pause;
end

function calculateAdm()
    global cx; global f; global l; global Y;
    Y = (2*pi*f*cx*l)*1i;

    fprintf("\nAdmetance calculated = %.2e Siemens\n", imag(Y));
    disp("Press enter to continue")
    pause;
end

function calculateImp()
    global lx; global f; global rac; global l; global Z;
    lx = 2*pi*f*l*lx*1e3;
    Z = rac*l*1e-3 + (lx)*1i;

    fprintf("\nImpedence calculated = %.3f%+.3fj ohm\n", real(Z), imag(Z));
    disp("Press enter to continue")
    pause;
end

function calculateCap()
    global cx; global gmd; global gmr; global f; global v;

    cx = (2*pi*8.854e-12) / (log(gmd/gmr));

    fprintf("\nCapacitance calculated = %.2e F/m", cx);
    fprintf("\nReactive power calculated = %.2e VAR", 2*pi*f*cx*(v*1e3)^2);
    fprintf("\nCharging current calculated = %.2e A/m \n", 2*pi*f*cx*(v*1e3)/(sqrt(3)));
    disp("Press enter to continue")
    pause;
end

function calculateInd()
    global lx; global gmd; global gmr; global k; global ph;

    if ph == 1
        lx = 4e-7 * log(gmd/(gmr*(0.7788)^(1/k)));
    elseif ph == 3
        lx = 2e-7 * log(gmd/(gmr*(0.7788)^(1/k)));
    end

    fprintf("\nInductance calculated = %f H/m \n", lx);
    disp("Press enter to continue")
    pause;
end

function calculateGMR() 
    fprintf("\n");
    global gmr; global k; global d;
    d = d*1e-3;
    if k==1
        gmr = d/2;
    elseif k==2
        s=input("Please enter distance between bundle in m >> ");
        while (s<=0)
            s=input("Please enter distance between bundle in m (>0) >> ");
        end
        gmr = ((d/2)*s^(k-1))^(1/k);
    else
        c = input("Are the TL bundles symmertical ? (y/n) >> ", "s");
        while (~(c=="y" || c=="n"))
            c = input("Are the TL bundles symmertical ? (y/n) >> ", "s");
        end

        if c=="y"
            s=input("Please enter distance between bundle in m >> ");
            while (s<=0)
                s=input("Please enter distance between bundle in m (>0) >> ");
            end
            if k==3
                gmr = ((d/2)*s^(k-1))^(1/k);
            else
                gmr = 1.091*((d/2)*s^(k-1))^(1/k);
            end
        else 
            disp("They need to be symmetrical to continue !");
            quit(1);
        end
    end
    fprintf("\nGMR calculated = %.2f m\n", gmr);
    disp("Press enter to continue")
    pause;
end

function calculateGMD()
    fprintf("\n");
    global gmd; global ph;
    if ph==1
        gmd = input("Please enter distance between phases in m >> ");
        while (gmd <= 0)
            gmd = input("Please enter distance between phases in m (>0) >> ");
        end
    else
        c = input("Are the TL wires symmertical ? (y/n) >> ", "s");
        while (~(c=="y" || c=="n"))
            c = input("Are the TL wires symmertical ? (y/n) >> ", "s");
        end
        if c=="y"
            gmd = input("Please enter distance between phases in m >> ");
            while (gmd <= 0)
                gmd = input("Please enter distance between phases in m (>0) >> ");
            end
        else
            s1 = input("Please enter distance between line a & b in m >> ");
            while (s1<=0)
                s1 = input("Please enter distance between line a & b in m (>0) >> ");
            end
            s2 = input("Please enter distance between line a & c in m >> ");
            while (s2<=0)
                s2 = input("Please enter distance between line a & c in m (>0) >> ");
            end
            s3 = input("Please enter distance between line b & c in m >> ");
            while (s3<=0)
                s3 = input("Please enter distance between line b & c in m (>0) >> ");
            end
            gmd = (s1*s2*s3)^(1/3);
        end
    end

    fprintf("\nGMD calculated = %.2f m\n", gmd);
    disp("Press enter to continue")
    pause;
end

function dispShort()
    clc;
    disp("=============================================");
    disp("================= Short TL ==================");
    disp("=============================================");
end

function dispMedium()
    global c;
    clc;
    disp("=============================================");
    if c==2
        disp("============= Medium TL PI-Model ============");
    else
        disp("============= Medium TL T-Model =============");
    end
    disp("=============================================");
end

function getK()
    fprintf("\n");
    global k;
    k=input("Please enter bundle number >> ");
    while (~(k<=4 && k>=1))
        k=input("Please enter bundle number (1-4) >> ");
    end

    fprintf("\nBundle number entered = %d\n", k);
    disp("Press enter to continue")
    pause;
end

function getR() 
    fprintf("\n");
    global rdc; global rac; global t; global d; global tw;
    tw = input("Please enter working temp in C >> ");
    c = input("Do you have resistivity (y/n)  >> ", "s");

    while (~(c=="y" || c=="n"))
            c = input("Do you have resistivity (y/n)  >> ", "s");
    end
    
    resis=0;

    if c=="y"
        resis = input("Please enter conductor resistivity @ 25C in mohm*km >> ");
        while (resis <= 0)
           resis = input("Please enter conductor resistivity @ 25C in mohm*km (>0) >> "); 
        end
        rdc = ((resis)/(pi*((d/2)*1e-6)^2));
    else
        rdc = input("Please enter conductor resistance @ 25C in mohm/km >> ");
        while (rdc <= 0)
           rdc = input("Please enter conductor resistance @ 25C in mohm/km (>0) >> "); 
        end
    end

    spir = input("Please enter spiraling effect in % >> ");
    while (~(spir<=100 && spir>=0))
        spir = input("Please enter spiraling effect in % (0-100) >> ");
    end

    rdc = (rdc*(t+tw)/(t+25))*(1+(spir/100));

    rac = input("Please enter Rac @ " + tw + "C in mohm/km >> ");
    while (rac <= 0)
        rac = input("Please enter Rac @ " + tw + "C in mohm/km (>0) >> ");
    end

    skin = rac/rdc;
    fprintf("Info: due to skin effect, resistance is increased by %.3f %% \n", skin);
    fprintf("Rdc calculated = %.7f mohm/km \n", rdc);
    fprintf("Rac calculated = %.7f mohm/km \n", rac);
    disp("Press enter to continue")
    pause;
end

function getPhase()
    fprintf("\n");
    global ph;
    ph=input("Please enter number of phases (1 or 3) >> ");
    while (~(ph == 3 || ph == 1))
        ph=input("Please enter number of phases (1 or 3) >> ");
    end

    fprintf("\nPhase entered = %d\n", ph);
    disp("Press enter to continue")
    pause;
end

function getPF()
    global pf;
    pf=input("Please input power factor >> "); 
    while (~(pf<=1 && pf>0))
        pf=input("Please input power factor (0-1) >> "); 
    end

    fprintf("\nPower factor entered = %.2f\n", pf);
    disp("Press enter to continue")
    pause;
end

function getPA()
    global pa;
    pa=input("Please input power angle in degree >> "); 
    while (~(pa<=90 && pa>0))
        pa=input("Please input power angle in degree (0-90) >> "); 
    end

    fprintf("\nPower angle entered = %.2f degree\n", pa);
    disp("Press enter to continue")
    pause;
end

function getTempC()    
    fprintf("\n");
    global t;
    t=input("Please enter conductor temprature constant in C' >> ");
    fprintf("\nTemprature constant entered = %d C'\n", t);
    disp("Press enter to continue")
    pause;
end

function getDiam()
    fprintf("\n");
    global d;
    d = input("Please enter conductor diameter in mm >> ");
    while (d<=0)
        d = input("Please enter conductor diameter in mm (>0) >> ");
    end

    fprintf("\nDiameter entered = %.2f mm\n", d);
    disp("Press enter to continue")
    pause;
end

function getFreq()
    fprintf("\n");
    global f;
    f=input("Please input operating frequency in Hz >> "); 
    while (f<=0)
        f=input("Please input operating frequency in Hz(>0) >> "); 
    end

    fprintf("\nFrequency entered = %d in Hz\n", f);
    disp("Press enter to continue")
    pause;
end

function getVol()
    fprintf("\n");
    global v;
    v=input("Please input line voltage in kV >> ");
    while (v<=0)
        v=input("Please input line voltage in kV (>0) >> "); 
    end

    fprintf("\nVoltage entered = %d kV\n", v);
    disp("Press enter to continue")
    pause;
end

function getLen(T)
    fprintf("\n");
    global l;
    l=input("Please input TL length in km >> "); 
    if T=="S"
        while (~(l<=80 && l>0))
            l=input("Please input TL length in km (0-80) >> "); 
        end
    else
        while (~(l<=240 && l>=80))
            l=input("Please input TL length in km (80-240) >> "); 
        end
    end

    fprintf("\nLength entered = %.2f km\n", l);
    disp("Press enter to continue")
    pause;
end

function getPow()
    fprintf("\n");
    global p;
    p=input("Please input load power in MW >> "); 
    while (p<=0)
        p=input("Please input load power in MW (>0) >> "); 
    end

    fprintf("\nPower entered = %d MW\n", p);
    disp("Press enter to continue")
    pause;
end