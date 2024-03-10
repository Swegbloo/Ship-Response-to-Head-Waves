clc
clear all;
% Load offset table data (assuming it is in CSV format with columns for station, waterline, and offsets)
offsetTable = readtable('wigley_new.csv'); % Update the filename accordingly
nw = (size(offsetTable, 2))-1;
ns = (size(offsetTable, 1))-1;
rho = 1025;
g = 9.81;

% Extract data from the table
ST = zeros(ns,1);
WL = zeros(nw,1);
for i = 2:ns+1
    ST(i-1) = table2array(offsetTable(i,1));
end
for i = 2:nw+1
    WL(i-1) = table2array(offsetTable(1,i));
end
% 7th waterline is considered for the design draft
h = (ST(ns) - ST(1)) / (ns-1); % Calculate the step size
offsets = table2array(offsetTable(2:end, end));
Awp = single(2*(h/3)*(offsets(1) + 4*sum(offsets(2:2:end-1)) + 2*sum(offsets(3:2:end-2)) + offsets(end)));
%Awp = single(2*h/3 * (offsetTable(2,nw+1) + 4*sum(offsetTable((3:2:ns),nw+1)) + 2*sum(offsetTable((4:2:ns-1),nw+1)) + offsetTable(ns+1,nw+1)));
C33 = rho * g * Awp;

fprintf('Waterplane Area using simpson rule: %.2f square units\n', Awp);
fprintf('C33 : %.2f\n', C33);

xintoy = zeros(ns,1);
for i = 1:ns 
    xintoy(i) = (offsets(i))*(ST(i));
end

xysum = single(2*(h/3)*(xintoy(1) + 4*sum(xintoy(2:2:end-1)) + 2*sum(xintoy(3:2:end-2)) + xintoy(end)));
lcf = xysum/Awp;

Istrip = zeros(ns,1);
for i = 1:ns 
    Istrip(i) = (offsets(i))*((ST(i)-lcf)^2);
end

GMlintodisplacement = single(2*(h/3)*(Istrip(1) + 4*sum(Istrip(2:2:end-1)) + 2*sum(Istrip(3:2:end-2)) + Istrip(end)));
C55 = GMlintodisplacement * rho * g;
fprintf('C55 : %.2f\n', C55);
C35 = -rho * g * xysum;
C53 = C35;
fprintf('C35 : %.2f\n', C35);
fprintf('C53 : %.2f\n', C53);