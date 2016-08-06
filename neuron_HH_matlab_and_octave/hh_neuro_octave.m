clc;
clear;

% Code is based on the set of governing equations of Hodgkin-Huxley mode

% See wiki 1) http://en.wikipedia.org/wiki/Hodgkin-Huxley_model
% Original work 2) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1392413/
% Book 3) http://nelson.beckman.illinois.edu/courses/physl317/part1/Lec3_HHsection.pdf

% Time interval
dt = 0.01; % Time Step [ms]
tmax = 7; % Max time [ms]

% Constants
Cm = 0.01; % Membrane capacity [uF/cm^2]
I_e = 0.5; % Input external current [nA]
% Reversal potentials: Na, K, leak [mV]
ENa = 55;
EK = -70;
EL = -50;
% Base conductances: Na, K, leak [ms/cm^2]
g_Na = 1.5;
g_K = 0.4;
g_L = 0.004;

% Allocate memory
t = 0:dt:tmax;
max = length(t);
Vm = zeros(1, max);
m = zeros(1, max);
n = zeros(1, max);
h = zeros(1, max);
gNa = zeros(1, max);
gK = zeros(1, max);
gL = zeros(1, max);
INa = zeros(1, max);
IK = zeros(1, max);
IL = zeros(1, max);
Ie = ones(1, max);
Ie = Ie .* I_e; % Constant current
%Ie = sin(10*t).^8 .* I_e; % Spiking current

% Set initial conditions
    
Vm(1) = -60; % Membrane potential
v = Vm(1);
an = 0.01 * (v + 50) / (1 - exp(-(v + 50)/10));
bn = 0.125 * exp(- (v+60) / 80);
am = 0.1 * (v + 35) / (1 - exp(-(v + 35)/10));
bm = 4 * exp(- (v+60) / 18);
ah = 0.07 * exp(- (v+60) / 20);
bh = 1 / (exp(-(v + 30)/10) + 1);
n(1) = an / (an + bn); % potassium channel activation
m(1) = am / (am + bm); % sodium channel activation
h(1) = ah / (ah + bh); % sodium channel inactivation

for i=1:1:max-1
    % Compute alpha and beta (ion counts) for n, m,h
    v = Vm(i);
    an = 0.01 * (v + 50) / (1 - exp(-(v + 50)/10));
    bn = 0.125 * exp(- (v+60) / 80);
    am = 0.1 * (v + 35) / (1 - exp(-(v + 35)/10));
    bm = 4 * exp(- (v+60) / 18);
    ah = 0.07 * exp(- (v+60) / 20);
    bh = 1 / (exp(-(v + 30)/10) + 1);
    
    % Forward Euler integration for n,m,h
    n(i+1) = n(i) + ( an * (1 - n(i)) - bn * n(i) ) * dt;
    m(i+1) = m(i) + (am * (1 - m(i)) - bm * m(i)) * dt;
    h(i+1) = h(i) + (ah * (1 - h(i)) - bh * h(i)) * dt;
    
    % Conductances: Na, K, leak
    gNa(i) = g_Na * m(i)^3 * h(i);
    gK(i) = g_K * n(i)^4;
    gL(i) = g_L;
    
    % Currents:  Na, K, leak
    INa(i) = gNa(i) * (Vm(i) - ENa);
    IK(i) = gK(i) * (Vm(i) - EK);
    IL(i) = gL(i) * (Vm(i) - EL);
    
    % Forward Euler integration for Vm
    Vm(i+1) = Vm(i) + (Ie(i) - (INa(i) + IK(i) + IL(i))) / Cm * dt;
end

% Plots settings
rows = 1;
cols = 4;
x_axis_name = 't [ms]';

% Voltage
subplot(cols, rows, 1);

plot(t,Vm); 
grid on;
legend('V');
xlabel(x_axis_name);
ylabel('Voltage [mV]');
title('The Hodgkin-Huxley Model');

% Conductances
subplot(cols, rows, 2);

plot(t,gNa, t,gK);
grid on;
legend('g_{Na}', 'g_{K}');
xlabel(x_axis_name);
ylabel('Conductances [ms/cm^2]');

% Currents
subplot(cols, rows, 3);

plot(t,INa, t,IK, t,Ie);
grid on;
legend('I_{Na}', 'I_{K}', 'I_{external}');
xlabel(x_axis_name);
ylabel('Currents [nA]');

% n, m, h
subplot(cols, rows, 4);

plot(t,n, t,m, t,h);
legend('n (Na act)', 'm (K act)', 'h (K inact)');
grid on;
xlabel(x_axis_name);

% Save an image
print('-dpng','-r600', 'hh_neuron_result.png');
