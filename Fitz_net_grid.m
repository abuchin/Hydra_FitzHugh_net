%% INTEGRATION PARAMETERs
clear;       % remove previous varibles
close all;

T=51;        % total time, mS
dt=0.01;     % time step, ms

Df=1/dt;            % Delta function approximation

Tframe=1;    % movie
dTframe=1; % ms
frame=1;

figure('units','normalized','outerposition',[0 0 0.8 0.8]); % show figure window
%%

%% NETWORK PARAMETERS
N=100;        % Number of neurons
%%

%% CONNECTIVITY
gEE_mean=0.3;  % mA/cm^2 (current synapses)
tau_EE=5.4;       % ms, AMPA current 5.4

%load('Adjacency_20X5.mat');  % load RP2 connectivity matrix
%S_EE=A_tube*gEE_mean;

A_tube=M_tor(sqrt(N));
S_EE=A_tube*gEE_mean;

% GAP junctions, different locations than synapses
g_GAP=0.1;
l=zeros(N,4);
%S_GAP=M_anytube(20,5);

%
S_GAP=S_EE;

for i=1:1:N 
   l(i,1:length(find(S_GAP(i,1:end)>0)))=find(S_GAP(i,1:end)>0);
end

gap_diff_coeff=sum(S_GAP);

%S_GAP=g_GAP*S_GAP;
%}

%%

%%  Graph plot details

A_graph=graph(A_tube);
% tube parameters
m=10;    % wideness
n=10;   % length

[Y,X]=meshgrid(1:n,1:m);        
x=reshape(X,[1,m*n]);
y=reshape(Y,[1,m*n]);
%%


%% NEURON PARAMETERS
% Fitz-Hugh model in the original form, see the implementation
% http://www.scholarpedia.org/article/FitzHugh-Nagumo_model

% dV/dt = V - V^3/3 -W + I
% dW/dt = 0.08(V+0.7-0.8W)


V_spike=0;     % mV, threshold for spike
t_ref=15;      % refractory period for spike recordings -> spike processing part
%%

%% ICs
%Constant ICs
%
V(1:N)=-1.1994;   % mV
W(1:N)=-0.6243;

%REPRESENTATIVE CELLS
repr=round(45);             % number of representative neuron
V_repr=zeros(1);    % mV

% SYNAPTIC VARIABLEs
I_EE=zeros(N,1);

% GAP-JUNCTION VARIABLE
I_GAP=zeros(N,1);

% FIRINGS
firings=[];           % spike timings
fired=[];             % indexes of spikes fired at t
V_sp=10*ones(N,1);   % vector of times of elements that did spike
fired_delta(1:N)=0;   % vector of fired delta-function
%%

%% STIMULATION
IN(1:N)=0;        %mA/cm^2
% stimulated neuron
IN(repr)=0.5;        %mA/cm^2

ts=5;           % end of stimulation, ms
%%

%% TIME INTEGRATION LOOP
for t=1:1:round(T/dt)
    
% Voltage
dVdt=V-(V.^3)./3-W+IN +I_EE' +I_GAP';
dWdt=0.08*(V+0.7-0.8.*W);

%dVdt=1/C*(-gL.*(V-VL) -gNa.*m.^3.*h.*(V-VNa)-gK*n.^4.*(V-VK) +IN +I_EE');
V=dVdt*dt + V;
W=dWdt*dt + W;

% Representative cell
V_repr(t)=V(repr);
I_repr_ext(t)=IN(repr);
I_repr_EE(t)=I_EE(repr);

% Synaptic variable
I_EE=(-I_EE/tau_EE + Df.*sum(S_EE(:,fired),2) )*dt + I_EE;

% GAP-junction current (no generalized for open borders)
I_GAP=g_GAP*(sum(V(l),2)-4*V');                            
%I_GAP=(g_GAP*(sum(repmat(V,N,1).*(S_GAP),1)-4*V))'; % (works for arbitrary graph, does not work propely)

% Stimulation

if t*dt>ts
    IN(repr)=0;
end

% Firings processing

% Debug
%if isempty(find(V>V_spike))==0
%   a=1;
%end

fired=[];
fired=find(V>V_spike);

INT_E=(intersect(find(V_sp<t_ref),fired))';   % intersection V_sp<2 and fired_E
fired=setxor(fired,INT_E);                    % remove fired_E elements that intersect with VSOMA_sp<2

V_sp=V_sp + dt;                               % update the t* vector
V_sp(fired)=0;
            
t_fired(1:length(fired))=t;
        
% check the format!!!
if isrow(fired)==0
     fired=fired';
end
    
if isrow(t_fired)==0
     t_fired=t_fired';
end
    
spikes=horzcat(t_fired',fired');   
firings=vertcat(firings,spikes);

t_fired=[];
spikes=[];

fired_delta(:)=0;           % vector of current delta functions of all neurons
fired_delta(fired)=Df;      

%}
%end
%%

%% FINAL PLOT / ONLINE PLOT

if (t-1)*dt==Tframe
    
subplot(2,2,1);
if isempty(firings)==0
plot(firings(:,1)*dt,firings(:,2),'.','MarkerSize',15);
end
ylabel('Cell index');
ylim([0 N]);
set(gca,'FontSize',20);             % set the axis with big font
title(sprintf('Fitz-Hugh population, T=%d ms',round(Tframe)));
set(gca,'FontSize',20);             % set the axis with big font
box off;

subplot(2,2,2);
plot((1:t)*dt,V_repr);
ylabel('Voltage (mV)');
set(gca,'FontSize',20);             % set the axis with big font
title('Representative neuron');
set(gca,'FontSize',20);             % set the axis with big font
box off;

subplot(2,2,3);
s=plot(A_graph,'Xdata',x,'Ydata',y);    %handle for the graph plot
set(gca, 'CLim', [-2, 2]);
s.NodeCData = V;
s.MarkerSize=7;
title('Voltage distribution');
set(gca,'Fontsize',15);
colormap('jet');
colorbar;
box off;

%Plot voltage on rectangular coordinates
%{
imagesc((reshape(V,sqrt(N),sqrt(N)))',[-2 2]); % Voltage distribution
set(gca,'Ydir','normal');
set(gca,'FontSize',20);             % set the axis with big font
title('Voltage distribution (mV)');
colorbar;
box off;
%}

subplot(2,2,4);
q=plot(A_graph,'Xdata',x,'Ydata',y);    %handle for the graph plot
set(gca, 'CLim', [0, 2*gEE_mean]);
q.NodeCData = IN +I_EE';
q.MarkerSize=7;
title('Synaptic input');
set(gca,'Fontsize',15);
colormap('jet');
colorbar;
box off;


% Plot input on rectangular coordinates
%{
imagesc((reshape((IN +I_EE'),sqrt(N),sqrt(N)))',[0 4*gEE_mean]); % Voltage distribution
set(gca,'Ydir','normal');
set(gca,'FontSize',20);             % set the axis with big font
colorbar;
title('Synaptic input (\muA/cm^2)');
box off;
%}

MOV(frame)=getframe(gcf);

frame=frame+1;          % counter for the movie
Tframe=Tframe + dTframe;
%%

end

end

% Save the movie
%movie2avi(MOV,'Fitz_net.avi','fps',10,'quality',1);

%%