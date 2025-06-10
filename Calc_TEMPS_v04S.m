% -------- Function Name ---------------------------------------------- 
% Calc_TEMPS_v04S.m
% -------- Purpose --------
% Finite-difference thermal solver with non-homogeneous thermal properties, 
% Pennes perfusion model, and heating source defined by Q variable.
% -------- Input Data --------
%   Modl                [nX,nY,nZ]. Segemented tissue model (ie.  1=water, 2=fat, 3...nTypes)
%   T0                  [nX,nY,nZ]. Initial temperature condition relative to baseline. [deg C].
%   Vox                 1x3 vector. Voxel dimensions [dy dx dz]. [m].
%   dt                  Scalar. Time step for thermal solver. [s].
%   HT                  1xnFZ vector. Heating time at each FZ location. [s].
%   CT                  1xnFZ vector. Cooling time at each FZ location. [s].
%   rho                 1xnTypes. Density. [kg/m^3]. 
%   k                   1xnTypes. Thermal conductivity. [W/(m*deg C)]. 
%   c                   1xnTypes. Specific heat capacity. [J/(kg*deg C)]. 
%   wType               '1'==Uniform perfusion; '2'==Local perfusion;
%   w                       if wType=1--> 1xnTypes. Tissue type uniform Pennes perfusion. [kg/(m^3*s)]. 
%                           if wType=2--> [nX,nY,nZ]. Local voxel by voxel Pennes perfusion. [kg/(m^3*s)].
%   Q                   [nX,nY,nZ,nFZ]. Q patterns for each focal zone. [W/m^3].
%                           NOTE: Q=rho*SAR. [W/m^3]=[kg/m^3]*[W/kg].
%   nFZ                 Scalar. Number of Focal Zone Locations.
%   tacq                Scalar. Final time resolution of temperatures. [s].
%                           NOTE: 'dt' should divide evenly into 'tacq'.
%   Tb                  Scalar. Arterial blood temperature for perfusion. [deg C].
%                           NOTE: Measured relative to baseline, typically 0 deg C.
%   BC                  '0' for zero temperature boundary condition. 
%                       '1' for adiabatic boundary condition. 
%                       '2' for matching slope boundary condition.
%   temp_file           (optional) A file for writing the temperatures to.
%                       For data sets which involve a large number of time
%                       points (i.e. enough that all the memory will be
%                       consumed), writing to a file will prevent memory
%                       thrashing and speed the program up.
% -------- Output Data --------
%   TEMPS               4D finite-difference temperature array
%   time                time vector associated with TEMPS array
% -------- Required Subfunctions --------
%       NA
% -------- Author Info --------
%       Christopher Dillon
%       Department of Bioengineering
%       University of Utah
%       Version 01: 28 September 2010
% -------- Updates --------
%       v02      30 Nov 2011        % Introduced 'timeratio' variable to save data every timeratio-th FD step
%       v03      12 Apr 2014        % Introduced 'T0', nonuniform initial condition
%                                   % Replaced 'timeratio' variable with 'tacq'
%       v04      18 Jul 2014        % Commented and cleaned up code
%                                   % Reintroduced stability criterion
%                                   % Replaced 'SAR' with 'Q'
%                                   % Introduced Adiabatic and Matching Slope boundary conditions
%                                   % Introduced option for local perfusion
%       v04S     05 Jan 2015        % (Scott Almquist) improved the run time by precalculating any
%                                      constant matrices and added option to write to a file (see
%                                      helper function read_temp_file.m)
function [TEMPS,time]=Calc_TEMPS_v04S(Modl,T0,Vox,dt,HT,CT,rho,k,cp,wType,w,Q,nFZ,tacq,Tb,BC,temp_file);

dx=Vox(1);dy=Vox(2);dz=Vox(3);      % Voxel dimensions
A=dx/dy; B=dx/dz;                   % Dimensionless increment
[nX,nY,nZ] = size(Modl);            % Identify number of voxels
t_final=sum(HT)+sum(CT);            % Total treatment time to be modeled. [s].
time=0:tacq:t_final;                % Time vector
NT=t_final/dt;                      % Total number of FD time steps
nntt=length(time);                  % Total number of temperature distributions in time to save
timeratio=tacq/dt;                  % NOTE: tacq/dt should be an integer

% Calculate the maximum time step for stability of the thermal model
w_max=max(w(:)); rho_min=min(rho); cp_min=min(cp); k_max=max(k);       % Required parameters for max time step calculation
dt_max=1/(w_max/rho_min+2*k_max*(1+A^2+B^2)/(rho_min*cp_min*dx^2));  % Maximum allowable time step before iterations become unstable (s)
if dt>dt_max
    errordlg('Time step ''dt'' is too large for stable finite-difference calculations. Reduce ''dt'' and try again.','ERROR!','Modal');return
end

% ----------------------------------------
% Create Matrices of Properties
% ----------------------------------------
k1=zeros(nX,nY,nZ,'single');
k1(:,:,:)=k(Modl(:,:,:));
inv_k1 = 1/k1;

rho_m=zeros(nX,nY,nZ,'single');
rho_m(:,:,:)=rho(Modl(:,:,:));

cp_m=zeros(nX,nY,nZ,'single');
cp_m(:,:,:)=cp(Modl(:,:,:));

if wType==1
    w_m=zeros(nX,nY,nZ,'single');
    w_m(:,:,:)=w(Modl(:,:,:));
elseif wType==2
    w_m=w; clear w;
end

rho_cp=rho_m.*cp_m;                       % Simplfies later equations by combining density and specific heat

% Shift k values for use in solver
inv_k2k1= 1/circshift(k1,[1 0 0 0]) + inv_k1;  % (m*degC/W)
inv_k3k1= 1/circshift(k1,[-1 0 0 0]) + inv_k1;
inv_k4k1= 1/circshift(k1,[0 1 0 0]) + inv_k1;
inv_k5k1= 1/circshift(k1,[0 -1 0 0]) + inv_k1;
inv_k6k1= 1/circshift(k1,[0 0 1 0]) + inv_k1;
inv_k7k1= 1/circshift(k1,[0 0 -1 0]) + inv_k1;

Coeff1 = 2*dt/(rho_cp*dx^2);                % (m*degC/W)
k8 = (1/(inv_k2k1)+1/(inv_k3k1)...          % x direction conduction (W/m/degC)
       +A^2/(inv_k4k1)+A^2/(inv_k5k1)...    % y direction conduction (W/m/degC)
       +B^2/(inv_k6k1)+B^2/(inv_k7k1));     % z direction conduction (W/m/degC)
Coeff2 = (1-(w_m*dt)./rho_m-2*dt/(rho_cp*dx^2).*k8); % Changes associated with this voxel's old temperature (Unitless)
Perf=(w_m*dt*Tb)./rho_m; % Precalculate perfusion term (degC)

J=nX+1;                                 % second to last voxel in direction X
K=nY+1;                                 % second to last voxel in direction Y
L=nZ+1;                                 % second to last voxel in direction Z

% ----------------------------------------
% Solver
% ----------------------------------------
tic   % Starts the stopwatch
h = waitbar(0,'Please wait... Model temperatures are being calculated');  % Initiate waitbar
c_old=1;    % counter
% Preallocate temperature arrays
T_old=zeros(nX+2,nY+2,nZ+2,'single');   % Define Old Temperatures ***NOTE: Expanded by 2 voxels in x y and z direction
% Initial Condition
T_old(2:J,2:K,2:L)=T0;              % Fill in initial condition temperatures.

% Boundary Condition
if BC==1;                           % Adiabatic boundary condition (zero-slope)
    T_old(1,   2:K, 2:L)=T0(1,:,:);
    T_old(end, 2:K, 2:L)=T0(end,:,:);
    T_old(2:J, 1,   2:L)=T0(:,1,:);
    T_old(2:J, end, 2:L)=T0(:,end,:);
    T_old(2:J, 2:K, 1  )=T0(:,:,1);
    T_old(2:J, 2:K, end)=T0(:,:,end);
elseif BC==2;                       % Matching Slope Boundary (slope between edge and 2nd voxel matches slope between 2nd and 3rd voxel)
    T_old(1,:,:)=T_old(2,:,:)-(T_old(3,:,:)-T_old(2,:,:));
    T_old(end,:,:)=T_old(end-1,:,:)-(T_old(end-2,:,:)-T_old(end-1,:,:));
    T_old(:,1,:)=T_old(:,2,:)-(T_old(:,3,:)-T_old(:,2,:));
    T_old(:,end,:)=T_old(:,end-1,:)-(T_old(:,end-2,:)-T_old(:,end-1,:));
    T_old(:,:,1)=T_old(:,:,2)-(T_old(:,:,3)-T_old(:,:,2));
    T_old(:,:,end)=T_old(:,:,end-1)-(T_old(:,:,end-2)-T_old(:,:,end-1));
end

T_new=T_old;                            % Define new temperatures

% Initial temperature profile
use_file = 0;

if exist('temp_file', 'var')
    use_file = 1;
    fid = fopen(temp_file, 'w+');
    fwrite(fid, nX, 'single');
    fwrite(fid, nY, 'single');
    fwrite(fid, nZ, 'single');
    fwrite(fid, nntt, 'single');
    fwrite(fid, T0, 'single');
    TEMPS = 0; %dummy value to return
else
    TEMPS=zeros(nX,nY,nZ,nntt,'single');    % Final exported array
    TEMPS(:,:,:,1) = T0; 
end



for mm=1:nFZ                                % Run Model for each focal zone location
    % Generate the PowerOn vector for each focal zone location (includes heating and cooling time)
        nt=ceil(HT(mm)/dt)+ceil(CT(mm)/dt); % Number of time steps at FZ location mm
        PowerOn=zeros(nt,1);                % Zero indicates no power.
        PowerOn(1:ceil(HT(mm)/dt))=1;       % 1 indicates power on. 
    Qmm(:,:,:)=Q(:,:,:,mm);                 % Power deposited at FZ location mm
    for nn=1:nt                             % Run Model for each timestep at FZ location mm
        cc=c_old;                           % Counter starts at 1 (line 120)
        c_old=cc+1;                         % Counter increments by 1 each iteration
        waitbar(cc/NT,h);                   % Increment the waitbar
        % Shift Temperatures for use in solver
            T2 = circshift(T_new,[1 0 0 0]);
            T3 = circshift(T_new,[-1 0 0 0]);
            T4 = circshift(T_new,[0 1 0 0]);
            T5 = circshift(T_new,[0 -1 0 0]);
            T6 = circshift(T_new,[0 0 1 0]);
            T7 = circshift(T_new,[0 0 -1 0]);
        % Solve for Temperature of Internal Nodes  (TEMPS_new = New Temperature)
            T_new(2:J,2:K,2:L) = squeeze (Coeff1.*...                                                       % Conduction associated with neighboring voxels
                                                    (T2(2:J,2:K,2:L)./(inv_k2k1)+T3(2:J,2:K,2:L)./(inv_k3k1)...       % x direction conduction
                                               +A^2*(T4(2:J,2:K,2:L)./(inv_k4k1)+T5(2:J,2:K,2:L)./(inv_k5k1))...      % y direction conduction
                                               +B^2*(T6(2:J,2:K,2:L)./(inv_k6k1)+T7(2:J,2:K,2:L)./(inv_k7k1)))...     % z direction conduction
                                          +Perf...                                                          % Perfusion associated with difference between baseline and Tb temperature
                                          +Qmm*PowerOn(nn)*dt./rho_cp...                                    % FUS power
                                          +T_old(2:J,2:K,2:L).*Coeff2);                                         % Temperature changes associated with this voxel's old temperature
                                                        
        % Make recently calculated temperature (T_new) the old temperature (T_old) for the next calculation
            T_old(2:J,2:K,2:L)= T_new(2:J,2:K,2:L);
            if BC==1;                           % Adiabatic Boundary
                T_old(1,:,:)=T_new(2,:,:);
                T_old(end,:,:)=T_new(end-1,:,:);
                T_old(:,1,:)=T_new(:,2,:);
                T_old(:,end,:)=T_new(:,end-1,:);
                T_old(:,:,1)=T_new(:,:,2);
                T_old(:,:,end)=T_new(:,:,end-1);
            elseif BC==2;                       % Matching Slope Boundary
                T_old(1,:,:)=T_new(2,:,:)-(T_new(3,:,:)-T_new(2,:,:));
                T_old(end,:,:)=T_new(end-1,:,:)-(T_new(end-2,:,:)-T_new(end-1,:,:));
                T_old(:,1,:)=T_new(:,2,:)-(T_new(:,3,:)-T_new(:,2,:));
                T_old(:,end,:)=T_new(:,end-1,:)-(T_new(:,end-2,:)-T_new(:,end-1,:));
                T_old(:,:,1)=T_new(:,:,2)-(T_new(:,:,3)-T_new(:,:,2));
                T_old(:,:,end)=T_new(:,:,end-1)-(T_new(:,:,end-2)-T_new(:,:,end-1));
            end
            
            if rem(cc,timeratio)==0                 % Determine whether to save the current TEMPS or not
                if use_file
                    fwrite(fid, T_new(2:J,2:K,2:L), 'single');
                else
                    TEMPS(:,:,:,cc/timeratio+1)=T_new(2:J,2:K,2:L);
                end
            end
    end
end

if use_file
    fclose(fid)
end

clear  T_new T_old T2 T3 T4 T5 T6 T7 k1 k2 k3 k4 k5 k6 k7 w_m rho_m cp_m rho_cp lambda
toc   % Stops the stopwatch
close(h);