% comp2DDM.m
% This script couples a complementarity algorithm and the two dimensional
%   displacement discontinuity method (DDM) to solve boundary value 
%   problems for isotropic, homogeneous, linear elastic solids containing 
%   fractures with user-defined frictional surface properties.
% --> Complementarity takes the form: f(Z)=M*Z+Q where f(Z)*Z=0 (e.g. 
%   Pang and Trinkle, 1996, Complementarity formulations and existence of 
%   solutions of dynamic multi-rigid-body contact problems with coulomb 
%   friction, Mathematical Programming, Springer).
% Complementarity scripts (for example, lemke.m, scoop.m, or pathlcp.m) 
%   and any other necessary files must reside in current directory.
% --> DDM is from Chap. 5 of Crouch and Starfield, 1983, Boundary Element 
%   Methods in Solid Mechanics, Unwin Hyman, London (abbreviated C&S).
% The number of elements making up each fracture and the element tip 
%   coordinates are stored in an Excel file comp2DDM_test.xls.
% Each fracture is digitized w/re to the (x,y)-global coordinate system
%   such that Ne(i)+1 points define the tips of Ne(i) straight boundary 
%   elements making up fracture i. 
% Code written by Ovunc (Uno) Mutlu (2006).
% Code modified by David D. Pollard (2007) & Elizabeth (Libby) Ritz (2011).

%% Clear variables, globals, functions, close figures
clear all; % Clear variables, globals, functions from memory.
close all; % Close all open figures
clc; % Clear the command window

%% Input parameters from Excel worksheet
% Read row vector with input parameters from Excel worksheet.
Parm = xlsread('comp2DDM_test.xls','parm');
% Elastic constants: Poisson's ratio, pr, and Young's modulus, ym.
% Note: units for ym and stress components are consistent (MPa).
pr = Parm(1);  ym = Parm(2);
% Three remote stresses in (x,y)-global coordinate system.
% Tension is taken as positive.
sxxr = Parm(3); syyr = Parm(4); sxyr = Parm(5);
% Number of fractures.
nf = Parm(6);
% Define the observation grid.
[XG,YG] = meshgrid(linspace(Parm(7),Parm(9),Parm(8)),...
    linspace(Parm(10),Parm(12),Parm(11))); 
% Read row vector with numbers of elements for each fracture.
Ne = xlsread('comp2DDM_test.xls','no');
ne = sum(Ne); % Total number of elements, including all fractures.
% Read rectangular arrays of the x- and y-coordinates of the digitized
%   points along all fractures in the global coordinate system.
PX = xlsread('comp2DDM_test.xls','px');
PY = xlsread('comp2DDM_test.xls','py');
% Read rectangular array of friction coefficients for all elements.
MU = xlsread('comp2DDM_test.xls','mu');
% Parse MU into column vector of length ne, omitting nans.
Mu = MU(isfinite(MU));
% Read rectangular array(s) of frictional strength, Sf, coefficients.
SF = xlsread('comp2DDM_test.xls','sf');
% Parse into column vector of length ne, omitting nans.
Sf = SF(isfinite(SF));
if ne == 2 % For 2 elements transpose Mu and Sf to column vectors.
    Mu = Mu'; Sf = Sf';
end
% Plot observation grid and fractures
plot(XG,YG,'.k'), hold on
for i = 1:nf
    plot(PX(:,i),PY(:,i))
end
axis equal, hold off
title('observation grid and fractures'), xlabel('x'), ylabel('y')

%% DDM formulation (from C&S chapter 5, p. 79-93)
% Parse coordinate data on element tips from padded arrays into single
%   column vectors of beginning and ending coordinates of all elements
%   in the (x,y)-global system.
% Note: the arrangement of elements into fractures is not relevant to the
%   procdure for solving the boundary value problem in which each
%   element is treated independently and simultaneously.
Xbeg = PX(1:Ne(1),1); Ybeg = PY(1:Ne(1),1);
Xend = PX(2:Ne(1)+1,1); Yend = PY(2:Ne(1)+1,1);
for i = 2:nf
    Xbeg = [Xbeg ; PX(1:Ne(i),i)]; Ybeg = [Ybeg ; PY(1:Ne(i),i)];
    Xend = [Xend ; PX(2:Ne(i)+1,i)]; Yend = [Yend ; PY(2:Ne(i)+1,i)];
end

% In what follows several square arrays are defined such that: 
%   the row number varies from i=1 to i=ne, and 
%   the column number varies from j=1 to j=ne.
% Define column vectors of projected element lengths onto the 
%   (x,y)-global coordinate axes.
Xpro = Xend-Xbeg; Ypro = Yend-Ybeg;
% Define column vector of actual half lengths of each element.
H = 0.5*sqrt(Xpro.*Xpro + Ypro.*Ypro);
% Define square array where each row contains all element half lengths
%   ordered from j=1 to j=ne.
HJ = repmat(H',ne,1);
% Define column vector of element midpoint coordinates in the
%   global (x,y)-system.				
Xmid = (Xbeg + Xend)/2; Ymid = (Ybeg + Yend)/2;
% Define square array where each column contains all midpoint coordinates.
XI=repmat(Xmid,1,ne); YI=repmat(Ymid,1,ne);
% Define square array where each row contains all midpoint coordinates.
XJ=repmat(Xmid',ne,1); YJ=repmat(Ymid',ne,1);

% Calculate orientations of elements in (x,y)-global coordinate system.
% Beta is (smaller) ccw angle (in radians) between Ox and element trace.
% This is the convention of C&S, Figures 4.7, p. 60 and 5.6, p. 91.
Beta = atan2(Ypro,Xpro);
% Define square array where each row contains all the angles Beta.
BJ = repmat(Beta',ne,1);
% Define square array of local coordinates for the ith element midpoint
%   with respect to an origin at the jth element midpoint and 
%   the Ox axis aligned with the jth element (Figure 4.7, C&S, p. 60). 
% These are eqs. 4.6.2a and b of C&S, p. 60.
XB = (XI-XJ).*cos(BJ) + (YI-YJ).*sin(BJ);
YB = -(XI-XJ).*sin(BJ) + (YI-YJ).*cos(BJ);

% Calculate derivatives of the function f(x,y), eq. 5.2.5 of C&S, p. 81, 
%   which are used to calculate the stress influence coefficients. 
% The following functions are used in these calculations.
% It is understood that X and Y refer to XB and YB.
Y2 = YB.^2;
XMH = XB-HJ; XPH = XB+HJ;
XMH2 = XMH.^2; XPH2 = XPH.^2; 
R1S = XMH2 + Y2; R1S2 = R1S.^2;
R2S = XPH2 + Y2; R2S2 = R2S.^2;
% The shear modulus, sm, is related to the prescribed elastic constants.
sm = ym/(2*(1+pr));
% Define material constant used in calculating influence coefficients.
con = 1/(4*pi*(1-pr)); 
% The following derivatives are eqs. 4.5.5c and d of C&S, p. 58.
FF4 = con*(YB./R1S - YB./R2S); 
FF5 = con*(XMH./R1S - XPH./R2S);
% The following derivatives are eqs. 5.5.3a and b of C&S, p. 91.
FF6 = con*((XMH2 - Y2)./R1S2 - (XPH2 - Y2)./R2S2);
FF7 = 2*con*YB.*(XMH./R1S2 - XPH./R2S2);

% Define square array with inclinations of ith element with respect to jth.
% Gamma is ccw angle between local Ox for jth and Ox for ith elements.
Gamma = repmat(Beta,1,ne) - BJ;
% Calculate trigonometric functions of the angle Gamma for each element, 
%   which are also used to calculate the stress influence coefficients. 
SG = sin(Gamma); CG = cos(Gamma); SG2 = SG.^2;
S2G = 2*SG.*CG; C2G = CG.^2 - SG2;
% Calculate the shear and normal stress influence coefficients for  
%   the midpoint of the ith element with respect to the orientation of 
%   the ith element plane.
% These are the stress components caused by unit shear and unit normal  
%   displacement discontinuities of the jth element.
% Sign conventions for element displacement discontinuity and stress
%   are from C&S Figure 5.4, p. 87.
% The stress influence coefficients are from eqs. 5.6.2 of C&S, p. 93.
% The element self-effects are given by eqs. 5.6.4 of C&S, p. 93.
Ass = 2*sm*(-S2G.*FF4-C2G.*FF5-YB.*(S2G.*FF6-C2G.*FF7));
Asn = 2*sm*(-YB.*(C2G.*FF6+S2G.*FF7));
Ans = 2*sm*(2*SG2.*FF4+S2G.*FF5-YB.*(C2G.*FF6+S2G.*FF7));
Ann = 2*sm*(-FF5+YB.*(S2G.*FF6-C2G.*FF7));
% Regroup the four stress influence coefficient arrays to form 
%   a square array of size 2*ne by 2*ne such that the shear and  
%   normal stress components, Ss and Sn, on the ith element are related
%   to the shear and normal displacement discontinuities, Ds and Dn,
%   on the jth element as given by eqs 5.4.3 of C&S, p. 88.
%       	  |Ss| =  |Ass Asn||Ds|	 
%             |Sn|    |Ans Ann||Dn|	  
A = [Ass, Asn; Ans, Ann]; % Aray of stress influence coefficients
   
% The 3 remote stress components are resolved into shear stress, Ssr, 
%   and normal stress, Snr, acting on the plane of each element.
% These components are referred to the element local coordinates (s,n)
%   where n is directed at angle beta + pi/2 to the global Ox axis and
%   s is directed to the right w/re to n (Figure 4.7, C&S, p. 60).
% Calculate trigonometric functions of the angle beta for each element.
Sinb = Ypro./(2*H); Cosb = Xpro./(2*H);	% sine and cosine beta
Sbcb = Sinb.*Cosb; % sine beta times cosine beta
Sinb2 = Sinb.^2; Cosb2 = Cosb.^2;	% sine and cosine beta squared
% The transformation equations are 2.8.8c and 2.8.8b from C&S, p. 24.
Ssr = -(sxxr-syyr).*Sbcb + sxyr.*(Cosb2-Sinb2);
Snr = sxxr.*Sinb2 - 2*sxyr.*Sbcb + syyr.*Cosb2;
% Regroup resolved remote stress components into a 2*ne column vector.
Sr=[Ssr;Snr];
    
% Given the stress components, S, and the influence coefficients, A, 
%   we have a system of 2*ne simultaneous linear equations with
%   2*ne unknowns, which are the displacement discontinuities, D.
%       	  |D| =  |C||S|, where C = inv(A)
% Invert the array of influence coefficients, A, to find the array C.
C = A\eye(size(A)); 
% Solve for the column vector of displacement discontinuities, D.
D = C*Sr;
% Signs of D appear to be opposite of S&C convention!!

%% Complementarity formulation
% Reformulate as a complementarity problem.
% Construct M and Q for the equation f(Z)=M*Z+Q, where
%  f(Z) = |-Dn|  and  Z = |-Sn|   --> DR,DL are the magnitudes of slip in
%         | DR|           | SR|   each direction, and SR,SL are the slack
%         | SL|           | DL|   variables for each direction.
% Rename each submatrix for easier assmbly into M.
CSS = C(1:ne,1:ne);
CSN = C(1:ne,ne+1:2*ne);
CNS = C(ne+1:2*ne,1:ne);
CNN = C(ne+1:2*ne,ne+1:2*ne);
% Form ne by ne array with coefficients of friction on the diagonal.
dMU = diag(Mu);
% Allocate ne by ne identity and zero matrices.
ID = eye(ne);
ZE = zeros(ne);
% Construct 3ne by 3ne matrix M. 
M = [(CNN-CNS*dMU)  CNS  ZE;
     (CSN-CSS*dMU)  CSS  ID;
      2*dMU         -ID  ZE];
% Construct 3ne by 1 column vector Q.  
Q = [D(ne+1:2*ne)-CNS*Sf;
     D(1:ne)-CSS*Sf;
     2*Sf]; 

% Solve the complementarity problem using the PATH algorithm.
% Source: http://pages.cs.wisc.edu/~ferris/path/
% pathlcp.m must be in the same Matlab directory.
Z=pathlcp(M,Q);

%% Calculate Ss, Sn, Ds, Dn on each element.
% Calculate complementary equation.
fZ = M*Z+Q;
% Re-calculate vectors Ss, Sn, Ds, Dn from fZ and Z.
Dn = -fZ(1:ne,1); % =-(-Dn)
Ds = Z(2*ne+1:3*ne,1)-fZ(ne+1:2*ne,1); % =DL-DR
Sn = -Z(1:ne,1); % =-(-Sn) 
Ss = -Z(ne+1:2*ne,1)-Mu.*Sn+Sf; % =-SR-Mu*Sn+Sf

% Plot displacement discontinuities on each fracture.
figure, plot(Xmid,Dn,'b.',Xmid,Ds,'r.');
title('displacement discontinuity'), xlabel('x (m)')
ylabel('Dn (blue) and Ds (red) (m)')

% Plot displacement discontinuities on each element.
figure, plot(1:ne,Dn,'b.',1:ne,Ds,'r.')
title('displacement discontinuity'), xlabel('element number')
ylabel('Dn (blue) and Ds (red) (m)')

%% Compute displacement and stress components on the observation grid.
% Initialize the displacement and stress components on the grid.
UX = zeros(size(XG)); UY = zeros(size(XG));
SXX = sxxr*ones(size(XG)); SYY = syyr*ones(size(XG)); 
SXY = sxyr*ones(size(XG));
% Loop over all elements, superimposing their contributions to each field.
for i = 1:ne
    Dxb = Ds(i); Dyb = Dn(i);
    sb = sin(Beta(i)); cb = cos(Beta(i));
    s2b = sin(2*Beta(i)); c2b = cos(2*Beta(i));
% Define array of local coordinates for the observation grid relative to
%   the midpoint and orientation of the ith element.
% Refer to (Figure 5.6, C&S, p. 91) and eqs. 4.5.1 of C&S, p. 57. 
    XB = (XG-Xmid(i))*cb + (YG-Ymid(i))*sb;
    YB = -(XG-Xmid(i))*sb + (YG-Ymid(i))*cb;

% Calculate derivatives of the function f(x,y), eq. 5.2.5 of C&S, p. 81. 
%   which are used to calculate the displacement and stress components. 
% It is understood that X and Y refer to XB and YB.
% First abbreviate repeated terms in the derivatives of f(x,y):
    Y2 = YB.^2;
    XMH = XB-H(i); XPH = XB+H(i); 
    XMH2 = XMH.^2; XPH2 = XPH.^2; 
    R1S = XMH2 + Y2; R1S2 = R1S.^2; 
    R2S = XPH2 + Y2; R2S2 = R2S.^2;
% The following derivatives are eqs. 4.5.5a thru d of C&S, p. 58.
    FF2 = con*(log(sqrt(R1S)) - log(sqrt(R2S)));
    FF3 = -con*(atan2(YB,XMH) - atan2(YB,XPH));
    FF4 = con*(YB./R1S - YB./R2S); 
    FF5 = con*(XMH./R1S - XPH./R2S);
% The following derivatives are eqs. 5.5.3a and b of C&S, p. 91.
    FF6 = con*((XMH2 - Y2)./R1S2 - (XPH2 - Y2)./R2S2);
    FF7 = 2*con*YB.*(XMH./R1S2 - XPH./R2S2);
    
% Define material constants used in calculating displacements.
    pr1 = 1-2*pr; pr2 = 2-2*pr;
% Calculate the displacement components using eqs. 5.5.4 of C&S, p. 91.
    UX = UX + Dxb*(-pr1*sb*FF2 + pr2*cb*FF3 + YB.*(sb*FF4 - cb*FF5))...
        +Dyb*(-pr1*cb*FF2 - pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5));
    UY = UY + Dxb*(+pr1*cb*FF2 + pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5))...
        +Dyb*(-pr1*sb*FF2 + pr2*cb*FF3 - YB.*(sb*FF4 - cb*FF5));

% Calculate the stress components using eqs. 5.5.5 of C&S, p. 92.
    SXX = SXX + 2*sm*Dxb*(2*cb*cb*FF4 + s2b*FF5 + YB.*(c2b*FF6-s2b*FF7))...
         +2*sm*Dyb*(-FF5 + YB.*(s2b*FF6 + c2b*FF7));
    SYY = SYY + 2*sm*Dxb*(2*sb*sb*FF4 - s2b*FF5 - YB.*(c2b*FF6-s2b*FF7))...
         +2*sm*Dyb*(-FF5 - YB.*(s2b*FF6 + c2b*FF7));
    SXY = SXY + 2*sm*Dxb*(s2b*FF4 - c2b*FF5 + YB.*(s2b*FF6+c2b*FF7))...
         +2*sm*Dyb*(-YB.*(c2b*FF6 - s2b*FF7));    
end

% %% Additional calculations
% % Calculate principal stresses, S1 and S2.
% % Equations from Pollard and Fletcher, 2005, Fundamentals of Structural 
% %   Geology, Cambridge University Press (abbreviated P&F)
% % P&F  page 220, equations 6.72.
% S1 = 0.5*(SXX+SYY)+sqrt(0.25*(SXX-SYY).^2+SXY.^2); 
% S2 = 0.5*(SXX+SYY)-sqrt(0.25*(SXX-SYY).^2+SXY.^2); 
% 
% % Calculate the maximum shear stress and mean normal stress.
% % P&F page 223, table 6.1
% SS = 0.5*(S1-S2); 
% % P&F page 223, table 6.1
% MS = 0.5*(S1+S2); 

%% Figures
% Plot vector arrows for displacements.
figure, hold on, axis equal tight
quiver(XG,YG,UX,UY, 'linewidth', 1.2), plot(PX,PY,'k', 'linewidth', 2)
axis([Parm(7), Parm(9), Parm(10), Parm(12)])
title('displacement field'), xlabel('x/w'), ylabel('y/w')

% Plot the magnitude of the Cartesian stress components.
figure, hold on, axis equal tight
contourf(XG,YG,SXX,18), plot(PX,PY,'k', 'linewidth', 2)
axis([Parm(7), Parm(9), Parm(10), Parm(12)])
title('SXX (MPa)'), xlabel('x (m)'), ylabel('y (m)'), colorbar

figure, hold on, axis equal tight
contourf(XG,YG,SYY,18), plot(PX,PY,'k', 'linewidth', 2)
axis([Parm(7), Parm(9), Parm(10), Parm(12)])
title('SYY (MPa)'), xlabel('x (m)'), ylabel('y (m)'), colorbar

figure, hold on, axis equal tight
contourf(XG,YG,SXY,18), plot(PX,PY,'k', 'linewidth', 2)
axis([Parm(7), Parm(9), Parm(10), Parm(12)])
title('SXY (MPa)'), xlabel('x (m)'), ylabel('y (m)'), colorbar

% % Plot the magnitude of the principal stresses.
% figure, hold on, axis equal tight
% contourf(XG,YG,S1,18), plot(PX,PY,'k', 'linewidth', 2)
% axis([Parm(7), Parm(9), Parm(10), Parm(12)])
% title('S1 (MPa)'), xlabel('x (m)'), ylabel('y (m)')
% 
% figure, hold on, axis equal tight
% contourf(XG,YG,S2,18), plot(PX,PY,'k', 'linewidth', 2)
% axis([Parm(7), Parm(9), Parm(10), Parm(12)])
% title('S2 (MPa)'), xlabel('x (m)'), ylabel('y (m)')
% 
% % Plot the maximum shear and mean normal stresses.
% figure, hold on, axis equal tight
% contourf(XG,YG,SS,18), plot(PX,PY,'k', 'linewidth', 2)
% axis([Parm(7), Parm(9), Parm(10), Parm(12)])
% title('maximum shear stress (MPa)'), xlabel('x (m)'), ylabel('y (m)')
% 
% figure, hold on, axis equal tight
% contourf(XG,YG,MS,18), plot(PX,PY,'k', 'linewidth', 2)
% axis([Parm(7), Parm(9), Parm(10), Parm(12)])
% title('mean normal stress (MPa)'), xlabel('x (m)'), ylabel('y (m)')