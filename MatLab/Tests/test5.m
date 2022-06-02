%% Example 6: Gaussian and mean curvature of a Klein bottle

%%

% initialize audi grid
[u,v] = ndgrid(linspace(0,1, 200));
[u,v] = ainit(u,v,2);

% define weight functions
r = 1;
w = 1e1*(u.^r*v.^r*(1-u).^r*(1-v)^r);
r = 1;
w1 = u.^r*(1-u).^r*(1-v)^r;
w2 = u.^r*v.^r*(1-v).^r;
w3 = u.^r*v.^r*(1-u)^r;
w4 = v.^r*(1-u).^r*(1-v)^r;

% w = 1;
% w1 = 0;
% w2 = 0;
% w3 = 0;
% w4 = 0;


%define base
xb = u;
yb = v;
zb = 3;
B0 = [xb;yb;zb];

%define ribbon 1
%kappa_1 is identity
p = u;
q = v;

x = p;
y = q;
z = 0;%q;
R1 = [x;y;z];

%define ribbon 2
%kappa 2
p = v;
q = 1-u;

x = (1-q);
y = p;
z = 0;%q;
R2 = [x;y;z];

%define ribbon 3
%kappa 3
p = 1-u;
q = 1-v;

x = (1-p);
y = 1-q;
z = 0;%q;
R3 = [x;y;z];

%define ribbon 4
%kappa 4
p = 1-v;
q = u;

x = q;
y = 1-p;
z = 0;%q;
R4 = [x;y;z];

% define ABC surface
S = 1/(w+w1+w2+w3+w4)*(w*B0+w1*R1+w2*R2+w3*R3+w4*R4);
% x = u;
% y = v;
% z = 0*u; %u.*v+u.*(1-u-v)+v.*(1-u-v);

% S = [x;y;z];

x = S(1,:);
y = S(2,:);
z = S(3,:);

% differential geometry
J = ajac(S);                              % Jacobian
N = cross(J(:,1),J(:,2));                 % normal vector
N = N/norm(N);                            % normalize
G = J'*J;                                 % first fundamental form
B = [N'*adiff(S,2,0) N'*adiff(S,1,1);...  % second fundamental form
     N'*adiff(S,1,1) N'*adiff(S,0,2)];
W = G\B;                                  % Weingarten map
K = det(W);                               % Gaussian curvature
H = trace(W)/2;                           % mean curvature

% color surface by Gaussian curvature
figure(1), clf, colormap(jet)
surf(x{0},y{0},z{0},K{0})
light, light, shading interp
axis tight, caxis([-1.5 .5]);
title('Gaussian curvature')

% color surface by mean curvature
% figure(2), clf, colormap(jet)
% surf(x{0},y{0},z{0},H{0})
% light, light, shading interp
% axis equal, caxis([-1 1])
% title('Mean curvature, discontinuous since surface is not orientable')