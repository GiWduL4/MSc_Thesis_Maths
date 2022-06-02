%% Example 6: Gaussian and mean curvature of a Klein bottle

%%

% initialize audi grid
[u,v] = ndgrid(linspace(0,1, 200));
[u,v] = ainit(u,v,2);

% define weight functions
r = 1;
w = 1e0*(u.^2*v.^3);
r = 1;
w1 = u.^2;
w2 = v.^3;

% w = 1;
% w1 = 0;
% w2 = 0;
% w3 = 0;
% w4 = 0;


%define base
xb = u;
yb = v;
zb = 0;
B0 = [xb;yb;zb];

%define ribbon 1
%kappa_1 is identity
p = u;
q = v;

x = p;
y = q;
z = p^2;%q;
R1 = [x;y;z];

%define ribbon 2
%kappa 2
p = 1-v;
q = u;

x = q;
y = 1-p;
z = 0;%q;
R2 = [x;y;z];

% define ABC surface
S = 1/(w+w1+w2)*(w*B0+w1*R1+w2*R2);

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