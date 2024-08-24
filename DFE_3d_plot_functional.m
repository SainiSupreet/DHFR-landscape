load('../DFE/3d_data_functional.mat');
xdata1 = x;
ydata1 = y;
zdata1 = z;

% Create figure
figure1 = figure;
colormap(jet);

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.14921253417468 0.794627437641294 0.77578746582532]);
hold(axes1,'on');

% Create mesh
meshz(xdata1,ydata1,zdata1,'Tag','meshz','Parent',axes1);

% Create zlabel
zlabel({'Probability density',''},'FontWeight','bold');

% Create ylabel
ylabel('Selection coefficient','FontWeight','bold','Rotation',-7);

% Create xlabel
xlabel('Background fitness','FontWeight','bold','Rotation',60);

xlim(axes1,[-0.428926404683203 1.22107359531679]);
ylim(axes1,[-2.40730053227758 2.09269946772242]);
zlim(axes1,[1e-06 3.46937812528149]);

view(axes1,[-251.996834978082 60.9607368671215]);
grid(axes1,'on');
hold(axes1,'off');

set(axes1,'Colormap',...
    [1 1 1;0.680555555555556 0.666666666666667 0.666666666666667;0.361111111111111 0.333333333333333 0.333333333333333;0.0416666666666666 0 0;0.0520833333333333 0 0;0.0625 0 0;0.0729166666666667 0 0;0.0833333333333333 0 0;0.09375 0 0;0.104166666666667 0 0;0.114583333333333 0 0;0.125 0 0;0.135416666666667 0 0;0.145833333333333 0 0;0.15625 0 0;0.166666666666667 0 0;0.177083333333333 0 0;0.1875 0 0;0.197916666666667 0 0;0.208333333333333 0 0;0.21875 0 0;0.229166666666667 0 0;0.239583333333333 0 0;0.25 0 0;0.260416666666667 0 0;0.270833333333333 0 0;0.28125 0 0;0.291666666666667 0 0;0.302083333333333 0 0;0.3125 0 0;0.322916666666667 0 0;0.333333333333333 0 0;0.34375 0 0;0.354166666666667 0 0;0.364583333333333 0 0;0.375 0 0;0.385416666666667 0 0;0.395833333333333 0 0;0.40625 0 0;0.416666666666667 0 0;0.427083333333333 0 0;0.4375 0 0;0.447916666666667 0 0;0.458333333333333 0 0;0.46875 0 0;0.479166666666667 0 0;0.489583333333333 0 0;0.5 0 0;0.510416666666667 0 0;0.520833333333333 0 0;0.53125 0 0;0.541666666666667 0 0;0.552083333333333 0 0;0.5625 0 0;0.572916666666667 0 0;0.583333333333333 0 0;0.59375 0 0;0.604166666666667 0 0;0.614583333333333 0 0;0.625 0 0;0.635416666666667 0 0;0.645833333333333 0 0;0.65625 0 0;0.666666666666667 0 0;0.677083333333333 0 0;0.6875 0 0;0.697916666666667 0 0;0.708333333333333 0 0;0.71875 0 0;0.729166666666667 0 0;0.739583333333333 0 0;0.75 0 0;0.760416666666667 0 0;0.770833333333333 0 0;0.78125 0 0;0.791666666666667 0 0;0.802083333333333 0 0;0.8125 0 0;0.822916666666667 0 0;0.833333333333333 0 0;0.84375 0 0;0.854166666666667 0 0;0.864583333333333 0 0;0.875 0 0;0.885416666666667 0 0;0.895833333333333 0 0;0.90625 0 0;0.916666666666667 0 0;0.927083333333333 0 0;0.9375 0 0;0.947916666666667 0 0;0.958333333333333 0 0;0.96875 0 0;0.979166666666667 0 0;0.989583333333333 0 0;1 0 0;1 0.0104166666666667 0;1 0.0208333333333333 0;1 0.03125 0;1 0.0416666666666667 0;1 0.0520833333333333 0;1 0.0625 0;1 0.0729166666666667 0;1 0.0833333333333333 0;1 0.09375 0;1 0.104166666666667 0;1 0.114583333333333 0;1 0.125 0;1 0.135416666666667 0;1 0.145833333333333 0;1 0.15625 0;1 0.166666666666667 0;1 0.177083333333333 0;1 0.1875 0;1 0.197916666666667 0;1 0.208333333333333 0;1 0.21875 0;1 0.229166666666667 0;1 0.239583333333333 0;1 0.25 0;1 0.260416666666667 0;1 0.270833333333333 0;1 0.28125 0;1 0.291666666666667 0;1 0.302083333333333 0;1 0.3125 0;1 0.322916666666667 0;1 0.333333333333333 0;1 0.34375 0;1 0.354166666666667 0;1 0.364583333333333 0;1 0.375 0;1 0.385416666666667 0;1 0.395833333333333 0;1 0.40625 0;1 0.416666666666667 0;1 0.427083333333333 0;1 0.4375 0;1 0.447916666666667 0;1 0.458333333333333 0;1 0.46875 0;1 0.479166666666667 0;1 0.489583333333333 0;1 0.5 0;1 0.510416666666667 0;1 0.520833333333333 0;1 0.53125 0;1 0.541666666666667 0;1 0.552083333333333 0;1 0.5625 0;1 0.572916666666667 0;1 0.583333333333333 0;1 0.59375 0;1 0.604166666666667 0;1 0.614583333333333 0;1 0.625 0;1 0.635416666666667 0;1 0.645833333333333 0;1 0.65625 0;1 0.666666666666667 0;1 0.677083333333333 0;1 0.6875 0;1 0.697916666666667 0;1 0.708333333333333 0;1 0.71875 0;1 0.729166666666667 0;1 0.739583333333333 0;1 0.75 0;1 0.760416666666667 0;1 0.770833333333333 0;1 0.78125 0;1 0.791666666666667 0;1 0.802083333333333 0;1 0.8125 0;1 0.822916666666667 0;1 0.833333333333333 0;1 0.84375 0;1 0.854166666666667 0;1 0.864583333333333 0;1 0.875 0;1 0.885416666666667 0;1 0.895833333333333 0;1 0.90625 0;1 0.916666666666667 0;1 0.927083333333333 0;1 0.9375 0;1 0.947916666666667 0;1 0.958333333333333 0;1 0.96875 0;1 0.979166666666667 0;1 0.989583333333333 0;1 1 0;1 1 0.015625;1 1 0.03125;1 1 0.046875;1 1 0.0625;1 1 0.078125;1 1 0.09375;1 1 0.109375;1 1 0.125;1 1 0.140625;1 1 0.15625;1 1 0.171875;1 1 0.1875;1 1 0.203125;1 1 0.21875;1 1 0.234375;1 1 0.25;1 1 0.265625;1 1 0.28125;1 1 0.296875;1 1 0.3125;1 1 0.328125;1 1 0.34375;1 1 0.359375;1 1 0.375;1 1 0.390625;1 1 0.40625;1 1 0.421875;1 1 0.4375;1 1 0.453125;1 1 0.46875;1 1 0.484375;1 1 0.5;1 1 0.515625;1 1 0.53125;1 1 0.546875;1 1 0.5625;1 1 0.578125;1 1 0.59375;1 1 0.609375;1 1 0.625;1 1 0.640625;1 1 0.65625;1 1 0.671875;1 1 0.6875;1 1 0.703125;1 1 0.71875;1 1 0.734375;1 1 0.75;1 1 0.765625;1 1 0.78125;1 1 0.796875;1 1 0.8125;1 1 0.828125;1 1 0.84375;1 1 0.859375;1 1 0.875;1 1 0.890625;1 1 0.90625;1 1 0.921875;1 1 0.9375;1 1 0.953125;1 1 0.96875;1 1 0.984375;1 1 1],...
    'FontSize',15,'FontWeight','bold','GridColor',...
    [0.149019607843137 0.149019607843137 0.149019607843137],'LineWidth',2,...
    'TickLength',[0.025 0.025],'XMinorGrid','on','YMinorGrid','on','ZMinorGrid',...
    'on');