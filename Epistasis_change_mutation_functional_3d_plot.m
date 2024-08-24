% Positive
figure1 = figure;

load('../Epistasis/Epistasis Change/3D_data_Functional_Positive_mutation.mat')

X1 = double(transpose(x));
axes1 = axes('Parent',figure1);
hold(axes1,'on');

bar3(X1);

xlabel('Number of neighbours',...
    'FontWeight','bold',...
    'Rotation',-18,...
    'Position',[8,12,0],...
    'VerticalAlignment','top','HorizontalAlignment','center');

ylabel(['Epistasis nature of';'background genotype'],'FontWeight','bold','Rotation',29, ...
    'FontWeight','bold',...
    'Position',[3.5,26,0],...
    'VerticalAlignment','top','HorizontalAlignment','center');

zlabel('Fraction') ;
ax = gca;
ax.FontSize = 15;
set(gca,'fontweight','bold','fontsize',15) 
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
grid minor;

xlim(axes1,[0.5 6.5]);
ylim(axes1,[0.5 21.5]);
view(axes1,[127.500003927411 38.4122626474681]);
grid(axes1,'on');
hold(axes1,'off');

set(axes1,'FontSize',15,'FontWeight','bold','LineWidth',2,...
    'MinorGridLineWidth',1,'TickLength',[0.025 0.025],...
    'XMinorGrid','on','XTick',1:1:6,...
    'XTickLabel',{'PE','NE','SSE','RSE','OSE','No'}, ...
    'YMinorGrid','on',...
    'YTick',1:5:21,...
    'YTickLabel',{0,5,10,15,20},...
    'ZMinorGrid','on');

saveas(figure1,'../Plots/Epistasis/Epistasis Change to Positive - Functional.png');
close;

% Negative
figure1 = figure;

load('../Epistasis/Epistasis Change/3D_data_Functional_Negative_mutation.mat')

X1 = double(transpose(x));
axes1 = axes('Parent',figure1);
hold(axes1,'on');

bar3(X1);

xlabel('Number of neighbours',...
    'FontWeight','bold',...
    'Rotation',-18,...
    'Position',[8,12,0],...
    'VerticalAlignment','top','HorizontalAlignment','center');

ylabel(['Epistasis nature of';'background genotype'],'FontWeight','bold','Rotation',29, ...
    'FontWeight','bold',...
    'Position',[3.5,26,0],...
    'VerticalAlignment','top','HorizontalAlignment','center');

zlabel('Fraction') ;
ax = gca;
ax.FontSize = 15;
set(gca,'fontweight','bold','fontsize',15) 
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
grid minor;

xlim(axes1,[0.5 6.5]);
ylim(axes1,[0.5 21.5]);
view(axes1,[127.500003927411 38.4122626474681]);
grid(axes1,'on');
hold(axes1,'off');

set(axes1,'FontSize',15,'FontWeight','bold','LineWidth',2,...
    'MinorGridLineWidth',1,'TickLength',[0.025 0.025],...
    'XMinorGrid','on','XTick',1:1:6,...
    'XTickLabel',{'PE','NE','SSE','RSE','OSE','No'}, ...
    'YMinorGrid','on',...
    'YTick',1:5:21,...
    'YTickLabel',{0,5,10,15,20},...
    'ZMinorGrid','on');

saveas(figure1,'../Plots/Epistasis/Epistasis Change to Negative - Functional.png')
close;

% Single Sign
figure1 = figure;

load('../Epistasis/Epistasis Change/3D_data_Functional_Single Sign_mutation.mat')

X1 = double(transpose(x));
axes1 = axes('Parent',figure1);
hold(axes1,'on');

bar3(X1);

xlabel('Number of neighbours',...
    'FontWeight','bold',...
    'Rotation',-18,...
    'Position',[8,12,0],...
    'VerticalAlignment','top','HorizontalAlignment','center');

ylabel(['Epistasis nature of';'background genotype'],'FontWeight','bold','Rotation',29, ...
    'FontWeight','bold',...
    'Position',[3.5,26,0],...
    'VerticalAlignment','top','HorizontalAlignment','center');

zlabel('Fraction') ;
ax = gca;
ax.FontSize = 15;
set(gca,'fontweight','bold','fontsize',15) 
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
grid minor;

xlim(axes1,[0.5 6.5]);
ylim(axes1,[0.5 21.5]);
view(axes1,[127.500003927411 38.4122626474681]);
grid(axes1,'on');
hold(axes1,'off');

set(axes1,'FontSize',15,'FontWeight','bold','LineWidth',2,...
    'MinorGridLineWidth',1,'TickLength',[0.025 0.025],...
    'XMinorGrid','on','XTick',1:1:6,...
    'XTickLabel',{'PE','NE','SSE','RSE','OSE','No'}, ...
    'YMinorGrid','on',...
    'YTick',1:5:21,...
    'YTickLabel',{0,5,10,15,20},...
    'ZMinorGrid','on');

saveas(figure1,'../Plots/Epistasis/Epistasis Change to Single Sign - Functional.png')
close;

% Reciprocal Sign
figure1 = figure;

load('../Epistasis/Epistasis Change/3D_data_Functional_Reciprocal Sign_mutation.mat')

X1 = double(transpose(x));
axes1 = axes('Parent',figure1);
hold(axes1,'on');

bar3(X1);

xlabel('Number of neighbours',...
    'FontWeight','bold',...
    'Rotation',-18,...
    'Position',[8,12,0],...
    'VerticalAlignment','top','HorizontalAlignment','center');

ylabel(['Epistasis nature of';'background genotype'],'FontWeight','bold','Rotation',29, ...
    'FontWeight','bold',...
    'Position',[3.5,26,0],...
    'VerticalAlignment','top','HorizontalAlignment','center');

zlabel('Fraction') ;
ax = gca;
ax.FontSize = 15;
set(gca,'fontweight','bold','fontsize',15) 
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
grid minor;

xlim(axes1,[0.5 6.5]);
ylim(axes1,[0.5 21.5]);
view(axes1,[127.500003927411 38.4122626474681]);
grid(axes1,'on');
hold(axes1,'off');

set(axes1,'FontSize',15,'FontWeight','bold','LineWidth',2,...
    'MinorGridLineWidth',1,'TickLength',[0.025 0.025],...
    'XMinorGrid','on','XTick',1:1:6,...
    'XTickLabel',{'PE','NE','SSE','RSE','OSE','No'}, ...
    'YMinorGrid','on',...
    'YTick',1:5:21,...
    'YTickLabel',{0,5,10,15,20},...
    'ZMinorGrid','on');

saveas(figure1,'../Plots/Epistasis/Epistasis Change to Reciprocal Sign - Functional.png')
close;

% Other Sign
figure1 = figure;

load('../Epistasis/Epistasis Change/3D_data_Functional_Other Sign_mutation.mat')

X1 = double(transpose(x));
axes1 = axes('Parent',figure1);
hold(axes1,'on');

bar3(X1);

xlabel('Number of neighbours',...
    'FontWeight','bold',...
    'Rotation',-18,...
    'Position',[8,12,0],...
    'VerticalAlignment','top','HorizontalAlignment','center');

ylabel(['Epistasis nature of';'background genotype'],'FontWeight','bold','Rotation',29, ...
    'FontWeight','bold',...
    'Position',[3.5,26,0],...
    'VerticalAlignment','top','HorizontalAlignment','center');

zlabel('Fraction') ;
ax = gca;
ax.FontSize = 15;
set(gca,'fontweight','bold','fontsize',15) 
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
grid minor;

xlim(axes1,[0.5 6.5]);
ylim(axes1,[0.5 21.5]);
view(axes1,[127.500003927411 38.4122626474681]);
grid(axes1,'on');
hold(axes1,'off');

set(axes1,'FontSize',15,'FontWeight','bold','LineWidth',2,...
    'MinorGridLineWidth',1,'TickLength',[0.025 0.025],...
    'XMinorGrid','on','XTick',1:1:6,...
    'XTickLabel',{'PE','NE','SSE','RSE','OSE','No'}, ...
    'YMinorGrid','on',...
    'YTick',1:5:21,...
    'YTickLabel',{0,5,10,15,20},...
    'ZMinorGrid','on');

saveas(figure1,'../Plots/Epistasis/Epistasis Change to Other Sign - Functional.png')
close;

% No Epistasis
figure1 = figure;

load('../Epistasis/Epistasis Change/3D_data_Functional_No Epistasis_mutation.mat')

X1 = double(transpose(x));
axes1 = axes('Parent',figure1);
hold(axes1,'on');

bar3(X1);

xlabel('Number of neighbours',...
    'FontWeight','bold',...
    'Rotation',-18,...
    'Position',[8,12,0],...
    'VerticalAlignment','top','HorizontalAlignment','center');

ylabel(['Epistasis nature of';'background genotype'],'FontWeight','bold','Rotation',29, ...
    'FontWeight','bold',...
    'Position',[3.5,26,0],...
    'VerticalAlignment','top','HorizontalAlignment','center');

zlabel('Fraction') ;
ax = gca;
ax.FontSize = 15;
set(gca,'fontweight','bold','fontsize',15) 
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
grid minor;

xlim(axes1,[0.5 6.5]);
ylim(axes1,[0.5 21.5]);
view(axes1,[127.500003927411 38.4122626474681]);
grid(axes1,'on');
hold(axes1,'off');

set(axes1,'FontSize',15,'FontWeight','bold','LineWidth',2,...
    'MinorGridLineWidth',1,'TickLength',[0.025 0.025],...
    'XMinorGrid','on','XTick',1:1:6,...
    'XTickLabel',{'PE','NE','SSE','RSE','OSE','No'}, ...
    'YMinorGrid','on',...
    'YTick',1:5:21,...
    'YTickLabel',{0,5,10,15,20},...
    'ZMinorGrid','on');

saveas(figure1,'../Plots/Epistasis/Epistasis Change to No Epistasis - Functional.png');
close;