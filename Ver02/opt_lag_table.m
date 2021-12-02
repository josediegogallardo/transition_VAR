%  The following routine writes a table with the results of Several
%  Non-Linearity Tests:

%%%%%%%%%%%%%%%%%%%%%%% Likelihood Ratio Portion %%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating Report
h=figure('units','normalize','Name','Report Table','Position',[0.0 0.05 1 1]);
axes('Position',[0.05 0.2 0.9 0.75],'Color','w');
x=0.5;
y=1;
% Title and Subtitles:
text(x,y,['\bfResults for the Likelihood Ratio Test for Non-Linearities'],'Horizontalalignment','center','Fontsize',14,'Color',[0.8 0 0.05]);
y=y-0.05;
text(x,y,['Likelihood Ratio (P-Value)'],'Horizontalalignment','center','Fontsize',12,'Color',[0 0 0],'FontWeight','bold');
% line([0.05 0.95],[0.85 0.85],'LineWidth',1,'Color','k'); 

% Adding Some Picky Lines
% Variable Names
y=0.825;
x=0.12;
text(x,y,['Variable'],'Horizontalalignment','center','Fontsize',11,'Color',[0.9 0 0.2],'FontWeight','bold');
step=(1-0.3)/neqs;
xs=[] ;
for ii=1:neqs
    x=x+step;
    xs=[xs x];
    text(x,y,[names(ii).name],'Horizontalalignment','center','Fontsize',11,'Color',[0.9 0 0.2],'FontWeight','bold');
end

% Model Names
x=0.12;
for ii=1:transilag
    y=y-0.05;
    text(x,y,['Lag ' num2str(ii)],'Horizontalalignment','center','Fontsize',11,'Color',[0.9 0 0.2],'FontWeight','bold');
end

% filling the variables
for jj=1:neqs
   y=0.825;
   x=xs(jj);
   for ii=1:transilag
       y=y-0.05;
       text(x,y,[num2str(round(LRtests(ii,jj)*100)/100)],'Horizontalalignment','center','Fontsize',11,'Color','k','FontWeight','bold'); 
   end
end;

axis off;

orient landscape;
print -dpsc2  -loose -r300 -f 'Likelihood Ratio Tests';

%%%%%%%%%%%%%%%%%%%%%% LM Test Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for qq=1:neqs    
    h=figure('units','normalize','Name','Report Table','Position',[0.0 0.05 1 1]);
    axes('Position',[0.05 0.2 0.9 0.75],'Color','w');
    x=0.5;
    y=1;
    % Title and Subtitles:
    text(x,y,['\bfResults for the Lagrange Multiplier Test for Non-Linearities'],'Horizontalalignment','center','Fontsize',14,'Color',[0.8 0 0.05]);
    y=y-0.05;
    text(x,y,[names(qq).name ' Equation (P-Value)'],'Horizontalalignment','center','Fontsize',12,'Color',[0 0 0],'FontWeight','bold');
    % line([0.05 0.95],[0.85 0.85],'LineWidth',1,'Color','k'); 

    % Adding Some Picky Lines
    % Variable Names
    y=0.825;
    x=0.12;
    text(x,y,['Variable'],'Horizontalalignment','center','Fontsize',11,'Color',[0.9 0 0.2],'FontWeight','bold');
    step=(1-0.3)/neqs
    xs=[];
    for ii=1:neqs
        x=x+step;
        xs=[xs x];
        text(x,y,[names(ii).name],'Horizontalalignment','center','Fontsize',11,'Color',[0.9 0 0.2],'FontWeight','bold');
    end

    % Model Names
    x=0.12;
    for ii=1:transilag
        y=y-0.05;
        text(x,y,['Lag ' num2str(ii)],'Horizontalalignment','center','Fontsize',11,'Color',[0.9 0 0.2],'FontWeight','bold');
    end

    % filling the variables

    for jj=1:neqs
       y=0.825 ; 
       x=xs(jj);
       for ii=1:transilag
           y=y-0.05;
           text(x,y,[num2str(round(fstat(ii,jj,qq)*100)/100)],'Horizontalalignment','center','Fontsize',11,'Color','k','FontWeight','bold'); 
       end
    end;

    axis off;

    orient landscape;
    eval(['print -dpsc2  -loose -r300 -f ''' names(qq).name ' Equation (P-Value)'''])
end