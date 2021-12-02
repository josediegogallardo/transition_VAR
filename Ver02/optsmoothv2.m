function results = optsmoothv2(Z, param, step_smooth, points_smooth, prop_dist, points_thresh)

% PURPOSE: Performs a grid search in order to find optimal threshold and 

mini=min(Z(:,param(2)))                                     ;
maxi=max(Z(:,param(2)))                                     ;
distance = (maxi-mini)                                      ;
start_thresh= mean(Z(:,param(2))) - (prop_dist/2)*distance  ;
step_thresh=(distance*prop_dist)/points_thresh              ;

for j = 1:points_smooth
                smooth_line(j)=1+(j-1)*step_smooth          ;
end;    
%good
for i = 1:points_thresh
            param(4)=start_thresh+i*step_thresh             ;
            threshold_line(i)=start_thresh+(i-1)*step_thresh;
            for j = 1:points_smooth
                param(5)=1+(j-1)*step_smooth                ;
                
                if nargin == 7 
                    stvar1 = stvar_tests_gfunc(y, param, x) ;
                else
                    stvar1 = stvar_tests_gfunc(y, param)    ;
                end

            temppval=LRpval(j,i)                    ;
            smooth=j*step_smooth                    ;
            threshold=start_thresh+i*step_thresh    ;
            end;
end;

results.smooth=smooth                       ;
results.smooth_line=smooth_line             ;
results.threshold=threshold                 ;
results.threshold_line=threshold_line       ;
results.LRpval=LRpval                       ;
results.pval=temppval                       ;
results.start_thresh=start_thresh           ;
results.min=mini                            ;
results.max=maxi                            ;

surf(threshold_line,smooth_line,LRpval)     ;
axis tight;
ylabel('smoothness parameter (\lambda)')    ;
xlabel('threshold value (\theta)')          ;    
zlabel('LR P-Value')                        ;

orient landscape;
print -dpsc2  -loose -r300 -f 'Parameter Grid Search'