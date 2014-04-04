function dist_plotting(F,c,dm,p,lim,err_on_point,print_type)

% Copyright (c) Jonathan E. Fieldsend 2013
        
    
figure
hold off

for j=1:max(F)
    scatter(c(F==j,1), c(F==j,2), 20, err_on_point(F==j), 'filled');
    hold on
    for i= 1:max(F)-1
        if print_type ==1
            I=[];
        elseif print_type ==2
            I = find(F(p)==i); % index of shell
        else
            I = find(F(p)<=i); % index of shell
        end
        II = find(F(p)==(i+1)); % index of next shell
        for k=1:length(I)
            for m = 1:length(II)
                if (dm(p(II(m)),p(I(k)))==1)
                    % 'drawing'
                    %domination, so draw line
                    plot([c(p(I(k)),1) c(p(II(m)),1)], [c(p(I(k)),2) c(p(II(m)),2)], 'k:');
                end
            end
        end
    end
end
axis(lim)
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
axis square
colorbar
refresh
drawnow


end