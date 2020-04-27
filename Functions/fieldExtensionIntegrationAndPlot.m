function integ_value = fieldExtensionIntegrationAndPlot(field_values, outu, outv, numberrepeats, k_u, k_v, domainwidth, domaindepth, slice_axis, ygrid, symmetric_color, show_figures)
field_comp = field_values;
field_comp_temp = field_values(2:end,:);
outu2 = outu;
outv2 = outv;
if numberrepeats > 1
    for m=2:numberrepeats
        field_comp = cat(1,field_comp, field_comp_temp.*exp(-1i*k_u*(m-1)*domainwidth));
        outu2 = cat(1, outu2, outu(2:end,:)+((m-1)*domainwidth));
        outv2 = cat(1, outv2, outv(2:end,:));
    end
end
if slice_axis == Axis.y && ygrid
    field_comp_temp = field_comp(:, 2:end);
    outu2_temp = outu2(:, 2:end);
    outv2_temp = outv2(:, 2:end);
    
    if numberrepeats > 1
        for n=2:numberrepeats
            field_comp = cat(2, field_comp, field_comp_temp);
            outu2 = cat(2, outu2, outu2_temp);
            outv2 = cat(2, outv2, outv2_temp + ((n-1)*domaindepth));
        end
    end
end

integ_value = trapz(outv2(1,:), trapz(outu2(:,1)', field_comp, 1));

if show_figures
    figure
    hold on
    axhand = gca;
    fighand = pcolor(outu2, outv2, real(field_comp));
    daspect(axhand, [1 1 1])
    set(fighand, 'EdgeColor', 'none')
    set(fighand, 'FaceColor', 'interp')
    xlabel('Location (nm)','Interpreter','latex')
    ylabel('Location (nm)','Interpreter','latex')
    colormap(axhand, b2r(4096))
    colorbar
    
    if symmetric_color
        maxabsvalue = max(max(abs(real(field_comp))));
        caxis([-maxabsvalue maxabsvalue]) %For 0 = white (symmetric color axis)
    else
        maxvalue = max(max(real(field_comp)));
        minvalue = min(min(real(field_comp)));
        caxis([minvalue maxvalue])
    end
    
    hold off
end