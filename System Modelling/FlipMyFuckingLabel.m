function FlipMyFuckingLabel(gca)
    ylh = get(gca,'ylabel');
    gyl = get(ylh);                                                         % Object Information
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
    gca.YLabel.PositionMode='auto';
end