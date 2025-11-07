function PlotColor = getPlotColor(MouseType)
    % Function to return specific plot color based on MouseType
    % Input: 
    %   MouseType - A string representing the type of mouse
    % Output:
    %   PlotColor - A string representing the hex color code for the plot

   switch MouseType
        case 'wtTRPV4l'
            PlotColor = '#484F98'; 
        case 'mtTRPV4'
            PlotColor = '#52BE80'; 
        case 'Control'
            PlotColor = '#808080'; % Grey
        otherwise
            error('Unknown MouseType: %s. Please use TRPV1, wtTRPV4, wtTRPV4l, mtTRPV4, or Control.', MouseType);
   end


end
