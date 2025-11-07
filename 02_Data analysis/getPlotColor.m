function PlotColor = getPlotColor(MouseType)
    % Function to return specific plot color based on MouseType
    % Input: 
    %   MouseType - A string representing the type of mouse
    % Output:
    %   PlotColor - A string representing the hex color code for the plot

    % Define the mapping of MouseType to color
%     switch MouseType
%         case 'TRPV1'
%             PlotColor = '#AC316E'; % Dark Pink
%         case 'wtTRPV4'
%             PlotColor = '#4569A8'; % Dark Blue
%         case 'mtTRPV4'
%             PlotColor = '#80BA7F'; % Green
%         case 'Control'
%             PlotColor = '#808080'; % Grey
%         otherwise
%             error('Unknown MouseType: %s. Please use TRPV1, wtTRPV4, mtTRPV4, or Control.', MouseType);
%     end

%    switch MouseType
%         case 'TRPV1'
%             PlotColor = '#F5B1B2'; %  Red
%         case 'wtTRPV4'
%             PlotColor = '#88B9E2'; %  Blue
%         case 'wtTRPV4l'
%             PlotColor = '#484F98'; %  Blue
%         case 'mtTRPV4'
%             PlotColor = '#81C383'; % Green
%         case 'Control'
%             PlotColor = '#A19F9E'; % Grey
%         otherwise
%             error('Unknown MouseType: %s. Please use TRPV1, wtTRPV4, wtTRPV4l, mtTRPV4, or Control.', MouseType);
%    end

   % switch MouseType
   %      case 'TRPV1'
   %          PlotColor = '#D5695D'; %  Red
   %      case 'wtTRPV4'
   %          PlotColor = '#8A7067'; %  Blue
   %      case 'wtTRPV4l'
   %          PlotColor = '#484F98'; %  Blue
   %      case 'mtTRPV4'
   %          PlotColor = '#52BE80'; % Green
   %      case 'Control'
   %          PlotColor = '#808080'; % Grey
   %      otherwise
   %          error('Unknown MouseType: %s. Please use TRPV1, wtTRPV4, wtTRPV4l, mtTRPV4, or Control.', MouseType);
   % end

  switch MouseType
        case 'TRPV1'
            PlotColor = '#D5695D'; %  Red
        case 'wtTRPV4'
            PlotColor = '#417FD2'; %  Blue
        case 'wtTRPV4l'
            PlotColor = '#417FD2'; %  Blue
        case 'mtTRPV4'
            PlotColor = '#EF4156'; % Green
        case 'Control'
            PlotColor = '#808080'; % Grey
        otherwise
            error('Unknown MouseType: %s. Please use TRPV1, wtTRPV4, wtTRPV4l, mtTRPV4, or Control.', MouseType);
   end

end
