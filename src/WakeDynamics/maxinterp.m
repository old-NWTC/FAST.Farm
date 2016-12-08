% Line search


xCoords = [10, 3, 4, 5, 3, 2, 4, 5, 7, 4, 5];
yCoords = [1, 2, 3, 4, 5,  6, 7, 8,  9, 10, 11];
yPt = 0.0;
xPt = 3.5;
yMin = 2.0;

for i=10:-1:1
    
    % need to test xPt = xCoords(i) for dx = 0
    if (( xPt >= xCoords(i) && xPt <= xCoords(i+1) ) || ( xPt <= xCoords(i) && xPt >= xCoords(i+1) ) )
        % Deal with case where dx = 0
        m =  (yCoords(i+1) - yCoords(i)) / (xCoords(i+1) - xCoords(i));
        yPt = m*(xPt - xCoords(i)) + yCoords(i)
        break;
    end
end
   
yPts = max(yPt,yMin)