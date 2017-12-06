function [s] = cartesian2Track(state,Track)
	% Extract components from the state vector
	x = state(1);
	y = state(3);
    
    % Do a rough cut
    roughDistLimit = 30;
	roughDist = sqrt((Track.cline(1,:)-x).^2+(Track.cline(2,:)-y).^2); 
    firstIdx = find(roughDist <= roughDistLimit,1,'First');
    lastIdx = find(roughDist <= roughDistLimit,1,'Last');
    sLims = Track.arc_s([firstIdx,lastIdx]);
    
    % Create three polynomials with the lines
    trackPoly = TrackPolyFit(Track,mean(Track.arc_s(firstIdx:lastIdx)));

    % Compute the distance polynomials
    [minLeftDist,   minLeftS]   = findMinimumLineDistance(trackPoly.leftXCoeffs,    trackPoly.leftYCoeffs,  sLims);    
    [minCenterDist, minCenterS] = findMinimumLineDistance(trackPoly.centerXCoeffs,  trackPoly.centerYCoeffs,sLims);    
    [minRightDist,  minRightS]  = findMinimumLineDistance(trackPoly.rightXCoeffs,   trackPoly.rightYCoeffs, sLims);    
    
    % Determine the minimum s value
    minSVals = [minLeftS,minCenterS,minRightS];
    [~, minDistIdx] = min([minLeftDist,minCenterDist,minRightDist]);
    s = minSVals(minDistIdx);
    
    function [minDist, minS] = findMinimumLineDistance(xPolyCoeffs,yPolyCoeffs,sLims)
        distXPoly = xPolyCoeffs + [0,0,-x];
        distYPoly = yPolyCoeffs + [0,0,-y];
        sqDistPoly = TrackPolyFit.square2OrderPoly(distXPoly) + TrackPolyFit.square2OrderPoly(distYPoly);
        sqDistDerPoly = polyder(sqDistPoly);
        distRoots =  roots(sqDistDerPoly);
        [~,realRootIdx] = min(abs(imag(distRoots)));
        
        % Check the endpoints to ensure minimum
        possibleSVals = [sLims,real(distRoots(realRootIdx))];
        [~,sIdx] = min(polyval(sqDistPoly,possibleSVals));
        minS = possibleSVals(sIdx);
        minDist = sqrt(polyval(sqDistPoly,minS));
    end
end
