function GetForces(PhiCutoff,Epsilon,sigma)
global x nAtoms Fx Phi

% n = 1;
for i = 1:nAtoms
    Fx(i) = 0;
    Phi(i) = 0;
%     if i == 10
%         plot(x(i),y(i),'o','markers',48);
%         hold on
%     end

    for j = 1:nAtoms
        if i == j, continue; end
        dx = x(i) - x(j);
        r = dx;

        if r > PhiCutoff, continue, end

        [aPhi dPhidr] = LJPot(r, Epsilon, sigma);

        
%         dx =  r*cos(ang)
%         dy =  r*sin(ang)

        dFx = - dPhidr;
%         if i == 10
%             plot(x(j),y(j),'ro','markers',48);
%         end

        Phi(i) = Phi(i) + aPhi;
        Fx(i) = Fx(i) + dFx;
    end
end

hold off

end

