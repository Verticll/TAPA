function ODMDGEN(LAtoms, X0, VX0, InitDist, Temp, Type)
global C
global x AtomSpacing
global nAtoms
global AtomType Vx Mass0 Mass1 Mass2

if Type == 0
    Mass = Mass0;
else if Type == 1
    Mass = Mass1;
    else
    Mass = Mass2;
    end
end

L = (LAtoms - 1) * AtomSpacing;


numAtoms = LAtoms;

xp(1, :) = linspace(0, L, LAtoms);


x(nAtoms + 1:nAtoms+LAtoms) = xp-L/2;


x(nAtoms * LAtoms + 1:nAtoms + LAtoms) = xp - L / 2;


x(nAtoms + 1:nAtoms + numAtoms) = x(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + X0;


AtomType(nAtoms + 1:nAtoms + numAtoms) = Type;

if Temp == 0
    Vx(nAtoms + 1:nAtoms + numAtoms) = 0;
else
    std0 = sqrt(C.kb * Temp / Mass);

    Vx(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
end

Vx(nAtoms + 1:nAtoms + numAtoms) = Vx(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vx(nAtoms + 1:nAtoms + numAtoms)) + VX0;

nAtoms = nAtoms + numAtoms;

end
