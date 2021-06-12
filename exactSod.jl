


function exactSod(gamma::Float64,rho1::Float64,u1::Float64,p1::Float64,rho4::Float64,u4::Float64,p4::Float64,tEnd::Float64,x::Array{Float64,1})
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##% Riemann Solver for solving shoc-tube problems
##%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##% This programs was modified by Manuel Diaz, and is based on the code of 
##% [1]  P. Wesseling. PRINCIPLES OF COMPUTATIONAL FLUID DYNAMICS
##% Springer-Verlag, Berlin etc., 2001. ISBN 3-540-67853-0
##% See http://dutita0.twi.tudelft.nl/nw/users/wesseling/
##%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##% NOTE:
##% A Cavitation Check is the is incorporated in the code. It further
##% prevents plotting for possible but physically unlikely case of expansion
##% shocks. 
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##% INPUT VARIABLES: 
##% Problem definition: Conditions at time t=0
##%   rho1, u1, p1
##%   rho4, u4, p4
##% 'tEnd' and 'n' are the final solution time and the gas DoFs.
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## return [x,rho,ux,p,e,t,Mach,entro]=


##% Gamma values

##gamma=(n+2)/n; 
alpha::Float64=(gamma+1.0)/(gamma-1.0);

##% Assumed structure of exact solution
##%
##%    \         /      |con |       |s|
##%     \   f   /       |tact|       |h|
##% left \  a  /  state |disc| state |o| right
##% state \ n /    2    |cont|   3   |c| state
##%   1    \ /          |tinu|       |k|   4
##%         |           |ity |       | |

PRL::Float64 = p4/p1;
cright::Float64 = sqrt(gamma*p4/rho4); 
cleft::Float64  = sqrt(gamma*p1/rho1);
CRL::Float64 = cright/cleft;
MACHLEFT::Float64 = (u1-u4)/cleft;

##% Basic shock tube relation equation (10.51)
##f = @(P) (1+MACHLEFT*(gamma-1)/2-(gamma-1)*CRL*(P-1)/sqrt(2*gamma*(gamma-1+(gamma+1)*P)))^(2*gamma/(gamma-1))/P-PRL;
f(P) = (1.0 + MACHLEFT*(gamma-1.0)/2.0-(gamma-1.0)*CRL*(P-1.0)/sqrt(2.0*gamma*(gamma-1.0+(gamma+1.0)*P)))^(2.0*gamma/(gamma-1.0))/P-PRL;

##% solve for P = p34 = p3/p4
p34::Float64 = fzero(f,3);

p3::Float64 = p34*p4;
rho3::Float64 = rho4*(1.0+alpha*p34)/(alpha+p34); 
rho2::Float64 = rho1*(p34*p4/p1)^(1.0/gamma);
u2::Float64 = u1-u4+(2.0/(gamma-1.0))*cleft*(1.0-(p34*p4/p1)^((gamma-1.0)/(2.0*gamma)));
c2::Float64 = sqrt(gamma*p3/rho2);
spos::Float64 = 0.5 + tEnd*cright*sqrt((gamma-1.0)/(2.0*gamma)+(gamma+1.0)/(2.0*gamma)*p34)+tEnd*u4;

x0::Float64 = 0.5;
conpos::Float64 =x0 + u2*tEnd+tEnd*u4;	##% Position of contact discontinuity
pos1::Float64 = x0 + (u1-cleft)*tEnd;	##% Start of expansion fan
pos2::Float64 = x0 + (u2+u4-c2)*tEnd;	##% End of expansion fan

##% Plot structures
##%x = 0:0.002:1; % <--------- now x is defined as an input!
p = zeros(Float64,size(x)); 
ux= zeros(Float64,size(x)); 
rho = zeros(Float64,size(x));
Mach = zeros(Float64,size(x));  
cexact = zeros(Float64,size(x));

N::Int64 = length(x);

solution = zeros(Float64,N,3);

for i = 1:N
    if x[i] <= pos1
        p[i] = p1;
        rho[i] = rho1;
        ux[i] = u1;
        cexact[i] = sqrt(gamma*p[i]/rho[i]);
        Mach[i] = ux[i]/cexact[i];
    elseif x[i] <= pos2
        p[i] = p1*(1.0+(pos1-x[i])/(cleft*alpha*tEnd))^(2.0*gamma/(gamma-1.0));
        rho[i] = rho1*(1.0+(pos1-x[i])/(cleft*alpha*tEnd))^(2.0/(gamma-1.0));
        ux[i] = u1 + (2.0/(gamma+1))*(x[i]-pos1)/tEnd;
        cexact[i] = sqrt(gamma*p[i]/rho[i]);
        Mach[i] = ux[i]/cexact[i];
    elseif x[i] <= conpos
        p[i] = p3;
        rho[i] = rho2;
        ux[i] = u2+u4;
        cexact[i] = sqrt(gamma*p[i]/rho[i]);
        Mach[i] = ux[i]/cexact[i];
    elseif x[i] <= spos
        p[i] = p3;
        rho[i] = rho3;
        ux[i] = u2+u4;
        cexact[i] = sqrt(gamma*p[i]/rho[i]);
        Mach[i] = ux[i]/cexact[i];
    else
        p[i] = p4;
        rho[i] = rho4;
        ux[i] = u4;
        cexact[i] = sqrt(gamma*p[i]/rho[i]);
        Mach[i] = ux[i]/cexact[i];
    end
end
##entro = log(p./rho.^gamma); ##	% entropy
##e = p./((gamma-1).*rho);    ##	% internal energy
##t = 2/n.*e;                 ##% temperature

	##return x,rho,ux,p  ,e,t,Mach,entro
	##return x,rho,ux,p 
	solution[:,1] = rho;
	solution[:,2] = ux;
	solution[:,3] = p;
	
	
	return x, solution;
end



function computeExactSod1D(N::Int64)

	println("Computing 1D Sod exact solution  ");

	nVar::Int64 = 3;
	xFirst::Float64 = 0.0;
	xLast::Float64 = 1.0;


	dx::Float64 = (xLast-xFirst)/N;	

	xProblem =(xFirst + dx/2.0):dx:(xLast - dx/2.0);

	xFigure = zeros(Float64,N+2);

	for i=2:N+1
		xFigure[i] = xProblem[i-1];
	end 	

	xFigure[1] = xProblem[1]- dx;
	xFigure[end] = xProblem[end] + dx;
	
	
	Gamma::Float64 = 1.4;
	Rgas::Float64 = 8314.3/29.0;

	## boundary conditions for SOD problem: 
	
	phs_leftBC = zeros(Float64, 3);
	phs_rightBC = zeros(Float64, 3);

	phs_leftBC[1] = 1.0;
	phs_leftBC[2] = 0.0;
	phs_leftBC[3] = 1.0e+2;
	
	phs_rightBC[1] = 0.125;
	phs_rightBC[2] = 0.0;
	phs_rightBC[3] = 1.0e+1;
	
	xExact,exactSolution = exactSod(
		Gamma,phs_leftBC[1], phs_leftBC[2], phs_leftBC[3], phs_rightBC[1],phs_rightBC[2],phs_rightBC[3],0.02, xFigure);

		
	return xExact, exactSolution;
		
end

