


#physLeft = [1.0 290.0 0.0 7143.0];
#physRight = [1.7 261.72 -51.62 15282];

#include("AUSMflux.jl");

#flux1 = compute_1D_ARBITRARY_INVISCID_AUSM_PLUS_FLUX_from_UPHYS(physLeft,physRight, -1.0, 0.0, 0.01 ,1.4);
#flux2 = compute_1D_ARBITRARY_INVISCID_AUSM_PLUS_FLUX_from_UPHYS(physLeft,physRight,  1.0, 0.0, 0.01 ,1.4);
#flux3 = compute_1D_ARBITRARY_INVISCID_AUSM_PLUS_FLUX_from_UPHYS(physLeft,physRight,  0.0,  1.0, 0.01 ,1.4);
#flux4 = compute_1D_ARBITRARY_INVISCID_AUSM_PLUS_FLUX_from_UPHYS(physLeft,physRight,  0.0, -1.0, 0.01 ,1.4);


#@time compute_1D_ARBITRARY_INVISCID_AUSM_PLUS_FLUX_from_UPHYS(physLeft,physRight, -1.0, 0.0, 0.01 ,1.4)

using PyPlot;
using Roots;
using BenchmarkTools;

include("auxFluxes.jl");
include("exactSod.jl");
include("ROE1d.jl");
include("RIEMANN1d.jl");
include("AUSM1d.jl");



function computeSod1D(N::Int64, fluxApproxMethod)

## N - number of grid points in x direction
## fluxApproximationMethod: 
##  0 - ROE
##  1 - AUSM+
##  2 - AUSM+Up
##  3 - RIEMANN

	#N = 198;
	NCells = N+2;
	nVar::Int64 = 3;
	xFirst::Float64 = 0.0;
	xLast::Float64 = 1.0;
	dx::Float64 = (xLast-xFirst)/N;	

	xProblem =(xFirst + dx/2.0):dx:(xLast - dx/2.0);
	##display(size(xProblem))
	xFigure = zeros(Float64, N+2);
	for i=2:N+1
	   xFigure[i] = xProblem[i-1];
	end 	
	xFigure[1] = xProblem[1]-dx;
	xFigure[end] = xProblem[end]+dx;
	
	
	## calculate left and right cell areas: 
	
	SLSide = zeros(Float64,NCells);
	SRSide = zeros(Float64,NCells);
	VRatio = zeros(Float64,NCells);

	for K = 1:NCells-1
		##SLSide[K]= pi*YSideCells[K]*YSideCells[K];
		SLSide[K]= 1.0;
	end

	for M = 2:NCells
		SRSide[M-1]= 1.0;
	end


	## Volume of cell
	VCells=(SRSide.+SLSide).*dx.*0.5;

	
	for i = 1:NCells
		VRatio[i] =1.0/VCells[i];
	end
	
	##  Step on x-directions
	HX=(xLast-xFirst)/NCells;
	
	CFL::Float64 = 0.5;
	N_iterations::Int64 = 1000;
	Total_time_duration::Float64 = 0.02; 
	
	
	Gamma::Float64 = 1.4;
	Rgas::Float64 = 8314.3/29.0;

	
	phs_leftBC = zeros(Float64,3);
	phs_rightBC = zeros(Float64,3);

	phs_leftBC[1] = 1.0;
	phs_leftBC[2] = 0.0;
	phs_leftBC[3] = 1.0e+2;
	
	phs_rightBC[1] = 0.125;
	phs_rightBC[2] = 0.0;
	phs_rightBC[3] = 1.0e+1;
	

	UconsCellsOld = zeros(Float64, N+2,nVar);
	UconsCellsNew = zeros(Float64, N+2,nVar);
	
	
	#RFLUX = zeros(N+2,nVar+1); ## tmp to hold fluxes
	UFLUX = zeros(Float64, N+2,nVar+1); ## tmp to hold fluxes
	##UPHS = zeros(N+2,nVar); ## tmp to hold fluxes
	
	aSound = zeros(Float64, N+2);
	vMaxU  = zeros(Float64, N+2);

	## initial conditions : t = 0.0; 
	## Vector of physical variables

	UphysCells = zeros(Float64, N+2,nVar);

	for i=1:Int32(N/2)-1
		UphysCells[i,1] = phs_leftBC[1];
		UphysCells[i,2] = phs_leftBC[2];
		UphysCells[i,3] = phs_leftBC[3];
	end

	for i= Int32(N/2):N+2	
		UphysCells[i,1] = phs_rightBC[1];
		UphysCells[i,2] = phs_rightBC[2];
		UphysCells[i,3] = phs_rightBC[3];
	end

	## convert physical variables -> conservative 
	## t = 0.0 
	UconsCellsOld = phs2cns(UphysCells,Gamma);
	#display(UconsCellsOld)
	UconsCellsNew = UconsCellsOld; 

	for i=1:N+2
	   aSound[i] = sqrt(Gamma*UphysCells[i,3]/UphysCells[i,1]);
	   vMaxU[i] = abs(UphysCells[i,2]) + aSound[i];	
	end


	VMax::Float64 = calcMagMaxVelocity(aSound, vMaxU, UphysCells, Gamma); 

	#Side::Float64  = 1.0/dx;
	#VCell::Float64 = dx*Side;
	#VRatio::Float64 = 1.0/VCell;
	
	betta::Float64 = 1.0;
	

	fluxName  = " ";
	if (fluxApproxMethod == 0)
		fluxName = "approximation by Roe ";
	elseif (fluxApproxMethod == 1)
		fluxName = "approximation by AUSM+ method ";
	elseif (fluxApproxMethod == 2)
		fluxName = "approximation by AUSM+Up method ";
	elseif (fluxApproxMethod == 3)
		fluxName = "approximation by Riemann ";
	else
		fluxName = "flux method is not defined";
	end
	
	#println("Solving 1D Sod problem using ", fluxName);
	#println("Starting calculations...");

	ni::Int64 = 0; ## iterations counter
	c::Int64 = 0; ## verbosity counter
	calcTime::Float64 = 0.0;  ## physical simulation time
	verbosity ::Int64= 0; 
	Tau::Float64 = CFL*dx/VMax;
	
	while( ni <= N_iterations && calcTime <= Total_time_duration )
	
		#if (c == verbosity)
		#	println("iteration:\t", ni, "\t tau = ", Tau, "\t total time: ", time); 
		# c = 0;
		#end	
		
		if (fluxApproxMethod == 0)
			UFLUX = ROE1dtmp(UconsCellsOld,Gamma);
		elseif (fluxApproxMethod == 1)
			UFLUX = AUSMplus1d(UphysCells,Gamma);
			#UFLUX = AUSM_DL1d(UphysCells,Gamma);
		elseif (fluxApproxMethod == 2)					
			UFLUX = AUSMplusUp1d(UphysCells,Gamma);
		elseif (fluxApproxMethod == 3)
			UFLUX  =  GOD1d(NCells,UphysCells,Gamma); #fluxes are PHS variables
		else
			#UFLUX = zeros(Float64,N+2,nVar);
		end
		

		for i = 2:N-1
		
			#HMid = (Gamma-1.0)*(UconsCellsOld[i,3] -UconsCellsOld[i,2]^2/2.0/UconsCellsOld[i,1])*(SRSide[i]-SLSide[i]);
			HMid = 0.0;
			UconsCellsNew[i,1]  =  UconsCellsOld[i,1] - betta * Tau* VRatio[i]* ( UFLUX[i,1]*SRSide[i] - UFLUX[i-1,1]*SLSide[i]) ;
			UconsCellsNew[i,2]  =  UconsCellsOld[i,2] - betta * Tau* VRatio[i]* ( UFLUX[i,2]*SRSide[i] - UFLUX[i-1,2]*SLSide[i]) + Tau*VRatio[i]*HMid ;
			UconsCellsNew[i,3]  =  UconsCellsOld[i,3] - betta * Tau* VRatio[i]* ( UFLUX[i,3]*SRSide[i] - UFLUX[i-1,3]*SLSide[i]) ;	

			
		end
	
		UconsCellsOld = UconsCellsNew; 

		UphysCells = cns2phs(UconsCellsOld,Gamma);
		UphysCells[1,1:3] = phs_leftBC;
		UphysCells[end,1:3] = phs_rightBC;
		UconsCellsOld = phs2cns(UphysCells,Gamma);
				
	
		for i=1:NCells
			aSound[i] = sqrt(Gamma*UphysCells[i,3]/UphysCells[i,1]);
			vMaxU[i] = abs(UphysCells[i,2]) + aSound[i];	
		end
		
		
		VMax = calcMagMaxVelocity(aSound, vMaxU, UphysCells, Gamma); 
		
		Tau  = CFL*dx/VMax;
		
		#display(Tau)
		
		
		calcTime = calcTime +  Tau;
		ni += 1;
		c += 1;
		
		
	end
	
	#println("done...");
	

	return xFigure, UphysCells;


end


function prime()



@btime computeSod1D(198, 0);
@btime computeSod1D(198, 1);
@btime computeSod1D(198, 2);
@btime computeSod1D(198, 3); 


@time x1, theorySolution = computeExactSod1D(198);
@time x2, RoeSolution = computeSod1D(198, 0);
@time x3, AUSMplusSolution = computeSod1D(198, 1);
@time x4, AUSMplusUpSolution = computeSod1D(198, 2);
@time x5, RiemannSolution = computeSod1D(198, 3); 

		
		
figure(2)
clf();

subplot(3,1,1)
plot(x1, theorySolution[:,1],"-k", label="theory");
plot(x1, RoeSolution[:,1],"-.g",label="Roe");
plot(x1, AUSMplusSolution[:,1],"--b",label="AUSM+");
plot(x1, AUSMplusUpSolution[:,1],"--c",label="AUSM+Up");
plot(x1, RiemannSolution[:,1],"-.r",label="Riemann");

xlabel("x")
ylabel("rho")
xlim(0.0,1.0);
legend();
subplot(3,1,2)
plot(x1, theorySolution[:,2],"-k", label="theory");
plot(x1, RoeSolution[:,2],"-.g",label="Roe");
plot(x1, AUSMplusSolution[:,2],"--b",label="AUSM+");
plot(x1, AUSMplusUpSolution[:,2],"--c",label="AUSM+Up");
plot(x1, RiemannSolution[:,2],"-.r",label="Riemann");

xlabel("x")
ylabel("U")
xlim(0.0,1.0);

legend();
subplot(3,1,3)
plot(x1, theorySolution[:,3],"-k", label="theory");
plot(x1, RoeSolution[:,3],"-.g",label="Roe");
plot(x1, AUSMplusSolution[:,3],"--b",label="AUSM+");
plot(x1, AUSMplusUpSolution[:,3],"--c",label="AUSM+Up");
plot(x1, RiemannSolution[:,3],"-.r",label="Riemann");

xlabel("x")
ylabel("P")
xlim(0.0,1.0);

legend();

end 
	

function auxTestFluxes()

	Gamma = 1.4;
	Rgas = 8314.3/29.0;

	N = 198;
	nVar = 3;
	NCells = 200;
	
	phs_leftBC = zeros(Float64,3);
	phs_rightBC = zeros(Float64,3);
	
	cns_leftBC = zeros(Float64,3);
	cns_rightBC = zeros(Float64,3);


	phs_leftBC[1] = 1.0;
	phs_leftBC[2] = 0.0;
	phs_leftBC[3] = 1.0e+2;
	
	phs_rightBC[1] = 0.125;
	phs_rightBC[2] = 0.0;
	phs_rightBC[3] = 1.0e+1;
	
	cns_leftBC = phs2cnsFlux3(phs_leftBC,Gamma);
	cns_rightBC = phs2cnsFlux3(phs_rightBC,Gamma);

	UconsCellsOld = zeros(Float64, N+2,nVar);
	UconsCellsNew = zeros(Float64, N+2,nVar);
	
	
	#RFLUX = zeros(N+2,nVar+1); ## tmp to hold fluxes
	UFLUX = zeros(Float64, N+2,nVar+1); ## tmp to hold fluxes
	##UPHS = zeros(N+2,nVar); ## tmp to hold fluxes
	
	aSound = zeros(Float64, N+2);
	vMaxU  = zeros(Float64, N+2);

	## initial conditions : t = 0.0; 
	## Vector of physical variables

	UphysCells = zeros(Float64, N+2,nVar);

	for i=1:Int32(N/2)-1
		UphysCells[i,1] = phs_leftBC[1];
		UphysCells[i,2] = phs_leftBC[2];
		UphysCells[i,3] = phs_leftBC[3];
	end

	for i= Int32(N/2):N+2	
		UphysCells[i,1] = phs_rightBC[1];
		UphysCells[i,2] = phs_rightBC[2];
		UphysCells[i,3] = phs_rightBC[3];
	end

	UconsCellsOld = phs2cns(UphysCells,Gamma);


	
	flux1 = RoeApproximation(cns_leftBC,cns_rightBC, Gamma);
	flux2 = computeRiemannFlux(phs_leftBC[1],  phs_leftBC[2], phs_leftBC[3], phs_rightBC[1], phs_rightBC[2], phs_rightBC[3], Gamma);
	flux3 = AUSMplusApproximation(phs_leftBC,phs_rightBC, Gamma);
	flux4 = AUSMplusUpApproximation(phs_leftBC,phs_rightBC, Gamma);
	
	
	flux5 = AUSMapproximationDL(phs_leftBC,phs_rightBC, Gamma);
	#flux5 = phs2cnsFlux3(flux5tmp,Gamma);


	println(flux1)
	println(flux2)
	println(flux3)
	println(flux4)
	println(flux5)
	
end


	







