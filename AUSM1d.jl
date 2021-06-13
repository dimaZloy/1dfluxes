

# rev 23-Sep-2018
# DL 
# AUSM+  flux splitting
# based on paperby Meng-Sing Liou "A Sequel to AUSM: AUSM+"
# JOURNAL OF COMPUTATIONAL PHYSICS 129, 364-382 (1996)
# there is a mistake in the article (I guess based on tests)
# in the equation 19b fro Mbetta - instead of 0.5 should use 0.25 !!!!!!!!!!!!!!!!!!!!
# find in this article:
# Azevedo, Korzenowski 
# An assessment of unstructured grid FV schemes for cold gas hypersonic flow calculations. 
# Journal of Aerospace Technology and Management, V1,n2,2009



##function get_gas_epsilon(p::Float64, rho::Float64, gamma::Float64)::Float64
## return p/rho/(gamma-1.0); 
##end

@inline function P_m(M::Float64,AUSM_ALFA::Float64)::Float64
	if (abs(M)>=1.0)
		return 0.5*(1.0-sign(M));
	else
		return Palfa_m(M,AUSM_ALFA);
	end
end

@inline function P_p(M::Float64,AUSM_ALFA::Float64)::Float64
	if (abs(M)>=1.0)
		return 0.5*(1.0+sign(M));
	else
		return Palfa_p(M,AUSM_ALFA);
	end
end

@inline function Palfa_m(M::Float64,AUSM_ALFA::Float64)::Float64
	return  0.25*(M-1.0)*(M-1.0)*(2.0+M)-AUSM_ALFA*M*(M*M-1.0)*(M*M-1.0);
end

@inline function Palfa_p(M::Float64,AUSM_ALFA::Float64)::Float64
	return  0.25*(M+1.0)*(M+1.0)*(2.0-M)+AUSM_ALFA*M*(M*M-1.0)*(M*M-1.0);
end

@inline function Mbetta_p(M::Float64,AUSM_BETTA::Float64)::Float64
	return  0.25*(M+1.0)*(M+1.0)+AUSM_BETTA*(M*M-1.0)*(M*M-1.0);
end

@inline function Mbetta_m(M::Float64,AUSM_BETTA::Float64)::Float64
	return -0.25*(M-1.0)*(M-1.0)-AUSM_BETTA*(M*M-1.0)*(M*M-1.0);
end

@inline function M_p(M::Float64,AUSM_BETTA::Float64)::Float64
	if (abs(M)>=1.0)
		return 0.5*(M+abs(M));
	else
		return Mbetta_p(M,AUSM_BETTA);
	end	
end

@inline function M_m(M::Float64,AUSM_BETTA::Float64)::Float64
	if (abs(M)>=1.0)
		return  0.5*(M-abs(M));
	else
		return Mbetta_m(M,AUSM_BETTA);
	end
end



@inline function AUSMapproximationDL( PhsLeft::Array{Float64,1}, PhsRight::Array{Float64,1}, gamma::Float64)::Array{Float64,1}

	rhoL::Float64 = 	PhsLeft[1];
	_UL::Float64 = PhsLeft[2];
	_VL::Float64 = 0.0;
	PL::Float64 = PhsLeft[3];
	
	rhoR::Float64 = PhsRight[1];
	_UR::Float64 = PhsRight[2];
	_VR::Float64 = 0.0;
	PR::Float64 = PhsRight[3];
	
	nx::Float64 = 1.0;
	ny::Float64 = 0.0;
	
	

	VL_tilda::Float64   = _UL*nx + _VL*ny;
	TLeft::Float64      = _UL*ny - _VL*nx;
	
	VR_tilda::Float64   = _UR*nx + _VR*ny;
	TRight::Float64     = _UR*ny - _VR*nx;

	UMAG_L2::Float64 = (_UL*_UL + _VL*_VL);
	UMAG_R2::Float64 = (_UR*_UR + _VR*_VR);

	AUSM_BETTA::Float64 = 1.0/8.0;
	AUSM_ALFA::Float64  = 3.0/16.0; 


	#htL =  get_gas_epsilon(PL,rhoL,gamma)+ 0.5*(UMAG_L2) +  PL/rhoL; 
	#htR =  get_gas_epsilon(PR,rhoR,gamma)+ 0.5*(UMAG_R2) +  PR/rhoR;

	htL::Float64 =  PL/rhoL/(gamma-1.0) + 0.5*(UMAG_L2) +  PL/rhoL; 
	htR::Float64 =  PR/rhoR/(gamma-1.0) + 0.5*(UMAG_R2) +  PR/rhoR;

	#display(htL)
	#display(htR)

	#htN::Float64  = 0.5*(htL + htR - 0.5*(TLeft*TLeft + TRight*TRight));



	aL_dot::Float64 = sqrt(2*( gamma-1.0)/(gamma+1.0)*htL);
	aL_tilda::Float64 = aL_dot*min(1.0,aL_dot/abs(VL_tilda));
	aR_dot::Float64 = sqrt(2*( gamma-1.0)/(gamma+1.0)*htR);
	aR_tilda::Float64 = aR_dot*min(1.0,aR_dot/abs(VR_tilda));
	a12::Float64 = min(aL_tilda,aR_tilda);

	


	MLeft::Float64  = VL_tilda/a12;
	MRight::Float64 = VR_tilda/a12;

	

	m_dot12::Float64 = M_p(MLeft,AUSM_BETTA) + M_m(MRight,AUSM_BETTA);
	p12::Float64 = P_p(MLeft,AUSM_ALFA)*PL + P_m(MRight,AUSM_ALFA)*PR;

	
	m_dot12_p::Float64 = 0.5*(m_dot12+abs(m_dot12));
	m_dot12_m::Float64 = 0.5*(m_dot12-abs(m_dot12));

	
	p = zeros(Float64,3);
	F_LEFT = zeros(Float64,3);
	F_RIGHT = zeros(Float64,3);
	flux = zeros(Float64,3);
	

	p[1] = 0.0;
	p[2] = p12*nx;
	##p[3] = p12*ny; 
	p[3] = 0.0;	
	
	
	F_LEFT[1] = rhoL;
	F_LEFT[2] = rhoL*_UL;
	##F_LEFT[3] = rhoL*_VL;
	F_LEFT[3] = rhoL*htL;

	F_RIGHT[1] = rhoR;
	F_RIGHT[2] = rhoR*_UR;
	##F_RIGHT[3] = rhoR*_VR;
	F_RIGHT[3] = rhoR*htR; 
	
	
	##flux[1]  = (0.5*(URelTilde*(rhoLeft +rhoRight) -magURelTilde*(rhoRight -rhoLeft)))*magSf;
	##flux[2]  = (0.5*(URelTilde*(rhoULeft+rhoURight)-magURelTilde*(rhoURight-rhoULeft))+pTilde*normalVector)*magSf;
    ##flux[3]  = (0.5*(URelTilde*(rhoHLeft+rhoHRight)-magURelTilde*(rhoHRight-rhoHLeft))+pTilde*(dotX * normalVector))*magSf;
	
	
	flux[1] = -( a12*(m_dot12_p*F_LEFT[1] + m_dot12_m*F_RIGHT[1]) + p[1] );
	flux[2] = -( a12*(m_dot12_p*F_LEFT[2] + m_dot12_m*F_RIGHT[2]) + p[2] );
	#flux[3] = -( a12*(m_dot12_p*F_LEFT[3] + m_dot12_m*F_RIGHT[3]) + p[3] )*side;
	flux[3] = -( a12*(m_dot12_p*F_LEFT[3] + m_dot12_m*F_RIGHT[3]) + p[3] );
	
	
	
	return flux; 
	#return flux.*-1.0; 
	
	
end


@inline function AUSMplus1d(UPHSCELLS::Array{Float64,2}, Gamma::Float64)::Array{Float64,2}



		phsLeft = zeros(Float64,3);
		phsRight = zeros(Float64,3);
		flux = zeros(Float64,3); 


		N = size(UPHSCELLS,1);
		nVar = size(UPHSCELLS,2);
		FLUXES = zeros(Float64,N,nVar+1);

		for i = 1:N-1
		
			phsLeft[1] = UPHSCELLS[i,1];
			phsLeft[2] = UPHSCELLS[i,2];
			phsLeft[3] = UPHSCELLS[i,3];
			
			phsRight[1] = UPHSCELLS[i+1,1];
			phsRight[2] = UPHSCELLS[i+1,2];
			phsRight[3] = UPHSCELLS[i+1,3];

			flux = AUSMplusApproximation(phsLeft, phsRight, Gamma);
			
			
			FLUXES[i,1] =  flux[1];
			FLUXES[i,2] =  flux[2];
			FLUXES[i,3] =  flux[3]; 
			FLUXES[i,4] =  0.0;

		end

		return FLUXES;
end

@inline function AUSMplusUp1d(UPHSCELLS::Array{Float64,2}, Gamma::Float64)::Array{Float64,2}



		phsLeft = zeros(Float64,3);
		phsRight = zeros(Float64,3);
		flux = zeros(Float64,3); 


		N = size(UPHSCELLS,1);
		nVar = size(UPHSCELLS,2);
		FLUXES = zeros(Float64,N,nVar+1);

		for i = 1:N-1
		
			phsLeft[1] = UPHSCELLS[i,1];
			phsLeft[2] = UPHSCELLS[i,2];
			phsLeft[3] = UPHSCELLS[i,3];
			
			phsRight[1] = UPHSCELLS[i+1,1];
			phsRight[2] = UPHSCELLS[i+1,2];
			phsRight[3] = UPHSCELLS[i+1,3];

			flux = AUSMplusUpApproximation(phsLeft, phsRight, Gamma);	
			
			FLUXES[i,1] =  flux[1];
			FLUXES[i,2] =  flux[2];
			FLUXES[i,3] =  flux[3]; 
			FLUXES[i,4] =  0.0;

		end

		return FLUXES;
end


@inline function AUSM_DL1d(UPHSCELLS::Array{Float64,2}, Gamma::Float64)::Array{Float64,2}



		phsLeft = zeros(Float64,3);
		phsRight = zeros(Float64,3);
		flux = zeros(Float64,3); 


		N = size(UPHSCELLS,1);
		nVar = size(UPHSCELLS,2);
		FLUXES = zeros(Float64,N,nVar+1);

		for i = 1:N-1
		
			phsLeft[1] = UPHSCELLS[i,1];
			phsLeft[2] = UPHSCELLS[i,2];
			phsLeft[3] = UPHSCELLS[i,3];
			
			phsRight[1] = UPHSCELLS[i+1,1];
			phsRight[2] = UPHSCELLS[i+1,2];
			phsRight[3] = UPHSCELLS[i+1,3];

			flux = AUSMapproximationDL(phsLeft, phsRight, Gamma);	
			
			
			FLUXES[i,1] =  flux[1];
			FLUXES[i,2] =  flux[2];
			FLUXES[i,3] =  flux[3]; 
			FLUXES[i,4] =  0.0;

		end

		return FLUXES;
end



@inline function sqr(a::Float64)::Float64
	return a*a;
end

@inline function AUSMplusUpApproximation(PhsLeft::Array{Float64,1}, PhsRight::Array{Float64,1},Gamma::Float64)::Array{Float64,1}
   
    # AUSM+Up flux splitting method. 
	# Uses primitive variables as input and gives back conservative numerical fluxes.
    #Luo, H.; Baum, Joseph D. and Löhner R. "On the computation of multi-material flows using ALE formulation."
    #Journal of Computational Physics 194 (2004): 304–328.    
    #Meng-Sing Liou, "A sequel to AUSM, PartII: AUSM+ -up for all speeds"
    #Journal of Computational Physics 214 (2006): 137-170


	rhoLeft::Float64 = PhsLeft[1];
	ULeft::Float64 = PhsLeft[2];
	pLeft::Float64 = PhsLeft[3];
	
	rhoRight::Float64 = PhsRight[1];
	URight::Float64 = PhsRight[2];
    pRight::Float64= PhsRight[3];
    
	Sf::Float64 = 1.0; 
	magSf::Float64 = 1.0; 
	dotX::Float64 = 0.0;

	## alpha::Float64 = 3.0/16.0;
    beta::Float64  = 1.0/8.0;
	MaInf::Float64 = 10.0;
	Kp::Float64    = 0.25;
	Ku::Float64    = 0.75;
	sigma::Float64 = 1.0;
	
	kLeft::Float64 = 0.0; #pLeft/rhoLeft;
	kRight::Float64 = 0.0; #pRight/rhoRight;
	
	
	#kappaLeft::Float64 = Gamma;
    #kappaRight::Float64 = Gamma;

    ## bounding variables
    rhoMin::Float64 = eps();

    ## normal vector
    normalVector::Float64 = Sf/magSf;

    ## speed of sound, for left and right side, assuming perfect gas
	#aLeft::Float64  = sqrt(max(eps(),Gamma*pLeft /max(rhoLeft, rhoMin)));
    #aRight::Float64 = sqrt(max(eps(),Gamma*pRight/max(rhoRight,rhoMin)));
	
	aLeft::Float64  = sqrt(Gamma*pLeft / rhoLeft);
    aRight::Float64 = sqrt(Gamma*pRight/ rhoRight);
	
		
    rhoHLeft::Float64 = pLeft *Gamma /(Gamma -1.0) + rhoLeft *(0.5*ULeft*ULeft + kLeft);
    rhoHRight::Float64 = pRight*Gamma/(Gamma-1.0) + rhoRight*(0.5*URight*URight + kRight);

    rhoULeft::Float64  = rhoLeft *ULeft;
    rhoURight::Float64 = rhoRight*URight;

	
	qLeft::Float64  = ((ULeft  - dotX) * normalVector);
    qRight::Float64 = ((URight - dotX) * normalVector);


    aTilde::Float64 = 0.5*(aLeft+aRight);
	rhoTilde::Float64 = 0.5*(rhoLeft+rhoRight);
	
	sqrMaDash::Float64 = (sqr(qLeft)+sqr(qRight))/(2.0*sqr(aTilde));
	
	sqrMaZero::Float64 = min(1.0,max(sqrMaDash,(MaInf*MaInf)));
	MaZero::Float64    = sqrt(sqrMaZero);
	
	fa::Float64 = MaZero*(2.0-MaZero);
	
	alpha::Float64 = 3.0/16.0*(-4.0+5.0*sqr(fa));

    MaRelLeft::Float64  = qLeft /aTilde;
    MaRelRight::Float64 = qRight/aTilde;

    magMaRelLeft::Float64  = abs(MaRelLeft);
    magMaRelRight::Float64 = abs(MaRelRight);

    Ma1PlusLeft::Float64   = 0.5*(MaRelLeft +magMaRelLeft );
    Ma1MinusRight::Float64 = 0.5*(MaRelRight-magMaRelRight);

     
    Ma2PlusLeft::Float64   =  0.25*sqr(MaRelLeft +1.0);
    Ma2PlusRight::Float64  =  0.25*sqr(MaRelRight+1.0);
    Ma2MinusLeft::Float64  = -0.25*sqr(MaRelLeft -1.0);
    Ma2MinusRight::Float64 = -0.25*sqr(MaRelRight-1.0);

    Ma4BetaPlusLeft::Float64   = ((magMaRelLeft  >= 1.0) ? Ma1PlusLeft   : (Ma2PlusLeft  *(1.0-16.0*beta*Ma2MinusLeft)));
    Ma4BetaMinusRight::Float64 = ((magMaRelRight >= 1.0) ? Ma1MinusRight : (Ma2MinusRight*(1.0+16.0*beta*Ma2PlusRight)));
	
	Mp::Float64 = -Kp/fa*max(1.0-sigma*sqrMaDash,0.0)*(pRight-pLeft)/(rhoTilde*sqr(aTilde));

    P5alphaPlusLeft::Float64   = ((magMaRelLeft  >= 1.0) ?
            (Ma1PlusLeft/MaRelLeft)    : (Ma2PlusLeft  *(( 2.0-MaRelLeft) -16.0*alpha*MaRelLeft *Ma2MinusLeft )));
			
    P5alphaMinusRight::Float64 = ((magMaRelRight >= 1.0) ?
            (Ma1MinusRight/MaRelRight) : (Ma2MinusRight*((-2.0-MaRelRight)+16.0*alpha*MaRelRight*Ma2PlusRight)));
	    
	pU::Float64= -Ku*P5alphaPlusLeft*P5alphaMinusRight*(rhoLeft+rhoRight)*(fa*aTilde)*(qRight-qLeft);
	    
    MaRelTilde::Float64 = Ma4BetaPlusLeft + Ma4BetaMinusRight + Mp;
    pTilde::Float64 = pLeft*P5alphaPlusLeft + pRight*P5alphaMinusRight + pU;

    URelTilde::Float64 = MaRelTilde*aTilde;
    magURelTilde::Float64 = abs(MaRelTilde)*aTilde;
	
	flux = zeros(Float64,3);
	
    flux[1]  = (0.5*(URelTilde*(rhoLeft +rhoRight) -magURelTilde*(rhoRight -rhoLeft)))*magSf;
    flux[2]  = (0.5*(URelTilde*(rhoULeft+rhoURight)-magURelTilde*(rhoURight-rhoULeft))+pTilde*normalVector)*magSf;
    flux[3]  = (0.5*(URelTilde*(rhoHLeft+rhoHRight)-magURelTilde*(rhoHRight-rhoHLeft))+pTilde*(dotX * normalVector))*magSf;
	
	return flux;

end

@inline function AUSMplusApproximation(PhsLeft::Array{Float64,1}, PhsRight::Array{Float64,1},Gamma::Float64)::Array{Float64,1}
   

    # AUSM+ flux splitting method. 
	# Uses primitive variables as input and gives back conservative numerical fluxes.
    # Luo, H.; Baum, Joseph D. and Löhner R. "On the computation of multi-material flows using ALE formulation."
    # Journal of Computational Physics 194 (2004): 304–328.
   

	rhoLeft::Float64 = PhsLeft[1];
	ULeft::Float64 = PhsLeft[2];
	pLeft::Float64 = PhsLeft[3];
	
	rhoRight::Float64 = PhsRight[1];
	URight::Float64 = PhsRight[2];
    pRight::Float64= PhsRight[3];
    
	Sf::Float64 = 1.0; 
	magSf::Float64 = 1.0; 
	dotX::Float64 = 0.0;
	
	alpha = 3.0/16.0;
    beta = 1.0/8.0;
       
    rhoMin = eps();

    normalVector = Sf/magSf;
	
	#kappaLeft::Float64 = Gamma;
    #kappaRight::Float64 = Gamma;
	kLeft::Float64 = 0.0;
	kRight::Float64 = 0.0;
	
	aLeft  = sqrt(Gamma *pLeft / rhoLeft);
    aRight = sqrt(Gamma*pRight/ rhoRight);
	        
    rhoHLeft =  pLeft *Gamma /(Gamma -1.0) + rhoLeft *(0.5*(ULeft*ULeft) +kLeft);
    rhoHRight =  pRight*Gamma/(Gamma-1.0) + rhoRight*(0.5*(URight*URight)+kRight);

    rhoULeft  = rhoLeft *ULeft;
    rhoURight = rhoRight*URight;
	
    qLeft  = ((ULeft  - dotX) * normalVector);
    qRight = ((URight - dotX) * normalVector);

    aTilde = 0.5*(aLeft+aRight);

    MaRelLeft  = qLeft /aTilde;
    MaRelRight = qRight/aTilde;

    magMaRelLeft  = abs(MaRelLeft);
    magMaRelRight = abs(MaRelRight);

    Ma1PlusLeft   = 0.5*(MaRelLeft +magMaRelLeft );
    Ma1MinusRight = 0.5*(MaRelRight-magMaRelRight);

    Ma2PlusLeft   =  0.25*sqr(MaRelLeft +1.0);
    Ma2PlusRight  =  0.25*sqr(MaRelRight+1.0);
    Ma2MinusLeft  = -0.25*sqr(MaRelLeft -1.0);
    Ma2MinusRight = -0.25*sqr(MaRelRight-1.0);

    Ma4BetaPlusLeft   = ((magMaRelLeft  >= 1.0) ? Ma1PlusLeft   : (Ma2PlusLeft  *(1.0-16.0*beta*Ma2MinusLeft)));
    Ma4BetaMinusRight = ((magMaRelRight >= 1.0) ? Ma1MinusRight : (Ma2MinusRight*(1.0+16.0*beta*Ma2PlusRight)));

    P5alphaPlusLeft   = ((magMaRelLeft  >= 1.0) ?
            (Ma1PlusLeft/MaRelLeft)    : (Ma2PlusLeft  *(( 2.0-MaRelLeft) -16.0*alpha*MaRelLeft *Ma2MinusLeft )));
    P5alphaMinusRight = ((magMaRelRight >= 1.0) ?
            (Ma1MinusRight/MaRelRight) : (Ma2MinusRight*((-2.0-MaRelRight)+16.0*alpha*MaRelRight*Ma2PlusRight)));

    MaRelTilde = Ma4BetaPlusLeft + Ma4BetaMinusRight;
    pTilde = pLeft*P5alphaPlusLeft + pRight*P5alphaMinusRight;

    URelTilde = MaRelTilde*aTilde;
    magURelTilde = abs(MaRelTilde)*aTilde;
	
	flux = zeros(Float64,3);
	
    flux[1]  = (0.5*(URelTilde*(rhoLeft +rhoRight) -magURelTilde*(rhoRight -rhoLeft)))*magSf;
    flux[2]  = (0.5*(URelTilde*(rhoULeft+rhoURight)-magURelTilde*(rhoURight-rhoULeft))+pTilde*normalVector)*magSf;
    flux[3]  = (0.5*(URelTilde*(rhoHLeft+rhoHRight)-magURelTilde*(rhoHRight-rhoHLeft))+pTilde*(dotX * normalVector))*magSf;
	
	return flux;

end



