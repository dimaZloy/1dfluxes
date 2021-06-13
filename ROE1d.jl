
function RoeApproximation(ConsLeft::Array{Float64,1}, ConsRight::Array{Float64,1}, gamma::Float64)::Array{Float64,1}

	##// based on P.L. ROE
	##// Approximate Riemann Solvers, Parameter Vectors and Difference Scheme
	##// Journal of Computational Physics , 43 , PP 357-372 (1981)
	## 
	##// version v1.0.0
	##// without FIXED ENTROPY CORRECTION !!!!!!!!!!!!!!!!!!!

	
	
	rhoLeft::Float64 = ConsLeft[0+1];
	uLeft::Float64   = ConsLeft[1+1]/ConsLeft[0+1];
	eLeft::Float64   = ConsLeft[2+1];

	pLeft::Float64   = (gamma-1.0)*(eLeft-0.5*rhoLeft*uLeft*uLeft);
	cSoundLeft::Float64  = sqrt(gamma*pLeft/rhoLeft);
	HLeft::Float64  = cSoundLeft*cSoundLeft/(gamma-1)+0.5*uLeft*uLeft;     ##   % H=CvT+p/rho

	rhoRight::Float64 = ConsRight[0+1];
	uRight::Float64   = ConsRight[1+1]/ConsRight[0+1];
	eRight::Float64   = ConsRight[2+1];
	pRight::Float64   = (gamma-1.0)*(eRight-0.5*rhoRight*uRight*uRight);
	cSoundRight::Float64 = sqrt(gamma*pRight/rhoRight);
	HRight::Float64  = cSoundRight*cSoundRight/(gamma-1)+0.5*uRight*uRight;



	dCons1::Float64 = rhoRight - rhoLeft;
	dCons2::Float64 = rhoRight*uRight - rhoLeft*uLeft;
	dCons3::Float64 = eRight - eLeft;

	FluxLeft1::Float64 = rhoLeft*uLeft;
	FluxLeft2::Float64 = rhoLeft*uLeft*uLeft + pLeft;
	FluxLeft3::Float64 = (eLeft+pLeft)*uLeft;
	
	FluxRight1::Float64 = rhoRight*uRight;
	FluxRight2::Float64 = rhoRight*uRight*uRight + pRight;
	FluxRight3::Float64 = (eRight+pRight)*uRight;

	D::Float64  = sqrt(rhoRight/rhoLeft);
	rho::Float64  = sqrt(rhoRight*rhoLeft);;
	u::Float64  =(D*uRight+uLeft)/(D+1);
	H::Float64  =(D*HRight+HLeft)/(D+1);

	cSound2::Float64 = (gamma-1)*(H-0.5*u*u);
	cSound::Float64  =  sqrt(cSound2);
	##cSound=(D.*cSoundRight+cSoundLeft)./(D+1);

	EV1::Float64 = abs(u-cSound);
	EV2::Float64 = abs(u);
	EV3::Float64 = abs(u+cSound);          
	
	## EV --- EigenValue
	## EV1=EV1+(EV1.^2+Epson^2)/(2*Epson).*(EV1<=Epson);
	## EV2=EV2+(EV2.^2+Epson^2)/(2*Epson).*(EV2<=Epson);
	## EV3=EV3+(EV3.^2+Epson^2)/(2*Epson).*(EV3<=Epson);
		


	r11::Float64=1.0;
    r12::Float64=1.0;
    r13::Float64=1.0;            
	r21::Float64=u-cSound;
    r22::Float64=u;
    r23::Float64=u+cSound;   
	r31::Float64=H-u*cSound;
    r32::Float64=0.5*u*u;
    r33::Float64=H+u*cSound;


	b2::Float64 = (gamma-1.0)/cSound/cSound;
	b1::Float64 = 0.5*b2*u*u;

	l11::Float64 = 0.5*(b1+u/cSound);
	l21::Float64 = 1.0 - b1;
	l31::Float64 = 0.5*(b1-u/cSound);

	l12::Float64 = 0.5*(-b2*u -1.0/cSound);
	l22::Float64 = b2*u;
	l32::Float64 = 0.5*(-b2*u +1.0/cSound);

	l13::Float64 = 0.5*b2;
	l23::Float64 = -b2;
	l33::Float64 = 0.5*b2;

   
	RELdCons1::Float64=(EV1*l11*r11 + EV2*l21*r12 + EV3*l31*r13)*dCons1 + (EV1*l12*r11 + EV2*l22*r12 + EV3*l32*r13)*dCons2 + (EV1*l13*r11 + EV2*l23*r12 + EV3*l33*r13)*dCons3;
    
	RELdCons2::Float64=(EV1*l11*r21 + EV2*l21*r22 + EV3*l31*r23)*dCons1 + (EV1*l12*r21 + EV2*l22*r22 + EV3*l32*r23)*dCons2 + (EV1*l13*r21 + EV2*l23*r22 + EV3*l33*r23)*dCons3;
   
	RELdCons3::Float64=(EV1*l11*r31 + EV2*l21*r32 + EV3*l31*r33)*dCons1 + (EV1*l12*r31 + EV2*l22*r32 + EV3*l32*r33)*dCons2 + (EV1*l13*r31 + EV2*l23*r32 + EV3*l33*r33)*dCons3;
      
	Flux = zeros(Float64,4);
	Flux[1] = 0.5*(FluxLeft1+FluxRight1-RELdCons1);   
	Flux[2] = 0.5*(FluxLeft2+FluxRight2-RELdCons2);  
	Flux[3] = 0.5*(FluxLeft3+FluxRight3-RELdCons3);
	Flux[4] = abs(u)+cSound;
	
	return Flux;

end




function ROE1d(UPHSCELLS::Array{Float64,2}, Gamma::Float64)::Array{Float64,2}



		consLeft = zeros(3);
		consRight = zeros(3);
		flux = zeros(4); 

		#CMatrix UCONSCELLS(get_num_cells()+2,get_num_variables());
		#UCONSCELLS.Fill();
		#phs2cns2(UPHSCELLS,&UCONSCELLS);

		UCONSCELLS = phs2cns(UPHSCELLS);
		N = size(UPHSCELLS,1);
		nVar = size(UPHSCELLS,2);
		FLUXES = zeros(N,nVar+1);

		@simd for i = 1:N-1
		
			consLeft[0+1] = UCONSCELLS[i,1];
			consLeft[1+1] = UCONSCELLS[i,2];
			consLeft[2+1] = UCONSCELLS[i,3];
			
			consRight[0+1] = UCONSCELLS[i+1,1];
			consRight[1+1] = UCONSCELLS[i+1,2];
			consRight[2+1] = UCONSCELLS[i+1,3];

			#consLeft  = UCONSCELLS[i,1:3];
			#consRight = UCONSCELLS[i+1,1:3];
			
			flux = RoeApproximation(consLeft,consRight, Gamma);
			
			FLUXES[i,1] =  flux[0+1];
			FLUXES[i,2] =  flux[1+1];
			FLUXES[i,3] =  flux[2+1]; 
			FLUXES[i,4] =  flux[3+1];  ## vMaxU[i] 

			#FLUXES[i,1:4] =  flux;

		end

		return FLUXES;
end

function ROE1dtmp(UCONSCELLS::Array{Float64,2}, Gamma::Float64)::Array{Float64,2}



		consLeft = zeros(Float64,3);
		consRight = zeros(Float64,3);
		flux = zeros(Float64,4); 

		#CMatrix UCONSCELLS(get_num_cells()+2,get_num_variables());
		#UCONSCELLS.Fill();
		#phs2cns2(UPHSCELLS,&UCONSCELLS);

		#UCONSCELLS = phs2cns(UPHSCELLS);
		N = size(UCONSCELLS,1);
		nVar = size(UCONSCELLS,2);
		FLUXES = zeros(Float64,N,nVar+1);

		for i = 1:N-1
		
			consLeft[0+1] = UCONSCELLS[i,1];
			consLeft[1+1] = UCONSCELLS[i,2];
			consLeft[2+1] = UCONSCELLS[i,3];
			
			consRight[0+1] = UCONSCELLS[i+1,1];
			consRight[1+1] = UCONSCELLS[i+1,2];
			consRight[2+1] = UCONSCELLS[i+1,3];

			#consLeft  = UCONSCELLS[i,1:3];
			#consRight = UCONSCELLS[i+1,1:3];
			
			flux = RoeApproximation(consLeft,consRight, Gamma);
			
			FLUXES[i,1] =  flux[0+1];
			FLUXES[i,2] =  flux[1+1];
			FLUXES[i,3] =  flux[2+1]; 
			FLUXES[i,4] =  flux[3+1];  ## vMaxU[i] 

			#FLUXES[i,1:4] =  flux;

		end

		return FLUXES;
end

