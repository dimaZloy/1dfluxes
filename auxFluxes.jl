
@inline function phs2cns(APhys::Array{Float64,2},gamma::Float64)::Array{Float64,2}

N = size(APhys,1);
ACons = zeros(Float64,N,3);

for i = 1:N
	ACons[i,1] = APhys[i,1];
	ACons[i,2] = APhys[i,1]*APhys[i,2];
	ACons[i,3] = APhys[i,3]/(gamma-1.0) + 0.5*APhys[i,1]*( APhys[i,2]*APhys[i,2] );
end 

return ACons;

end




@inline function phs2cnsFlux(APhys::Array{Float64,2},gamma::Float64)::Array{Float64,2}

N = size(APhys,1);
ACons = zeros(Float64,N,4);

for i = 1:N
	ACons[i,1] = APhys[i,1];
	ACons[i,2] = APhys[i,1]*APhys[i,2];
	ACons[i,3] = APhys[i,3]/(gamma-1.0) + 0.5*APhys[i,1]*( APhys[i,2]*APhys[i,2] );
	ACons[i,4] = APhys[i,4];
end 

return ACons;

end


	
@inline function cns2phs(ACons::Array{Float64,2},gamma::Float64)::Array{Float64,2}
N = size(ACons,1);
APhys = zeros(Float64,N,3);

for i = 1:N
	APhys[i,1] = ACons[i,1];
	APhys[i,2] = ACons[i,2]/ACons[i,1];
	APhys[i,3] = (gamma-1.0)*( ACons[i,3] - 0.5*( ACons[i,2]*ACons[i,2])/ACons[i,1] );
end #for

return APhys;

end


@inline function cns2phsFlux(ACons::Array{Float64,2},gamma::Float64)::Array{Float64,2}
N = size(ACons,1);
APhys = zeros(Float64,N,4);

for i = 1:N
	APhys[i,1] = ACons[i,1];
	APhys[i,2] = ACons[i,2]/ACons[i,1];
	APhys[i,3] = (gamma-1.0)*( ACons[i,3] - 0.5*( ACons[i,2]*ACons[i,2])/ACons[i,1] );
	APhys[i,4] =  ACons[i,4];
end #for

return APhys;

end


@inline function calcMagMaxVelocity(
	cSound::Array{Float64,1}, 
	magMaxU::Array{Float64,1}, 
	UphysCells::Array{Float64,2},
	Gamma::Float64
	)::Float64
	
	cSound= (Gamma*UphysCells[:,3]./UphysCells[:,1]).^(0.5);
	magMaxU = abs.( UphysCells[:,2] .+ cSound);
	
	VMax,id = findmax( magMaxU  );
	
	return VMax; 

end



