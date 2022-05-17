#spectralRadiusOfMappingFast2,# TODO: deprecated
#PhiLeft,PhiRight,PhiSteps# TODO: deprecated

#fixPointOfMappingLR, spectralRadiusOfMappingLR,  #TODO: renamed to ("LR" is removed)



struct DiscreteMapping_LR{tT,mxT,vT}
    ts::Vector{tT}
    LmappingMX::mxT
    RmappingMX::mxT
    mappingVs::Vector{vT}
end


###############################################################################
############################## Discrete Mapping ###############################
###############################################################################
function DiscreteMapping_LR(LDDEP::AbstractLDDEProblem, method::DiscretizationMethod, DiscretizationLength::Real; n_steps::Int64=nStepOfLength(DiscretizationLength, method.Δt), calculate_additive::Bool=false)
    result = calculateResults(LDDEP, method, DiscretizationLength,n_steps = n_steps, calculate_additive = calculate_additive);
    DiscreteMapping_LR(DiscreteMappingSteps_LR(result)...)
end



function DiscreteMappingSteps_LR(rst::AbstractResult{d}) where d
    NT=size(rst.subMXs[1])[1];
    Nτ=rst.n ÷ d;
    Nmaximal=max(Nτ,NT);

    ΦL=spdiagm( 0 => ones(d * Nmaximal));
    if Nτ>NT
        ΦR=spdiagm( -d * (NT) => ones(d * (Nτ-NT)));
    else
        ΦR=spzeros(d * Nmaximal, d * Nmaximal );
    end

   for smx_t in rst.subMXs
        for (it,smx) in enumerate(smx_t) 
           # pos=(it-1)*d
           # addSubmatrixToL!(ΦL,smx,pos,pos,NT*d)
           # addSubmatrixToR!(ΦR,smx,pos,pos-NT*d,Nmaximal*d)
            addSubmatrixToL!(ΦL,smx,(it-1)*d,it*d,NT*d)
            addSubmatrixToR!(ΦR,smx,(it-1)*d,(it-NT)*d,Nmaximal*d)
        end
    end
    
    mappingVs = [spzeros(size(ΦL,1))];

    for (it,subV) in enumerate(rst.subVs) 
       # pos=(it-1)*d
       # mappingVs[1][(1+pos) : (d+pos)] .+= subV.V
        mappingVs[1][(1+(it-1)*d) : (d+(it-1)*d)] .+= subV.V
    end
    ([rst.ts[1],rst.ts[end]],ΦL,ΦR,mappingVs)
end

function addSubmatrixToL!(Φ::AbstractArray, subMX::SubMX, it_shift1, it_shift2, nlimit)
    for i in eachindex(subMX.ranges, subMX.MXs)
        positionrange=LR_positionshift(subMX.ranges[i],it_shift1,it_shift2);
        if (positionrange[2][end])<=nlimit && (positionrange[2][1])>0 
            Φ[positionrange...] .+= -subMX.MXs[i];
        end
    end
end

function addSubmatrixToR!(Φ::AbstractArray, subMX::SubMX, it_shift1, it_shift2, nlimit)
    for i in eachindex(subMX.ranges, subMX.MXs)
        positionrange=LR_positionshift(subMX.ranges[i],it_shift1,it_shift2);
        if (positionrange[2][1])>0 
            Φ[positionrange...] .+= subMX.MXs[i];
        end
    end
end


function LR_positionshift(range::Tuple{UnitRange{Int64},UnitRange{Int64}},d)
    ((range[1][1]+d):(range[1][end]+d),(range[2][1]+d):(range[2][end]+d))
end

function LR_positionshift(range::Tuple{UnitRange{Int64},UnitRange{Int64}},d1,d2)
    ((range[1][1]+d1):(range[1][end]+d1),(range[2][1]+d2):(range[2][end]+d2))
end
      

function spectralRadiusOfMapping(mappLR::DiscreteMapping_LR; args...)
    return abs(eigs(mappLR.RmappingMX,mappLR.LmappingMX; args...)[1][1])
end

function fixPointOfMapping(mappLR::DiscreteMapping_LR)
    (mappLR.LmappingMX-mappLR.RmappingMX) \ Vector(mappLR.mappingVs[1])
end





############ TODO: deprecated



abstract type PhiLR{T} end




struct PhiSteps{Ti}
    dm::DiscreteMapping
    NT::Ti
    Nτ::Ti
    BS::Ti#BlockSize
end

function PhiSteps(dm::DiscreteMapping,BS::Integer)
    NT=size(dm.mappingMXs)[1];
    Nτ=size(dm.mappingMXs[1])[2] ÷ BS;
    PhiSteps(dm,NT,Nτ,BS);
 end
 Base.size(Φs::PhiSteps)  = (Φs.NT*Φs.BS,(Φs.Nτ+Φs.NT)*Φs.BS,);
 Base.size(Φs::PhiSteps, indDim::Integer) = indDim==1 ? Φs.NT*Φs.BS : (Φs.Nτ+Φs.NT)*Φs.BS;

 function Base.getindex(Φs::PhiSteps, idx1,idx2) 
    idxNT,subidx1=divrem(idx1-1, Φs.BS);
    idxNT += 1;
    subidx1 += 1;
    subidx2 = idx2 - (idxNT-1) * Φs.BS;
    if subidx2<=Φs.BS
        return subidx1==subidx2 ? 1.0 : 0.0
    elseif (subidx2-2)>Φs.Nτ * Φs.BS
        return 0.0
    else
        return -Φs.dm.mappingMXs[idxNT][subidx1,subidx2-Φs.BS]
    end
end


struct PhiLeft <:AbstractMatrix{Float64}
    properties
    Ps::PhiSteps
end

function Base.getindex(ΦL::PhiLeft, idx1,idx2)
    if idx1 <= (ΦL.Ps.BS * ΦL.Ps.NT) && idx2 <= (ΦL.Ps.BS * ΦL.Ps.NT)
        return ΦL.Ps[idx1,idx2]
    elseif idx1==idx2
        return 1.0
    else
        return 0.0
    end
    return 0.0
end
Base.size(ΦL::PhiLeft)  = (maximum([ΦL.Ps.Nτ,ΦL.Ps.NT])*ΦL.Ps.BS,maximum([ΦL.Ps.Nτ,ΦL.Ps.NT])*ΦL.Ps.BS,);
Base.size(ΦL::PhiLeft, Icucc::Integer) = maximum([ΦL.Ps.Nτ,ΦL.Ps.NT])*ΦL.Ps.BS;
Base.zero(ΦL::PhiLeft)= zeros(typeof(ΦL.Ps.dm.mappingMXs[1][1]),size(ΦL)...);


struct PhiRight <:AbstractMatrix{Float64}
    Ps::PhiSteps
end

function Base.getindex(ΦR::PhiRight, idx1,idx2)
    if idx1 <= (ΦR.Ps.BS * ΦR.Ps.NT)
        return -ΦR.Ps[idx1,idx2+ΦR.Ps.BS * ΦR.Ps.NT]
    elseif (idx1-ΦR.Ps.BS * ΦR.Ps.NT)==idx2
        return 1.0
    else
        return 0.0
    end
    return 0.0
end
Base.size(ΦR::PhiRight)  = (maximum([ΦR.Ps.Nτ,ΦR.Ps.NT])*ΦR.Ps.BS,maximum([ΦR.Ps.Nτ,ΦR.Ps.NT])*ΦR.Ps.BS,);
Base.size(ΦR::PhiRight, Icucc::Integer) = maximum([ΦR.Ps.Nτ,ΦR.Ps.NT])*ΦR.Ps.BS;
Base.zero(ΦR::PhiRight)= zeros(typeof(ΦR.Ps.dm.mappingMXs[1][1]),size(ΦR)...);



function spectralRadiusOfMappingFast2(dm::DiscreteMapping, BlockSize::Integer; args...)
  
    #@show  BlockSize=size(mathieu_lddep.A(0))[2];
    NT=size(dm.mappingMXs)[1];
    Nτ=size(dm.mappingMXs[1])[2] ÷ BlockSize;
    #if 0>1
    #    Φsteps = spzeros(BlockSize * NT, BlockSize * (NT + Nτ ));
    #    for k in 1:NT
    #    Φsteps[(1:BlockSize) .+ (k-1) * BlockSize,(1:BlockSize) .+ (k-1) * BlockSize] .= diagm(ones(BlockSize));
    #    Φsteps[(1:BlockSize) .+ (k-1) * BlockSize,(1:(Nτ * BlockSize)) .+ (k) * BlockSize] .= -dm.mappingMXs[k][1:BlockSize,:];
    #    end
    #else
        Φsteps=spdiagm(BlockSize * NT, BlockSize * (NT + Nτ ), 0 => ones(BlockSize * NT))
        #dropzeros!.(Φsteps)
        for k in 1:NT
            Φsteps[(1:BlockSize) .+ (k-1) * BlockSize,(1:(Nτ * BlockSize)) .+ (k) * BlockSize] .= -dm.mappingMXs[k][1:BlockSize,:];
            # (I,J,V)=findnz(dm.mappingMXs[k][1:BlockSize,:])
            # Φsteps += sparse(I, J .+ k*BlockSize, V,BlockSize * NT, BlockSize * (NT + Nτ ));
        end
        dropzeros!(Φsteps)
    #end

    Nmaximal=maximum([Nτ,NT]);
    ΦL=spdiagm( 0 => ones(BlockSize * Nmaximal));
    ΦL[1:BlockSize * NT, 1:BlockSize * NT ] .= Φsteps[1:BlockSize * NT, 1:BlockSize * NT ];
    dropzeros!(ΦL)
    ΦR=spzeros(BlockSize * Nmaximal, BlockSize * Nmaximal );
    ΦR[1:BlockSize * NT,1:BlockSize * Nτ]= -Φsteps[1:BlockSize * NT,BlockSize * NT + 1 : end];
    if Nτ>NT
        ΦR[BlockSize * NT+1:end, 1:BlockSize * (Nτ-NT)] .= diagm(ones(BlockSize * (Nτ-NT)));
    end  
    #return 0.0
    #return abs(eigs(ΦR,ΦL; args...)[1][1]),ΦR,ΦL,Φsteps
    return ΦR,ΦL#,Φsteps
end

