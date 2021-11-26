using Polynomial, LinearAlgebra



function MFDFA(data::Array, scales::Array, q::Array)
	"""
    data ---------> signal or data array 
    scales -------> array of scales 
    q ------------> array of moments 
    """
    
    Profile = cumsum( data .- mean(data) ) ; 
    SP1_Ord = sqrt( mean(Profile.^2) ) ;
    
    X = Profile ;
    X = X' ;
    function polyfit(xVals,yVals)
        n = length(xVals)
        xBar, yBar = @fastmath mean(xVals), mean(yVals)
        sXX, sXY = @fastmath ones(n)'*(xVals.-xBar).^2 , dot(xVals.-xBar,yVals.-yBar)
        b1A = @fastmath sXY/sXX
        b0A = @fastmath yBar - b1A*xBar
        return b0A, b1A
    end
    """ Analisis Multifractal sin Tendencia  """
    
    segments = zeros(Int64, (1,length(scales))) ;
    global qRMS = zeros( length(q) ,length(scales)  ) ;
    global Fq = zeros( length(q) , length(scales) ) ;
    global RMScell = Array{Float64}[] ;
    global qRMScell =[] ;
    global segmentsFq = [] ;
    @inbounds for ns = 1:length(scales)
        global segments[ns] = Int(floor( length(X)/scales[ns] ) ) ;  
        global ft = zeros(Float64, (segments[ns], length(scales) ) ) ;
        global RMS = zeros(Float64, segments[ns]);
        @inbounds  for v=1:segments[ns]
            global Index = ( (v-1)*scales[ns] ) + 1: v*scales[ns] ;
            global C = polyfit( Index, X[Index]) ;
            global p = Polynomial(C) ;
            ft =p.(Index) ;
            RMS[v] = sqrt(mean((X[Index] .- ft).^2))  ;
            end
            l = deepcopy(RMS)
            push!(RMScell,l)
            global IndexFq = ((ns-1)*length(q) ) + 1 : ns*length(q) ;
            push!(segmentsFq, IndexFq) ;
        @inbounds for nq = 1:length(q)
            l = RMScell[ns].^q[nq]
            r = deepcopy(l) ;
            push!(qRMScell, r) ;
            end
        @inbounds for nq = 1: length(scales)
            Fq[nq,ns] = mean( qRMScell[segmentsFq[ns]][nq] ).^(1/q[nq] ) ;
            end
        Fq[findall(x->x==0, q)[1], ns] = exp( 0.5*mean(log.(RMScell[ns].^2) ) ) ;
    end
    
    Hq = zeros( Float64,length(q) ) ;
    global qRegLine = Array{Float64}[] ;
   @inbounds  for  nq = 1:length(q)
        global C = polyfit( log2.(scales),log2.(Fq[nq,:]) ) ;
        Hq[nq] = C[2] ;
        global p = Polynomial(C) ;
        push!( qRegLine, p.( log2.(scales) ) )
    end
    
    
    tq = Hq.*q .- 1 ;
    hq = diff(tq)./(q[2]-q[1]) ;
    
    Dq = ( q[1:end-1].*hq ) - tq[1:end-1] ;
    return Hq, tq, Dq, Fq, hq
    
    end