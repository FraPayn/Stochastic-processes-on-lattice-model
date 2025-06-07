function linsearch(v,C) # v is my random number, while c is the cumulative vector  
    j=0
    for i in C
        if i>v
            #print(j)
            return j-1
        end
        j+=1
    end
end

function binsearch1(v,C)
    #print("0")
    if(length(C)<=2)
        print("fine")
        A = C[1]
        return C[1]
    end
    #print(length(C))
    A=0
    if(length(C)%2!==0)
        print("zaaa")
        if(length(C)<=2)
            print("fine")
            A = C[1]
            return C[1]
        end
        if (v>C[Int(((length(C)-1)/2)+1)])
            print("b")
            binsearch(v,C[Int(((length(C)-1)/2)+2):Int(length(C))])
            
        else
            print("c")
            binsearch(v,C[1:Int((length(C)-1)/2)+1])
            
        end
    end
    #print(" ",length(C)," ",C[Int(((length(C))/2))]," ",v,"\n")
    if(length(C)%2==0)
        print("\n nu ",v,"   ",C,"\n")
        if(length(C)<=2)
            A=C[1]
            return C[1]
        end
        if (v>C[Int(((length(C))/2))])
            binsearch(v,C[Int(((length(C))/2)+1):Int(length(C))])
            print("d")
            
        end
        if(v<=C[Int(((length(C))/2))])
            binsearch(v,C[2:Int(((length(C))/2))])
            print("e")
        end
    end
    return A
end
# 3494334231 dumber

function sampler(weights::Vector, nsamples::Int)
    sum(weights) ≈ 1 || error("non normalized weigths")
    cumvec = vcat(0,cumsum(weights))
    #print([linsearch(rand(),cumvec) for i in 1:nsamples])
    return [linsearch(rand(),cumvec) for i in 1:nsamples]
end

function sampler_bin(weights::Vector, nsamples::Int)
    sum(weights) ≈ 1 || error("non normalized weigths")
    cumvec = vcat(0,cumsum(weights))
    return binsearch(0.223,cumvec)
    #return [binsearch(rand(),cumvec) for i in 1:nsamples]
end

function binw(n,p)
    W = [binomial(n,k)*(p^k)*((1-p)^(n-k)) for k in 0:n]
    return W
end
