function lagrange_polynomial(x, vecx, vecf)
    l = length(vecx)
    count=0
    for i in 1:l
        count2 = 1
        for j in 1:l
            if j==i
                continue 
            end
            count2 *= ((x-vecx[j])/(vecx[i]-vecx[j]))
        end
        count += vecf[i]*count2
    end
    return count
end

# the idea here is to pass from 1,2 and 3 point for every little ck - ck+1 interval

midpoint(a,b,f::Function) = f((a+b)/2)*(b-a)
trapezoid(a,b,f::Function) = ((f(a)+f(b))/2)*(b-a)
cavalieri(a,b,f::Function) = ((f(a)+f(b)+4*f((a+b)/2))/6)*(b-a)

function numint(f::Function, a::Number,b::Number; nstep::Integer=100, nlagrange::Integer)
    I = 0
    #print(a,"\n",b,"\n",nstep,"\n")
    base = LinRange(a,b,nstep)
   
    if nlagrange == 1
        for j in 2:length(base)
            I += midpoint(Number(base[j-1]),base[j],f)
        end
    end
    if nlagrange == 2
        for j in 2:length(base)
            I += trapezoid(base[j-1],base[j],f)
        end
    end
    if nlagrange == 3
        for j in 2:length(base)
            I += cavalieri(base[j-1],base[j],f)
        end
    end
    return I


end
