
function monod=f_monod(Nnew,Nreg,u_max,kN,n_fig)

monod.u = ((Nnew+Nreg)./(Nnew+Nreg+kN))*u_max;

if length(monod.u)>1
    if n_fig==1
    scatter(Nnew,monod.u,5)
    elseif n_fig==2
        scatter(Nnew+Nreg,monod.u,5)  
    end
end
    
return