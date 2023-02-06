

function monod2=f_monod2(PO4,u_max,kN)

monod2.u = (PO4./(PO4+kN))*u_max;

    scatter(PO4,monod2.u,5)
    
return