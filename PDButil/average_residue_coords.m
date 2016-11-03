function[resout] = average_residue_coords(pdbin)

    allchains = unique(pdbin.chainID);
    resout = [];
    
    for(ch = allchains)
       this_chain = get_chain(pdbin,ch);
       allres = unique(this_chain.resNum);
       this_chain_avgres = zeros(length(allres),3);
       cc = 1;
       for( resi = allres ) 
           rr = get_residue(this_chain,resi);
           this_chain_avgres(cc,1) = mean(rr.X);
           this_chain_avgres(cc,2) = mean(rr.Y);
           this_chain_avgres(cc,3) = mean(rr.Z);
           cc = cc + 1;
       end
       resout = [resout; this_chain_avgres];
    end

end