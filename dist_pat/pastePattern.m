function [realization, frozenRealiz] = pastePattern(Pattern, wx, wy, wz, realization, frozenRealiz, par, outT)

% This function will paste the pattern on the realization knowing the
% previous frozen nodes and and also updates the inner patch to become 
% frozen.

Pat = par.Pat;
Patz= par.Patz;
innerPatch = par.innerPatch;
innerPatchz= par.innerPatchz;

% paste the dual template first. 
% ---------------------------------------------
if par.bUseDualTemplate
    % only for the first multipleGrid.
    if par.m1 == 2^(par.multipleGrid-1)
        if par.bPasteDualOnFrozen
            realization(wx(1):wx(end),wy(1):wy(end),wz(1):wz(end)) = outT;
        else
            frozenTemplate      = reshape(frozenRealiz(wx(1):wx(end),wy(1):wy(end),wz(1):wz(end)),1,[]);
            realizationTemplate = reshape(realization (wx(1):wx(end),wy(1):wy(end),wz(1):wz(end)),1,[]);
            outT = reshape(outT, 1, []);
            outT(1,frozenTemplate == 1) = realizationTemplate(1,frozenTemplate == 1);
            realization(wx(1):wx(end),wy(1):wy(end),wz(1):wz(end)) = ...
                reshape(outT, wx(end)-wx(1)+1, wy(end)-wy(1)+1, wz(end)-wz(1)+1);
        end
    end
end

frozenTemplate      = reshape(frozenRealiz(wx, wy, wz),1,Pat^2*Patz);
realizationTemplate = reshape(realization (wx, wy, wz),1,Pat^2*Patz);


% change the pattern values for frozen nodes
%---------------------------------------------------
Pattern(1,frozenTemplate == 1) = realizationTemplate(1,frozenTemplate == 1);



% Paste the pattern
%---------------------
% (one technique)
realization(wx, wy, wz) = reshape(Pattern, Pat, Pat, Patz);
% (second technique)
% for i = 1:Pat
%     for j = 1:Pat
%         for k = 1:Patz
%           realization(wx(i),wy(j),wz(k)) = Pattern(1, i+(j-1)*Pat+(k-1)*Pat^2);
%         end
%     end
% end



% Update the frozen nodes matrix
%------------------------------------
middleIdx  = (Pat +1)/2; 
middleIdxz = (Patz+1)/2;
boundary   = (innerPatch -1)/2;
boundaryz  = (innerPatchz-1)/2;
wxInner = wx(middleIdx- boundary  : middleIdx +boundary);
wyInner = wy(middleIdx- boundary  : middleIdx +boundary);
wzInner = wz(middleIdxz-boundaryz : middleIdxz+boundaryz);
frozenRealiz(wxInner, wyInner, wzInner) = 1;




end
