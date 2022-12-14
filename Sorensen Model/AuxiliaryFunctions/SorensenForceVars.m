%========================================================================================================================================
% SCRIPT SorensenForceVars.m: force the CurrentY state variables values within their prescribed limits
% Sorensen V01.01.41 20190724 (Gemini 13.01.06, BMLib 10.0.2, Autocoder 02.11.09, coded 24-Jul-2019 16:06:54)
%========================================================================================================================================


for (v = 1:nDepVars)
   if (CurrentY(v) < VARMIN(v+1))
      CurrentY(v) = VARMIN(v+1);
   elseif (CurrentY(v) > VARMAX(v+1))
      CurrentY(v) = VARMAX(v+1);
   end
end

%========================================================================================================================================