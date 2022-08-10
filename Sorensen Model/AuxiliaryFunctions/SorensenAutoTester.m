%========================================================================================================================================
% SCRIPT SorensenAutoTester.m: automatic pure-Matlab (no C++) Tester
% Sorensen V01.01.41 20190724 (Gemini 13.01.06, BMLib 10.0.2, Autocoder 02.11.09, coded 24-Jul-2019 16:06:54)
%========================================================================================================================================


SorensenInitializeParvals;      %  calls DetermineParamCLeters and produces the initialization values for the variables (parameters initialized BEFORE variables)
SorensenInitializeStateVars;
TIME = Tzero:Tdelta:Tend;
nTimes = length(TIME);
STATEVARS = zeros(nTimes,nDepVars);
STATEVARS(1,:) = CurrentY;

% 20190709 Andrea De Gaetano, Marcello Pompa: Runge Kutta 4
% NOTE: this MATLAB implementation of RK4 uses vector assignments and therefore does NOT replicate exactly the corresponding C++ RK routine.
% By vector-assigning the increments to the (differentially-defined) Ys, the intermediate and final computation of the derivatives of one Y depend
% on the intermediate and final values of the other Ys. Conversely, when using loops, the intermediate and final computation of the derivatives of
% one Y depend on the interval INITIAL values of the other Ys, same as in C++. This last arrangement is consistent with the fact that the algebraic
% variables too are updated only after the computation of all the differential variables and that therefore the algebraic variable values used
% during all of the intermediate computations are the interval initial values.
% For smooth functions the differences in the implementations ought to be minimal. In case of differences, the C++ version is the logically
% consistent one. In any case, for non-smooth functions a variable-step integrator such as Runge-Kutta-Fehlberg ought to be used.

fprintf(1,'\n.');
Time = TIME(1);
for (tk=1:(nTimes-1))
   CurrentY0 = CurrentY;

   SorensenCurrentY2NamedVars;
   SorensenComputeDerivatives;
   dY1 = out1;

   Time = TIME(tk) + Tdelta/2;
   CurrentY = CurrentY0 + Tdelta/2 * dY1;
   SorensenCurrentY2NamedVars;
   SorensenComputeDerivatives;
   dY2 = out1;

   CurrentY = CurrentY0 + Tdelta/2 * dY2;
   SorensenCurrentY2NamedVars;
   SorensenComputeDerivatives;
   dY3 = out1;

   Time = TIME(tk+1);
   CurrentY = CurrentY0 + Tdelta * dY3;
   SorensenCurrentY2NamedVars;
   SorensenComputeDerivatives;
   dY4 = out1;

   CurrentY = CurrentY0 + Tdelta * (dY1/6 + dY2/3 + dY3/3 + dY4/6);
   SorensenForceVars;
   SorensenCurrentY2NamedVars; % needed to do the next computations
   SorensenComputeDiracs;
   SorensenNamedVars2CurrentY; % rebuild CurrentY
   SorensenForceVars;
   SorensenCurrentY2NamedVars; % needed to do the next computations
   SorensenComputeAlgebraic;
   SorensenNamedVars2CurrentY; % rebuild CurrentY
   SorensenForceVars;
   STATEVARS(tk+1,:) = CurrentY;
   fprintf(1,'.'); if (mod(tk,50)==0) fprintf(1,'\n'); end  % show progress
end

SorensenLoadNames;
stdxlabel = sprintf('Time [%s]',TimeUnit);
for (v=1:nDepVars)
   prompt = sprintf('Do you want to plot %s? (y/n) [y]: ',char(DepVarName(v)));
   str = lower(input(prompt,'s'));
   if isempty(str)
       str = 'y';
   end
   str = str(1);  % use only first character of input
   if (strcmp('y',str))
      newfig = figure;
      plot(TIME, STATEVARS(:,v))
      ylabel(char(DepVarLongName(v)));
      xlabel(stdxlabel);
      str = input('ENTER to continue (s to save)','s');
      if strcmp('s',str)
         figname = sprintf('Sorensen.Var%03d(%s).tif',v,char(DepVarName(v)));
         print (newfig,'-dtiff',figname);
      end
      close (newfig);
   end
end
%========================================================================================================================================