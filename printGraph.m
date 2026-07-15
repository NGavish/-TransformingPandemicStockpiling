%% Save figure
%% -------------------
% This function saves figure data into .fig and .eps files.  
% In addition it generates a latex statment for the figure.
%
% Input: fileName - determines the name of saved files
% Output:
%       fileName.fig    - Matlab figure file
%       fileName.eps  - eps file
%       fileName.tif  -    tiff file
%       fileName.txt     - Latex statment
%%
function printGraph(fileName)
% Save .fig file
savefig([fileName,'.fig'])
% Save .eps file
print([fileName,'.eps'],'-depsc');
% Save .eps file
print([fileName,'.pdf'],'-dpdf','-bestfit');
% Save .tif file
if ispc
print([fileName,'.tiff'],'-bestfit','-dtiff');
end
% Generate latex text
latexFigureStatment='\\begin{figure}[ht!]\n\\begin{center}\n\\scalebox{1}{\\includegraphics{%s.eps}}\n\\caption{}\n\\label{fig:%s}\n\\end{center}\n%%%s\n\\end{figure}';
currDir=cd;
fid=fopen([fileName,'.txt'],'w');
fprintf(fid,latexFigureStatment,fileName,fileName,currDir)
fclose(fid);
return