function saveThesisFigure(fig, pngPath)
%SAVETHESISFIGURE Εξαγωγή PNG 300 dpi στο μέγεθος του σχήματος (χωρίς αυτόματο
% περικόψιμο, που έκοβε τίτλους) και κλείσιμο.
set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
print(fig, pngPath, '-dpng', '-r300');
close(fig);
end
