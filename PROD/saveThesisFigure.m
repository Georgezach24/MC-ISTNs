function saveThesisFigure(fig, pngPath)
%SAVETHESISFIGURE PNG 300 dpi χωρίς περικοπή, και κλείσιμο.
set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
print(fig, pngPath, '-dpng', '-r300');
close(fig);
end
