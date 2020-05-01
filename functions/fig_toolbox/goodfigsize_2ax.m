function figpos = goodfigsize_2ax(FIG)

% get the perfect figure size for large paneled figures (SEML,
% plot_campbell, etc...)

scrsiz = get(0,'screensize');
xsiz = .7*scrsiz(3);
ysiz = .7*scrsiz(4);
xpos = (scrsiz(3)-xsiz)/2;
ypos = (scrsiz(4)- ysiz)/1.5;
figpos = [xpos ypos xsiz ysiz];