% This script is written and read by pdetool and should NOT be edited.
% There are two recommended alternatives:
 % 1) Export the required variables from pdetool and create a MATLAB script
 %    to perform operations on these.
 % 2) Define the problem completely using a MATLAB script. See
 %    http://www.mathworks.com/help/pde/examples/index.html for examples
 %    of this approach.
function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',1);
set(ax,'DataAspectRatio',[1 1 1]);
set(ax,'PlotBoxAspectRatio',[1.5 1 1]);
set(ax,'XLim',[-1.5 1.5]);
set(ax,'YLim',[-1 1]);
set(ax,'XTickMode','auto');
set(ax,'YTickMode','auto');

% Geometry description:
pdepoly([ 1.0533596837944663,...
 1.0217391304347831,...
 0.96245059288537593,...
 0.92292490118577142,...
 0.85177865612648285,...
 0.80434782608695699,...
 0.74110671936758932,...
 0.67391304347826164,...
 0.63833992094861713,...
 0.54743083003952631,...
 0.48023715415019819,...
 0.43675889328063278,...
 0.40118577075098871,...
 0.31027667984189788,...
 0.26284584980237202,...
 0.21936758893280661,...
 0.17193675889328119,...
 0.12450592885375533,...
 0.081027667984190366,...
 0.033596837944664504,...
 -0.053359683794465873,...
 -0.14031620553359625,...
 -0.19960474308300347,...
 -0.29841897233201564,...
 -0.34980237154150151,...
 -0.42094861660079008,...
 -0.44071146245059278,...
 -0.46047430830039504,...
 -0.47628458498023685,...
 -0.4723320158102764,...
 -0.43280632411067144,...
 -0.41304347826086918,...
 -0.40909090909090873,...
 -0.41699604743082963,...
 -0.45652173913043459,...
 -0.49999999999999956,...
 -0.54743083003952531,...
 -0.58695652173913015,...
 -0.62648221343873489,...
 -0.66600790513833963,...
 -0.71739130434782572,...
 -0.75691699604743057,...
 -0.85968379446640297,...
 -0.91501976284584952,...
 -0.93083003952569154,...
 -0.99407114624505899,...
 -1.0494071146245059,...
 -1.1086956521739126,...
 -1.1719367588932803,...
 -1.2193675889328062,...
 -1.2391304347826084,...
 -1.2549407114624505,...
 -1.2391304347826084,...
 -1.2193675889328062,...
 -1.1640316205533594,...
 -1.0928853754940708,...
 -1.00197628458498,...
 -0.93873517786561245,...
 -0.91897233201580997,...
 -0.89920948616600771,...
 -0.89130434782608681,...
 -0.89920948616600771,...
 -0.92687747035573109,...
 -0.97826086956521718,...
 -1.0375494071146245,...
 -1.0889328063241104,...
 -1.1324110671936758,...
 -1.1679841897233199,...
 -1.2075098814229248,...
 -1.2312252964426875,...
 -1.2391304347826084,...
 -1.2233201581027666,...
 -1.1877470355731223,...
 -1.160079051383399,...
 -1.0691699604743081,...
 -0.96245059288537538,...
 -0.7885375494071144,...
 -0.65810276679841873,...
 -0.50395256916996001,...
 -0.42885375494071099,...
 -0.310276679841897,...
 -0.20750988142292437,...
 -0.0098814229249009067,...
 0.11660079051383443,...
 0.27865612648221383,...
 0.3498023715415024,...
 0.41699604743083052,...
 0.48814229249011909,...
 0.55928853754940766,...
 0.61067193675889397,...
 0.67786561264822209,...
 0.70553359683794525,...
 0.72924901185770796,...
 0.78853754940711518,...
 0.8715415019762851,...
 0.94664031620553413,...
 1.0533596837944672,...
 1.1284584980237162,...
 1.1798418972332021,...
 1.2035573122529653,...
 1.2312252964426884,...
 1.2747035573122538,...
 1.3418972332015815,...
 1.3814229249011865,...
 1.4090909090909096,...
 1.4051383399209492,...
 1.3735177865612656,...
 1.3181818181818188,...
 1.2588932806324116,...
 1.1719367588932812,...
 1.0928853754940717,...
 0.98221343873517863,...
 0.91501976284585052,...
 0.87944664031620601,...
 0.85177865612648285,...
 0.87549407114624556,...
 0.93083003952569232,...
 0.99011857707509954,...
 1.0533596837944672,...
 1.1324110671936767,...
 1.2075098814229257,...
 1.227272727272728,...
 1.2351778656126489,...
 1.2233201581027675,...
 1.1758893280632416,...
 1.144268774703558,...
 1.1324110671936767,...
 1.112648221343874,...
 1.0849802371541508,...
],...
[ 0.52766798418972316,...
 0.59090909090909127,...
 0.64229249011857736,...
 0.67786561264822165,...
 0.71739130434782639,...
 0.73320158102766819,...
 0.75296442687747067,...
 0.75691699604743112,...
 0.75691699604743112,...
 0.72529644268774729,...
 0.689723320158103,...
 0.63438735177865624,...
 0.59486166007905172,...
 0.58695652173913082,...
 0.61067193675889353,...
 0.67786561264822165,...
 0.72529644268774729,...
 0.78853754940711474,...
 0.83992094861660105,...
 0.86758893280632443,...
 0.88735177865612691,...
 0.88735177865612691,...
 0.88735177865612691,...
 0.87154150197628488,...
 0.85968379446640331,...
 0.8438735177865615,...
 0.83201581027668015,...
 0.79249011857707541,...
 0.72924901185770774,...
 0.68577075098814255,...
 0.63833992094861669,...
 0.59881422924901218,...
 0.55533596837944676,...
 0.5276679841897236,...
 0.5118577075098818,...
 0.5118577075098818,...
 0.54347826086956541,...
 0.57114624505928879,...
 0.60276679841897263,...
 0.63438735177865624,...
 0.67786561264822165,...
 0.69762845849802391,...
 0.71343873517786571,...
 0.72924901185770774,...
 0.77667984189723338,...
 0.81620553359683834,...
 0.82015810276679879,...
 0.80830039525691721,...
 0.77667984189723338,...
 0.71343873517786571,...
 0.65415019762845872,...
 0.50000000000000022,...
 0.42094861660079075,...
 0.3300395256916997,...
 0.21936758893280661,...
 0.1679841897233203,...
 0.11264822134387376,...
 0.037549407114624733,...
 -0.0019762845849800037,...
 -0.06126482213438722,...
 -0.13636363636363624,...
 -0.19565217391304346,...
 -0.2351778656126482,...
 -0.26679841897233181,...
 -0.31027667984189722,...
 -0.35375494071146241,...
 -0.40909090909090895,...
 -0.45256916996047414,...
 -0.51581027667984181,...
 -0.58300395256916993,...
 -0.65019762845849804,...
 -0.70158102766798436,...
 -0.73320158102766797,...
 -0.74901185770750978,...
 -0.76877470355731203,...
 -0.79249011857707519,...
 -0.80434782608695654,...
 -0.82015810276679835,...
 -0.88339920948616601,...
 -0.90711462450592872,...
 -0.92687747035573143,...
 -0.92292490118577053,...
 -0.88339920948616601,...
 -0.82806324110671925,...
 -0.67391304347826075,...
 -0.57114624505928857,...
 -0.51581027667984181,...
 -0.51976284584980226,...
 -0.52371541501976271,...
 -0.5671936758893279,...
 -0.61067193675889309,...
 -0.67391304347826075,...
 -0.72529644268774707,...
 -0.78063241106719383,...
 -0.81225296442687744,...
 -0.8162055335968379,...
 -0.80039525691699609,...
 -0.74505928853754932,...
 -0.6620553359683794,...
 -0.56324110671936745,...
 -0.43675889328063233,...
 -0.31818181818181812,...
 -0.22332015810276662,...
 -0.16798418972332008,...
 -0.073122529644268575,...
 -0.01383399209486158,...
 0.029644268774703608,...
 0.061264822134387442,...
 0.084980237154150373,...
 0.10079051383399218,...
 0.11660079051383421,...
 0.14031620553359692,...
 0.20750988142292504,...
 0.24308300395256932,...
 0.32213438735177879,...
 0.38142292490118601,...
 0.42094861660079075,...
 0.43675889328063278,...
 0.40909090909090917,...
 0.40513833992094872,...
 0.42094861660079075,...
 0.44466403162055368,...
 0.49604743083003977,...
 0.53557312252964451,...
 0.55138339920948631,...
 0.53952569169960496,...
 0.51581027667984225,...
 0.50000000000000022,...
 0.50395256916996067,...
],...
 'P1');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','P1')

% Boundary conditions:
pdetool('changemode',0)
pdesetbd(129,...
'dir',...
1,...
'1',...
'0')
pdesetbd(128,...
'dir',...
1,...
'1',...
'0')
pdesetbd(127,...
'dir',...
1,...
'1',...
'0')
pdesetbd(126,...
'dir',...
1,...
'1',...
'0')
pdesetbd(125,...
'dir',...
1,...
'1',...
'0')
pdesetbd(124,...
'dir',...
1,...
'1',...
'0')
pdesetbd(123,...
'dir',...
1,...
'1',...
'0')
pdesetbd(122,...
'dir',...
1,...
'1',...
'0')
pdesetbd(121,...
'dir',...
1,...
'1',...
'0')
pdesetbd(120,...
'dir',...
1,...
'1',...
'0')
pdesetbd(119,...
'dir',...
1,...
'1',...
'0')
pdesetbd(118,...
'dir',...
1,...
'1',...
'0')
pdesetbd(117,...
'dir',...
1,...
'1',...
'0')
pdesetbd(116,...
'dir',...
1,...
'1',...
'0')
pdesetbd(115,...
'dir',...
1,...
'1',...
'0')
pdesetbd(114,...
'dir',...
1,...
'1',...
'0')
pdesetbd(113,...
'dir',...
1,...
'1',...
'0')
pdesetbd(112,...
'dir',...
1,...
'1',...
'0')
pdesetbd(111,...
'dir',...
1,...
'1',...
'0')
pdesetbd(110,...
'dir',...
1,...
'1',...
'0')
pdesetbd(109,...
'dir',...
1,...
'1',...
'0')
pdesetbd(108,...
'dir',...
1,...
'1',...
'0')
pdesetbd(107,...
'dir',...
1,...
'1',...
'0')
pdesetbd(106,...
'dir',...
1,...
'1',...
'0')
pdesetbd(105,...
'dir',...
1,...
'1',...
'0')
pdesetbd(104,...
'dir',...
1,...
'1',...
'0')
pdesetbd(103,...
'dir',...
1,...
'1',...
'0')
pdesetbd(102,...
'dir',...
1,...
'1',...
'0')
pdesetbd(101,...
'dir',...
1,...
'1',...
'0')
pdesetbd(100,...
'dir',...
1,...
'1',...
'0')
pdesetbd(99,...
'dir',...
1,...
'1',...
'0')
pdesetbd(98,...
'dir',...
1,...
'1',...
'0')
pdesetbd(97,...
'dir',...
1,...
'1',...
'0')
pdesetbd(96,...
'dir',...
1,...
'1',...
'0')
pdesetbd(95,...
'dir',...
1,...
'1',...
'0')
pdesetbd(94,...
'dir',...
1,...
'1',...
'0')
pdesetbd(93,...
'dir',...
1,...
'1',...
'0')
pdesetbd(92,...
'dir',...
1,...
'1',...
'0')
pdesetbd(91,...
'dir',...
1,...
'1',...
'0')
pdesetbd(90,...
'dir',...
1,...
'1',...
'0')
pdesetbd(89,...
'dir',...
1,...
'1',...
'0')
pdesetbd(88,...
'dir',...
1,...
'1',...
'0')
pdesetbd(87,...
'dir',...
1,...
'1',...
'0')
pdesetbd(86,...
'dir',...
1,...
'1',...
'0')
pdesetbd(85,...
'dir',...
1,...
'1',...
'0')
pdesetbd(84,...
'dir',...
1,...
'1',...
'0')
pdesetbd(83,...
'dir',...
1,...
'1',...
'0')
pdesetbd(82,...
'dir',...
1,...
'1',...
'0')
pdesetbd(81,...
'dir',...
1,...
'1',...
'0')
pdesetbd(80,...
'dir',...
1,...
'1',...
'0')
pdesetbd(79,...
'dir',...
1,...
'1',...
'0')
pdesetbd(78,...
'dir',...
1,...
'1',...
'0')
pdesetbd(77,...
'dir',...
1,...
'1',...
'0')
pdesetbd(76,...
'dir',...
1,...
'1',...
'0')
pdesetbd(75,...
'dir',...
1,...
'1',...
'0')
pdesetbd(74,...
'dir',...
1,...
'1',...
'0')
pdesetbd(73,...
'dir',...
1,...
'1',...
'0')
pdesetbd(72,...
'dir',...
1,...
'1',...
'0')
pdesetbd(71,...
'dir',...
1,...
'1',...
'0')
pdesetbd(70,...
'dir',...
1,...
'1',...
'0')
pdesetbd(69,...
'dir',...
1,...
'1',...
'0')
pdesetbd(68,...
'dir',...
1,...
'1',...
'0')
pdesetbd(67,...
'dir',...
1,...
'1',...
'0')
pdesetbd(66,...
'dir',...
1,...
'1',...
'0')
pdesetbd(65,...
'dir',...
1,...
'1',...
'0')
pdesetbd(64,...
'dir',...
1,...
'1',...
'0')
pdesetbd(63,...
'dir',...
1,...
'1',...
'0')
pdesetbd(62,...
'dir',...
1,...
'1',...
'0')
pdesetbd(61,...
'dir',...
1,...
'1',...
'0')
pdesetbd(60,...
'dir',...
1,...
'1',...
'0')
pdesetbd(59,...
'dir',...
1,...
'1',...
'0')
pdesetbd(58,...
'dir',...
1,...
'1',...
'0')
pdesetbd(57,...
'dir',...
1,...
'1',...
'0')
pdesetbd(56,...
'dir',...
1,...
'1',...
'0')
pdesetbd(55,...
'dir',...
1,...
'1',...
'0')
pdesetbd(54,...
'dir',...
1,...
'1',...
'0')
pdesetbd(53,...
'dir',...
1,...
'1',...
'0')
pdesetbd(52,...
'dir',...
1,...
'1',...
'0')
pdesetbd(51,...
'dir',...
1,...
'1',...
'0')
pdesetbd(50,...
'dir',...
1,...
'1',...
'0')
pdesetbd(49,...
'dir',...
1,...
'1',...
'0')
pdesetbd(48,...
'dir',...
1,...
'1',...
'0')
pdesetbd(47,...
'dir',...
1,...
'1',...
'0')
pdesetbd(46,...
'dir',...
1,...
'1',...
'0')
pdesetbd(45,...
'dir',...
1,...
'1',...
'0')
pdesetbd(44,...
'dir',...
1,...
'1',...
'0')
pdesetbd(43,...
'dir',...
1,...
'1',...
'0')
pdesetbd(42,...
'dir',...
1,...
'1',...
'0')
pdesetbd(41,...
'dir',...
1,...
'1',...
'0')
pdesetbd(40,...
'dir',...
1,...
'1',...
'0')
pdesetbd(39,...
'dir',...
1,...
'1',...
'0')
pdesetbd(38,...
'dir',...
1,...
'1',...
'0')
pdesetbd(37,...
'dir',...
1,...
'1',...
'0')
pdesetbd(36,...
'dir',...
1,...
'1',...
'0')
pdesetbd(35,...
'dir',...
1,...
'1',...
'0')
pdesetbd(34,...
'dir',...
1,...
'1',...
'0')
pdesetbd(33,...
'dir',...
1,...
'1',...
'0')
pdesetbd(32,...
'dir',...
1,...
'1',...
'0')
pdesetbd(31,...
'dir',...
1,...
'1',...
'0')
pdesetbd(30,...
'dir',...
1,...
'1',...
'0')
pdesetbd(29,...
'dir',...
1,...
'1',...
'0')
pdesetbd(28,...
'dir',...
1,...
'1',...
'0')
pdesetbd(27,...
'dir',...
1,...
'1',...
'0')
pdesetbd(26,...
'dir',...
1,...
'1',...
'0')
pdesetbd(25,...
'dir',...
1,...
'1',...
'0')
pdesetbd(24,...
'dir',...
1,...
'1',...
'0')
pdesetbd(23,...
'dir',...
1,...
'1',...
'0')
pdesetbd(22,...
'dir',...
1,...
'1',...
'0')
pdesetbd(21,...
'dir',...
1,...
'1',...
'0')
pdesetbd(20,...
'dir',...
1,...
'1',...
'0')
pdesetbd(19,...
'dir',...
1,...
'1',...
'0')
pdesetbd(18,...
'dir',...
1,...
'1',...
'0')
pdesetbd(17,...
'dir',...
1,...
'1',...
'0')
pdesetbd(16,...
'dir',...
1,...
'1',...
'0')
pdesetbd(15,...
'dir',...
1,...
'1',...
'0')
pdesetbd(14,...
'dir',...
1,...
'1',...
'0')
pdesetbd(13,...
'dir',...
1,...
'1',...
'0')
pdesetbd(12,...
'dir',...
1,...
'1',...
'0')
pdesetbd(11,...
'dir',...
1,...
'1',...
'0')
pdesetbd(10,...
'dir',...
1,...
'1',...
'0')
pdesetbd(9,...
'dir',...
1,...
'1',...
'0')
pdesetbd(8,...
'dir',...
1,...
'1',...
'0')
pdesetbd(7,...
'dir',...
1,...
'1',...
'0')
pdesetbd(6,...
'dir',...
1,...
'1',...
'0')
pdesetbd(5,...
'dir',...
1,...
'1',...
'0')
pdesetbd(4,...
'dir',...
1,...
'1',...
'0')
pdesetbd(3,...
'dir',...
1,...
'1',...
'0')
pdesetbd(2,...
'dir',...
1,...
'1',...
'0')
pdesetbd(1,...
'dir',...
1,...
'1',...
'0')

% Mesh generation:
setappdata(pde_fig,'Hgrad',1.3);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
setappdata(pde_fig,'MesherVersion','preR2013a');
pdetool('initmesh')
pdetool('refine')
pdetool('refine')
pdetool('refine')

% PDE coefficients:
pdeseteq(1,...
'1.0',...
'0.0',...
'10.0',...
'1.0',...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['1.0 ';...
'0.0 ';...
'10.0';...
'1.0 '])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
char('0','117024','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[1 1 1 1 1 1 1 1 0 0 0 1 1 0 0 0 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');